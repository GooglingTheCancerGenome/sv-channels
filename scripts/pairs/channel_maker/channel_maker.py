# Imports

import argparse
import errno
import gzip
import logging
import os
import pickle
import statistics
from collections import defaultdict
from time import time

import bz2file
import numpy as np
import pyBigWig
import pysam
#from functions import get_one_hot_sequence_by_list, load_clipped_read_positions
from functions import *
from candidate_pairs import *

# import matplotlib.pyplot as plt

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = False

# Only clipped read positions supported by at least min_cr_support clipped reads are considered
min_cr_support = 3
# Window half length
win_hlen = 100
# Window size
win_len = win_hlen * 2


def get_chr_len(ibam, chrName):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    return chrLen


def create_dir(directory):
    '''
    Create a directory if it does not exist. Raises an exception if the directory exists.
    :param directory: directory to create
    :return: None
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def count_clipped_read_positions(cpos_cnt):
    '''

    :param cpos_cnt: dictionary of clipped read positions (keys) and counts of clipped reads per position (values) as
    returned by the clipped_read_pos.py script
    :return: None. Prints the number of clipped read positions with clipped read support greater than the integers
    specified in the range
    '''
    for i in range(0, 5):
        logging.info('Number of positions with at least %d clipped reads: %d' %
                     (i + 1, len([k for k, v in cpos_cnt.items() if v > i])))


def get_mappability_bigwig():
    mappability_file = "/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Mappability/GRCh37.151mer.bw" if HPC_MODE \
        else "/Users/lsantuari/Documents/Data/GEM/GRCh37.151mer.bw"
    bw = pyBigWig.open(mappability_file)

    return bw


def load_bam(ibam):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    return pysam.AlignmentFile(ibam, "rb")


def get_chr_len_dict(ibam):
    bamfile = load_bam(ibam)
    # Extract chromosome length from the BAM header
    header_dict = bamfile.header

    chrLen = {i['SN']: i['LN'] for i in header_dict['SQ']}
    return chrLen


def load_channels(sample, chr_list):

    prefix = ''
    channel_names = ['candidate_pairs', 'clipped_reads', 'clipped_read_distance',
                     'coverage', 'split_read_distance']

    channel_data = defaultdict(dict)

    for chrom in chr_list:
        logging.info('Loading data for Chr%s' % chrom)
        for ch in channel_names:
            logging.info('Loading data for channel %s' % ch)
            suffix = '.npy.bz2' if ch == 'coverage' else '.pbz2'
            if HPC_MODE:
                filename = os.path.join(prefix, sample, ch, '_'.join([chrom, ch + suffix]))
            else:
                filename = ch + suffix
            assert os.path.isfile(filename)

            logging.info('Reading %s for Chr%s' % (ch, chrom))
            with bz2file.BZ2File(filename, 'rb') as f:
                if suffix == '.npy.bz2':
                    channel_data[chrom][ch] = np.load(f)
                else:
                    channel_data[chrom][ch] = pickle.load(f)
            logging.info('End of reading')

        # unpack clipped_reads
        channel_data[chrom]['read_quality'], channel_data[chrom]['clipped_reads'], \
        channel_data[chrom]['clipped_reads_inversion'], channel_data[chrom]['clipped_reads_duplication'], \
        channel_data[chrom]['clipped_reads_translocation'] = channel_data[chrom]['clipped_reads']

        # unpack split_reads
        channel_data[chrom]['split_read_distance'], \
        channel_data[chrom]['split_reads'] = channel_data[chrom]['split_read_distance']

    return channel_data


def channel_maker(chrom, sampleName, outFile):

    n_channels = 29
    bp_padding = 10

    channel_data = load_channels(sampleName, [chrom])

    bw_map = get_mappability_bigwig()

    candidate_pairs_chr = [sv for sv in channel_data[chrom]['candidate_pairs']
                           if sv.tuple[0].chr == sv.tuple[1].chr and sv.tuple[0].chr == chrom]

    channel_windows = np.zeros(shape=(len(candidate_pairs_chr),
                                      win_len * 2 + bp_padding, n_channels), dtype=np.uint32)

    # dictionary of key choices
    direction_list = {'clipped_reads': ['left', 'right', 'D_left', 'D_right', 'I'],
                      'split_reads': ['left', 'right'],
                      'split_read_distance': ['left', 'right'],
                      'clipped_reads_inversion': ['before', 'after'],
                      'clipped_reads_duplication': ['before', 'after'],
                      'clipped_reads_translocation': ['opposite', 'same'],
                      'clipped_read_distance': ['forward', 'reverse']
                      }

    # Consider a single sample
    # sample_list = sampleName.split('_')

    # for sample in sample_list:

    positions = []
    for sv in candidate_pairs_chr:
        bp1, bp2 = sv.tuple
        positions.extend(list(range(bp1.pos - win_hlen, bp1.pos + win_hlen)) +
                         list(range(bp2.pos - win_hlen, bp2.pos + win_hlen)))
    positions = np.array(positions)

    idx = np.arange(win_len)
    idx2 = np.arange(start=win_len + bp_padding, stop=win_len * 2 + bp_padding)
    idx = np.concatenate((idx, idx2), axis=0)

    channel_index = 0

    for current_channel in ['coverage', 'read_quality',
                            'clipped_reads', 'split_reads',
                            'clipped_reads_inversion', 'clipped_reads_duplication',
                            'clipped_reads_translocation',
                            'clipped_read_distance', 'split_read_distance']:

        logging.info("Adding channel %s" % current_channel)

        if current_channel == 'coverage' or current_channel == 'read_quality':

            payload = channel_data[chrom][current_channel][positions]
            payload.shape = channel_windows[:, idx, channel_index].shape
            channel_windows[:, idx, channel_index] = payload
            channel_index += 1

        elif current_channel in ['clipped_reads', 'split_reads',
                               'clipped_reads_inversion', 'clipped_reads_duplication',
                               'clipped_reads_translocation']:
            for split_direction in direction_list[current_channel]:

                channel_pos = set(positions) & set(channel_data[chrom][current_channel][split_direction].keys())
                payload = [ channel_data[chrom][current_channel][split_direction][pos] if pos in channel_pos else 0 \
                 for pos in positions ]
                payload = np.array(payload)
                payload.shape = channel_windows[:, idx, channel_index].shape
                channel_windows[:, idx, channel_index] = payload
                channel_index += 1

        elif current_channel == 'clipped_read_distance':
            for split_direction in direction_list[current_channel]:
                for clipped_arrangement in ['left', 'right', 'all']:

                    channel_pos = set(positions) & \
                                  set(channel_data[chrom][current_channel][split_direction][clipped_arrangement].keys())
                    payload = [ statistics.median(
                        channel_data[chrom][current_channel][split_direction][clipped_arrangement][pos]) \
                                    if pos in channel_pos else 0 for pos in positions ]
                    payload = np.array(payload)
                    payload.shape = channel_windows[:, idx, channel_index].shape
                    channel_windows[:, idx, channel_index] = payload
                    channel_index += 1

        elif current_channel == 'split_read_distance':
            for split_direction in direction_list[current_channel]:

                channel_pos = set(positions) & \
                              set(channel_data[chrom][current_channel][split_direction].keys())
                payload = [ statistics.median(
                    channel_data[chrom][current_channel][split_direction][pos]) \
                                if pos in channel_pos else 0 for pos in positions ]
                payload = np.array(payload)
                payload.shape = channel_windows[:, idx, channel_index].shape
                channel_windows[:, idx, channel_index] = payload
                channel_index += 1

    current_channel = 'one_hot_encoding'
    logging.info("Adding channel %s" % current_channel)

    nuc_list = ['A', 'T', 'C', 'G', 'N']

    payload = get_one_hot_sequence_by_list(chrom, positions, HPC_MODE)
    payload.shape = channel_windows[:, idx, channel_index:channel_index+len(nuc_list)].shape
    channel_windows[:, idx, channel_index:channel_index+len(nuc_list)] = payload
    channel_index += len(nuc_list)

    current_channel = 'mappability'
    logging.info("Adding channel %s" % current_channel)

    payload = []
    for sv in candidate_pairs_chr:
        bp1, bp2 = sv.tuple
        payload.extend(bw_map.values(chrom, bp1.pos - win_hlen, bp1.pos + win_hlen) +
                         bw_map.values(chrom, bp2.pos - win_hlen, bp2.pos + win_hlen))
    payload = np.array(payload)
    payload.shape = channel_windows[:, idx, channel_index].shape
    channel_windows[:, idx, channel_index] = payload

    logging.info("channel_windows shape: %s" % str(channel_windows.shape))

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=channel_windows)
    f.close()


def inspect_windows(outFile):

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "r") as f:
        channel_windows = np.load(f)
    f.close()

    for i in range(29):
        print(channel_windows[0,:,i])


def main():
    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/'+
    #   'artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output")
    parser.add_argument('-s', '--sample', type=str, default='NA12878',
                        help="Specify sample")
    parser.add_argument('-l', '--logfile', default='channel_maker.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    channel_maker(chrom=args.chr, sampleName=args.sample, outFile=args.out)

    # inspect_windows(outFile=args.out)

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    print('Elapsed time channel_maker_real = %f' % (time() - t0))


if __name__ == '__main__':
    main()
