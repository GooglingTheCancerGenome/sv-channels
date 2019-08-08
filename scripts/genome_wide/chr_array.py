# Imports

import argparse
import errno
import gzip
import logging
import os
import json
import statistics
from collections import defaultdict
from time import time
import numpy as np
import pyBigWig
import pysam
import dask.array as da
from functions import *

config = get_config_file()
HPC_MODE = config["DEFAULT"]["HPC_MODE"]
REF_GENOME = config["DEFAULT"]["REF_GENOME"]


def get_chr_len(ibam, chrName):
    # check if the BAM file exists
    assert os.path.isfile(ibam), ibam + " file not existent!"
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


def load_candidate_positions(datadir, sampleName, chrom1, chrom2):
    vec_type = 'split_read_pos'

    print('Loading candidate positions for Chr%s to Chr%s' % (chrom1, chrom2))

    # Load files
    fn = os.path.join(datadir, sampleName, vec_type, chrom1 + '_' + vec_type + '.json.gz')

    with gzip.GzipFile(fn, 'r') as fin:
        positions, locations = json.loads(fin.read().decode('utf-8'))

    locations = [(c1, p1, c2, p2) for (c1, p1, c2, p2) in locations if c1 == chrom1 and c2 == chrom2]

    print('%d candidate positions for Chr%s to Chr%s' % (len(locations), chrom1, chrom2))

    return locations


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
    mappability_file = os.path.join("/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Mappability",
                                    REF_GENOME, REF_GENOME+".151mer.bw") if HPC_MODE \
        else os.path.join("/Users/lsantuari/Documents/Data/GEM", REF_GENOME, REF_GENOME+".151mer.bw")
    bw = pyBigWig.open(mappability_file)

    return bw


def load_bam(ibam):
    # check if the BAM file exists
    assert os.path.isfile(ibam), ibam + " file not existent!"
    # open the BAM file
    return pysam.AlignmentFile(ibam, "rb")


def get_chr_len_dict(ibam):
    bamfile = load_bam(ibam)
    # Extract chromosome length from the BAM header
    header_dict = bamfile.header

    chrLen = {i['SN']: i['LN'] for i in header_dict['SQ']}
    return chrLen


def load_channels(sample, chr_list, outDir):
    channel_names = ['coverage', 'clipped_reads', 'clipped_read_distance',
                     'snv', 'split_read_distance']

    channel_data = defaultdict(dict)

    for chrom in chr_list:
        logging.info('Loading data for Chr%s' % chrom)

        for ch in channel_names:

            logging.info('Loading data for channel %s' % ch)
            suffix = '.npy.gz' if ch in ['snv', 'coverage'] else '.json.gz'
            filename = os.path.join(outDir, sample, ch, '_'.join([chrom, ch + suffix]))
            assert os.path.isfile(filename), filename + " does not exists!"

            logging.info('Reading %s for Chr%s' % (ch, chrom))

            if suffix == '.npy.gz':
                with gzip.GzipFile(filename, 'r') as fin:
                    channel_data[chrom][ch] = np.load(fin)
                channel_data[chrom][ch] = np.swapaxes(channel_data[chrom][ch], 0, 1)
            else:
                with gzip.GzipFile(filename, 'r') as fin:
                    channel_data[chrom][ch] = json.loads(fin.read().decode('utf-8'))
            logging.info('End of reading')

        # unpack clipped_reads
        channel_data[chrom]['clipped_reads'], \
        channel_data[chrom]['clipped_reads_inversion'], channel_data[chrom]['clipped_reads_duplication'], \
        channel_data[chrom]['clipped_reads_translocation'] = channel_data[chrom]['clipped_reads']

        # unpack split_reads
        channel_data[chrom]['split_read_distance'], \
        channel_data[chrom]['split_reads'] = channel_data[chrom]['split_read_distance']

    return channel_data


def create_hdf5(sampleName, ibam, chrom, outDir, cmd_name):
    chrlen = get_chr_len(ibam, chrom)
    n_channels = 33

    channel_data = load_channels(sampleName, [chrom], outDir)
    chr_array = np.zeros(shape=(chrlen, n_channels), dtype=np.float32)

    bw_map = get_mappability_bigwig()

    # dictionary of key choices
    direction_list = {'clipped_reads': ['left', 'right', 'D_left', 'D_right', 'I'],
                      'split_reads': ['left', 'right'],
                      'split_read_distance': ['left', 'right'],
                      'clipped_reads_inversion': ['before', 'after'],
                      'clipped_reads_duplication': ['before', 'after'],
                      'clipped_reads_translocation': ['opposite', 'same'],
                      'clipped_read_distance': ['forward', 'reverse']
                      }

    channel_index = 0
    for current_channel in ['coverage', 'snv',
                            'clipped_reads', 'split_reads',
                            'clipped_reads_inversion', 'clipped_reads_duplication',
                            'clipped_reads_translocation',
                            'clipped_read_distance', 'split_read_distance']:

        logging.info("Adding channel %s at index %d" % (current_channel, channel_index))

        if current_channel == 'coverage' or current_channel == 'snv':

            ch_num = channel_data[chrom][current_channel].shape[1]
            chr_array[:, channel_index:channel_index + ch_num] = channel_data[chrom][current_channel][:chrlen, :]
            channel_index += ch_num
            del channel_data[chrom][current_channel]

        elif current_channel in ['clipped_reads',
                                 'split_reads',
                                 'clipped_reads_inversion',
                                 'clipped_reads_duplication',
                                 'clipped_reads_translocation']:

            for split_direction in direction_list[current_channel]:

                if len(channel_data[chrom][current_channel][split_direction]) > 0:

                    # print(split_direction)
                    idx = np.array(
                        list(
                            map(
                                int,
                                channel_data[chrom][current_channel][split_direction].keys()
                                )
                        )
                    )

                    vals = np.array(
                        list(
                            channel_data[chrom][current_channel][split_direction].values()
                        )
                    )

                    chr_array[idx, channel_index] = vals

                channel_index += 1
                del channel_data[chrom][current_channel][split_direction]

        elif current_channel == 'clipped_read_distance':
            for split_direction in direction_list[current_channel]:
                for clipped_arrangement in ['left', 'right', 'all']:

                    idx = np.array(
                        list(
                            map(
                                int,
                                channel_data[chrom][current_channel][split_direction][
                                    clipped_arrangement].keys()
                            )
                        )
                    )
                    vals = np.array(
                        list(
                            map(
                                statistics.median,
                                             channel_data[chrom][current_channel][split_direction][
                                                 clipped_arrangement].values()
                            )
                        )
                    )

                    chr_array[idx, channel_index] = vals

                    channel_index += 1
                    del channel_data[chrom][current_channel][split_direction][clipped_arrangement]

        elif current_channel == 'split_read_distance':
            for split_direction in direction_list[current_channel]:

                idx = np.array(
                    list(
                        map(
                            int,
                            channel_data[chrom][
                                          current_channel][split_direction].keys()
                        )
                    )
                )
                vals = np.array(
                    list(
                        map(
                            statistics.median,
                                         channel_data[chrom][
                                             current_channel][split_direction].values()
                        )
                    )
                )

                chr_array[idx, channel_index] = vals

                channel_index += 1
                del channel_data[chrom][current_channel][split_direction]

    current_channel = 'mappability'
    logging.info("Adding channel %s at index %d" % (current_channel, channel_index))

    chr_array[:, channel_index] = np.array(bw_map.values(chrom, 0, chrlen), dtype=np.float32)
    channel_index += 1

    current_channel = 'one_hot_encoding'
    logging.info("Adding channel %s at index %d" % (current_channel, channel_index))

    nuc_list = ['A', 'T', 'C', 'G', 'N']

    chr_array[:, channel_index:channel_index + len(nuc_list)] = get_one_hot_sequence_by_list(
        chrom, list(np.arange(chrlen)), HPC_MODE)
    channel_index += len(nuc_list)

    logging.info("chr_array shape: %s" % str(chr_array.shape))

    dask_array = da.from_array(chr_array, chunks=("auto", -1))

    outfile = os.path.join(outDir, sampleName, cmd_name, sampleName + '_' + chrom + '.hdf5')

    logging.info("Writing HDF5...")

    da.to_hdf5(outfile, '/' + 'chr' + chrom, dask_array)  # , compression='lzf', shuffle=False)


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
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-s', '--sample', type=str, default='T1',
                        help="Specify sample")
    parser.add_argument('-l', '--logfile', default='window_maker.log',
                        help='File in which to write logs.')
    parser.add_argument('-w', '--window', type=str, default=200,
                        help="Specify window size")

    args = parser.parse_args()

    cmd_name = 'chr_array'
    output_dir = os.path.join(args.outputpath, args.sample, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, args.logfile)
    # output_file = os.path.join(output_dir, args.out)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    create_hdf5(sampleName=args.sample,
                ibam=args.bam,
                chrom=args.chr,
                outDir=args.outputpath,
                cmd_name=cmd_name
                )

    # inspect_windows(outFile=args.out)
    # chr_list = [args.chr1]
    # load_channels(sample=args.sample, chr_list=chr_list, outDir=args.outputpath)

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    logging.info('Elapsed time channel_maker_real = %f mins' % (time() - t0))


if __name__ == '__main__':
    main()
