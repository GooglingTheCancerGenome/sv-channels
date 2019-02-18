# Imports

import argparse
import errno
import gzip
import logging
import os
import pickle
import statistics
from collections import Counter, defaultdict
from itertools import chain
from time import time

import bz2file
import numpy as np
import pyBigWig
import pysam
from functions import get_one_hot_sequence, is_outlier
from candidate_pairs import *

# import matplotlib.pyplot as plt

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = True

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


def load_clipped_read_positions(sampleName, chrName):
    channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'

    vec_type = 'clipped_read_pos'
    print('Loading CR positions for Chr %s' % chrName)
    # Load files
    if HPC_MODE:
        fn = '/'.join((sampleName, vec_type, chrName + '_' + vec_type + '.pbz2'))
    else:
        fn = '/'.join((channel_dir, sampleName, vec_type, chrName + '_' + vec_type + '.pbz2'))
    with bz2file.BZ2File(fn, 'rb') as f:
        cpos = pickle.load(f)

    cr_pos = [elem for elem, cnt in cpos.items() if cnt >= min_cr_support]

    return cr_pos


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


def get_gc_bigwig():
    bw = pyBigWig.open("/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/UCSC/hg19/hg19.gc5Base.bw")
    return bw


def get_mappability_bigwig():
    bw = pyBigWig.open("/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Mappability/GRCh37.151mer.bw")
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
    channel_names = ['candidate_pairs', 'clipped_reads_tuple', 'clipped_read_distance',
                     'coverage', 'split_read_distance_tuple']

    channel_data = defaultdict(dict)

    for chrom in chr_list:
        logging.info('Loading data for Chr%s' % chrom)
        for ch in channel_names:
            logging.info('Loading data for channel %s' % ch)
            suffix = '.npy.bz2' if ch == 'coverage' else '.pbz2'
            filename = os.path.join(prefix, sample, ch, '_'.join([chrom, ch + suffix]))
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
        channel_data[chrom]['clipped_reads_translocation'] = channel_data[chrom]['clipped_reads_tuple']
        del channel_data[chrom]['clipped_reads_tuple']

        # unpack split_reads
        channel_data[chrom]['split_read_distance'], \
        channel_data[chrom]['split_reads'] = channel_data[chrom]['split_read_distance_tuple']
        del channel_data[chrom]['split_read_distance_tuple']

    return channel_data


def channel_maker(ibam, chrom, sampleName, outFile):

    def check_progress(i, n_r, last_t):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d candidate pairs processed (%f pairs / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

    n_channels = 29
    bp_padding = 10
    channel_index = 28

    channel_data = load_channels(sampleName, [chrom])

    bw_map = get_mappability_bigwig()

    candidate_pairs_chr = [(bp1, bp2) for (bp1, bp2) in channel_data[chrom]['candidate_pairs']
                           if bp1.chr == bp2.chr and bp1.chr == chrom]

    channel_windows = np.zeros(shape=(len(candidate_pairs_chr),
                                      win_len * 2 + bp_padding, n_channels), dtype=np.uint32)

    # Consider a single sample
    sample_list = sampleName.split('_')

    for sample in sample_list:

        # Log info every n_r times
        n_r = 10 ** 3
        # print(n_r)
        last_t = time()

        for i, sv in enumerate(candidate_pairs_chr, start=0):
            check_progress(i, n_r, last_t)

            bp1, bp2 = sv

            channel_windows[i, :win_len, channel_index] = bw_map.values(bp1.chr,
                                                             bp1.pos - win_hlen, bp1.pos + win_hlen)
            channel_windows[i, win_len + bp_padding:, channel_index] = bw_map.values(bp2.chr,
                                                                          bp2.pos - win_hlen, bp2.pos + win_hlen)

    logging.info("channel_windows shape: %s" % channel_windows.shape)

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=channel_windows)
    f.close()


def channel_maker_test(ibam, chrList, sampleName, SVmode, trainingMode, outFile):
    '''
    This function loads the channels and for each clipped read position with at least min_cr_support clipped reads, it
    creates a vstack with 22 channel vectors (11 for the Tumor and 11 for the Normal sample) of width equal to
    twice the window length (win_hlen*2).
    :param ibam: BAM file used only to get the length of the chromosome from the header
    :param chrName: chromosome to consider
    :param sampleName: name of the sample pair (HMF, GiaB, or other)
    :param trainingMode: set to True only if it is used to generate the NoSV category for the Training data
    :param outFile: main output file for the list of channel vstacks
    :return: None. It saves a list of channel vstacks. If in 'NoSV' mode, it also saves a list of labels.
    '''

    # List where to store the channel vstacks
    ch_list = []

    chr_list = get_chr_len_dict(ibam).keys()

    # Get Mappability BigWig
    bw_map = get_mappability_bigwig()

    # Set the correct prefix for the path
    if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
        prefix_train = 'Training_' + SVmode + '/'
    else:
        prefix_train = ''
        # only for ovarian cancer
        # prefix_train = 'OC/'

    # Check for file existence
    if not HPC_MODE:

        # Local BAM file for testing
        sample_list = ['T0']

        clipped_read_pos_file = 'clipped_read_pos.pbz2'
        clipped_read_distance_file = 'clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads.pbz2'
        coverage_file = 'coverage.npy.gz'
        split_read_distance_file = 'split_read_distance.pbz2'

        # Check file existence
        assert os.path.isfile(clipped_read_pos_file)
        assert os.path.isfile(clipped_read_distance_file)
        assert os.path.isfile(clipped_reads_file)
        assert os.path.isfile(coverage_file)
        assert os.path.isfile(split_read_distance_file)

    else:

        # Consider a single sample
        sample_list = sampleName.split('_')

        clipped_read_pos_file = dict()
        clipped_read_distance_file = dict()
        clipped_reads_file = dict()
        coverage_file = dict()
        split_read_distance_file = dict()
        clipped_pos_cnt = dict()

        for chrName in chrList:

            # File with clipped read positions, output of the clipped_read_pos script
            # clipped_read_pos_file[chrName] = prefix_train + sample_list[0] + \
            #                                  '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
            clipped_read_pos_file[chrName] = 'clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
            # File with the clipped read distances, output of the clipped_read_distance script
            clipped_read_distance_file[chrName] = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
            # File with the clipped reads, output of the clipped_reads script
            clipped_reads_file[chrName] = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
            # File with the coverage array, output of the coverage script
            coverage_file[chrName] = 'coverage/' + chrName + '_coverage.npy.bz2'
            # File with the split reads and split read distance, output of the split_read_distance script
            split_read_distance_file[chrName] = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

            # Check file existence
            # print('Checking file: %s' % clipped_read_pos_file[chrName])
            # assert os.path.isfile(clipped_read_pos_file[chrName])

            for sample in sample_list:
                # Check file existence
                logging.info('Checking file: %s => %s' % (sample, clipped_read_pos_file[chrName]))
                assert os.path.isfile(prefix_train + sample + '/' + clipped_read_pos_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + clipped_read_distance_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + clipped_reads_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + coverage_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + split_read_distance_file[chrName])

            logging.info('Chromosome %s' % chrName)

            if len(sample_list) > 1:

                clipped_pos_cnt_per_sample = dict()
                clipped_pos = dict()

                for sample in sample_list:
                    logging.info('Reading clipped read positions for sample %s' % sample)
                    with bz2file.BZ2File(prefix_train + sample + '/' +
                                         clipped_read_pos_file[chrName], 'rb') as f:
                        clipped_pos_cnt_per_sample[sample] = pickle.load(f)
                    logging.info('End of reading')
                    logging.info('Length of clipped_pos_cnt_per_sample for sample %s: %d' % (sample,
                                                                                             len(
                                                                                                 clipped_pos_cnt_per_sample[
                                                                                                     sample])))

                    # Count the number of clipped read positions with a certain minimum number of clipped reads
                    count_clipped_read_positions(clipped_pos_cnt_per_sample[sample])

                    # if sample == sample_list[0]:
                    #     cr_support = min_cr_support
                    # else:
                    #     cr_support = 1

                    clipped_pos[sample] = [k for k, v in clipped_pos_cnt_per_sample[sample].items()
                                           if v >= min_cr_support]
                    logging.info('Length of clipped_pos_cnt_per_sample ' +
                                 ' for sample %s after min support = %d: %d' %
                                 (sample, min_cr_support, len(clipped_pos[sample])))

                clipped_pos_keep = set(clipped_pos[sample_list[0]]) - set(clipped_pos[sample_list[1]])
                logging.info('Length of cr_pos_keep: %d' % len(clipped_pos_keep))

                sample = sample_list[0]
                # clipped_pos_cnt[chrName] = {k: v for (k, v) in clipped_pos_cnt_per_sample[sample_list[0]]
                #                             if k in clipped_pos_keep}
                logging.info('Length of clipped_pos_cnt keys: %d, intersection size: %d' %
                             (len(clipped_pos_cnt_per_sample[sample].keys()),
                              len(set(clipped_pos_cnt_per_sample[sample].keys()) & clipped_pos_keep)))

                clipped_pos_cnt[chrName] = {k: clipped_pos_cnt_per_sample[sample][k]
                                            for k in clipped_pos_cnt_per_sample[sample].keys()
                                            if k in clipped_pos_keep}

                # Count the number of clipped read positions with a certain minimum number of clipped reads
                logging.info('Clipped read positions with support only in the Tumor:')
                count_clipped_read_positions(clipped_pos_cnt[chrName])

                clipped_pos[sample] = [pos for pos in clipped_pos[sample]
                                       if win_hlen <= pos <= (chrLen[chrName] - win_hlen)]
                logging.info('Length of cr_pos for sample %s after extremes removed: %d' % (sample,
                                                                                            len(clipped_pos[sample])))

            else:

                logging.info('Reading clipped read positions')
                with bz2file.BZ2File(prefix_train + sample + '/' + clipped_read_pos_file[chrName], 'rb') as f:
                    clipped_pos_cnt[chrName] = pickle.load(f)
                logging.info('End of reading')

                # Count the number of clipped read positions with a certain minimum number of clipped reads
                count_clipped_read_positions(clipped_pos_cnt[chrName])

    # Load channel data
    # Dictionaries where to load the channel data
    clipped_read_distance = dict()
    read_quality = dict()
    clipped_reads = dict()
    clipped_reads_inversion = dict()
    clipped_reads_duplication = dict()
    clipped_reads_translocation = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    outliers = dict()

    for sample in sample_list:

        prefix = prefix_train + sample + '/' if HPC_MODE else ''

        clipped_pos = dict()

        clipped_read_distance[sample] = dict()
        read_quality[sample] = dict()
        clipped_reads[sample] = dict()
        clipped_reads_inversion[sample] = dict()
        clipped_reads_duplication[sample] = dict()
        clipped_reads_translocation[sample] = dict()
        coverage[sample] = dict()
        split_reads[sample] = dict()
        split_read_distance[sample] = dict()

        outliers[sample] = dict()

        logging.info('Considering %s' % sample)

        for chrName in chrList:
            logging.info('Considering %s' % chrName)

            logging.info('Reading clipped read distances')
            with bz2file.BZ2File(prefix + clipped_read_distance_file[chrName], 'rb') as f:
                clipped_read_distance[sample][chrName] = pickle.load(f)
            logging.info('End of reading')

            logging.info('Reading clipped reads')
            with bz2file.BZ2File(prefix + clipped_reads_file[chrName], 'rb') as f:
                read_quality[sample][chrName], clipped_reads[sample][chrName], \
                clipped_reads_inversion[sample][chrName], \
                clipped_reads_duplication[sample][chrName], \
                clipped_reads_translocation[sample][chrName] = pickle.load(
                    f)
            logging.info('End of reading')

            logging.info('Reading coverage')
            with bz2file.BZ2File(prefix + coverage_file[chrName], 'rb') as f:
                coverage[sample][chrName] = np.load(file=f)
            logging.info(
                'End of reading, coverage length: %d out of %d' % (len(coverage[sample][chrName]), chrLen[chrName]))

            logging.info('Reading split read distances')
            with bz2file.BZ2File(prefix + split_read_distance_file[chrName], 'rb') as f:
                split_read_distance[sample][chrName], split_reads[sample][chrName] = pickle.load(f)
            logging.info('End of reading')

            logging.info('Finding outliers')
            outliers[sample][chrName] = dict()
            for direction in ['forward', 'reverse']:
                outliers[sample][chrName][direction] = dict()
                for clipped_arrangement in ['left', 'right', 'all']:
                    points = np.array(list(chain.from_iterable(
                        clipped_read_distance[sample][chrName][direction][clipped_arrangement].values()
                    )))
                    outlier_vec = is_outlier(points)
                    outliers[sample][chrName][direction][clipped_arrangement] = \
                        set(points[np.where(outlier_vec)].flatten())
                    # print(outliers)
            logging.info('Outliers found')

    for chrName in chrList:

        # If in 'NoSV' mode, consider all the clipped read positions (minimum clipped read support equal to 0)
        if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
            clipped_pos[chrName] = [k for k, v in clipped_pos_cnt[chrName].items()]
        else:
            clipped_pos[chrName] = [k for k, v in clipped_pos_cnt[chrName].items() if v >= min_cr_support]

    # print(clipped_pos)

    # Dictionaries where to store the channel arrays as generated from the dictionaries
    clipped_read_distance_array = dict()
    clipped_read_distance_num = dict()
    clipped_read_distance_median = dict()
    clipped_read_distance_outlier = dict()

    read_quality_array = dict()
    clipped_reads_array = dict()
    clipped_reads_inversion_array = dict()
    clipped_reads_duplication_array = dict()
    clipped_reads_translocation_array = dict()

    coverage_array = dict()

    split_read_distance_array = dict()
    split_read_distance_num = dict()
    split_read_distance_median = dict()

    split_reads_array = dict()

    gc_array = dict()
    mappability_array = dict()

    # Log info every n_r times
    n_r = 10 ** 3
    # print(n_r)
    last_t = time()

    # See how to proceed from here

    total_clipped_pos = []
    total_chr_for_clipped_pos = []

    for chrName in chrList:
        total_clipped_pos.extend(clipped_pos[chrName])
        total_chr_for_clipped_pos.extend([chrName] * len(clipped_pos[chrName]))
    print('Count of total_chr_for_clipped_pos: %s' % Counter(total_chr_for_clipped_pos))

    for i, outzipped in enumerate(zip(total_chr_for_clipped_pos, total_clipped_pos), start=1):

        chrName = outzipped[0]
        center_pos = outzipped[1]

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d clipped positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Check if center_pos is within the chromosome boundaries
        if win_hlen <= center_pos <= (chrLen[chrName] - win_hlen):

            start_win = center_pos - win_hlen
            end_win = center_pos + win_hlen

            # Build arrays for the numpy vstack
            for sample in sample_list:

                # clipped read distance
                clipped_read_distance_array[sample] = dict()
                clipped_read_distance_num[sample] = dict()
                clipped_read_distance_median[sample] = dict()
                clipped_read_distance_outlier[sample] = dict()

                for direction in ['forward', 'reverse']:
                    clipped_read_distance_array[sample][direction] = dict()
                    clipped_read_distance_num[sample][direction] = dict()
                    clipped_read_distance_median[sample][direction] = dict()
                    clipped_read_distance_outlier[sample][direction] = dict()

                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['left', 'right', 'all']:
                        clipped_read_distance_array[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                       dtype=np.uint32)
                        clipped_read_distance_num[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                     dtype=np.uint32)
                        clipped_read_distance_median[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                        dtype=np.uint32)
                        clipped_read_distance_outlier[sample][direction][clipped_arrangement] = np.zeros(
                            win_len, dtype=np.uint32)

                        for pos in range(start_win, end_win):
                            if pos in clipped_read_distance[sample][chrName][direction][clipped_arrangement].keys():
                                clipped_read_distance_array[sample][direction][clipped_arrangement][pos - start_win] = \
                                    sum(clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                clipped_read_distance_num[sample][direction][clipped_arrangement][pos - start_win] = \
                                    len(clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                clipped_read_distance_median[sample][direction][clipped_arrangement][pos - start_win] = \
                                    statistics.median(
                                        clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])

                                if direction == 'forward':

                                    clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                    (pos - start_win):] = \
                                        clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                        (pos - start_win):] + \
                                        len(set(
                                            clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                            & outliers[sample][chrName][direction][clipped_arrangement])

                                elif direction == 'reverse':

                                    clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                    :(pos - start_win)] = \
                                        clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                        :(pos - start_win)] + \
                                        len(set(
                                            clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                            & outliers[sample][chrName][direction][clipped_arrangement])

                        # print(clipped_read_distance_array[direction][clipped_arrangement])

                # read quality
                read_quality_array[sample] = read_quality[sample][chrName][start_win:end_win]

                # clipped reads
                clipped_reads_array[sample] = dict()
                for split_direction in ['left', 'right', 'D_left', 'D_right', 'I']:
                    clipped_reads_array[sample][split_direction] = np.zeros(win_len, dtype=np.uint32)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads[sample][chrName][split_direction].keys():
                            clipped_reads_array[sample][split_direction][pos - start_win] = \
                                clipped_reads[sample][chrName][split_direction][pos]

                # clipped reads inversions
                clipped_reads_inversion_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_inversion_array[sample][mate_position] = np.zeros(win_len, dtype=np.uint32)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_inversion[sample][chrName][mate_position].keys():
                            clipped_reads_inversion_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_inversion[sample][chrName][mate_position][pos]

                # clipped reads duplication
                clipped_reads_duplication_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_duplication_array[sample][mate_position] = np.zeros(win_len, dtype=np.uint32)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_duplication[sample][chrName][mate_position].keys():
                            clipped_reads_duplication_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_duplication[sample][chrName][mate_position][pos]

                # clipped reads traslocation
                clipped_reads_translocation_array[sample] = dict()
                for orientation in ['opposite', 'same']:
                    clipped_reads_translocation_array[sample][orientation] = np.zeros(win_len, dtype=np.uint32)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_translocation[sample][chrName][orientation].keys():
                            clipped_reads_translocation_array[sample][orientation][pos - start_win] = \
                                clipped_reads_translocation[sample][chrName][orientation][pos]

                # coverage
                coverage_array[sample] = coverage[sample][chrName][start_win:end_win]
                assert len(coverage_array[sample]) == win_len

                # split read distance
                split_read_distance_array[sample] = dict()
                split_read_distance_num[sample] = dict()
                split_read_distance_median[sample] = dict()

                for split_direction in ['left', 'right']:
                    split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=np.uint32)
                    split_read_distance_num[sample][split_direction] = np.zeros(win_len, dtype=np.uint32)
                    split_read_distance_median[sample][split_direction] = np.zeros(win_len, dtype=np.uint32)

                    if pos in split_read_distance[sample][chrName][split_direction].keys():
                        split_read_distance_array[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][chrName][split_direction][pos])
                        split_read_distance_num[sample][split_direction][pos - start_win] = \
                            len(split_read_distance[sample][chrName][split_direction][pos])
                        split_read_distance_median[sample][split_direction][pos - start_win] = \
                            statistics.median(split_read_distance[sample][chrName][split_direction][pos])

                # split reads
                split_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_reads_array[sample][split_direction] = np.zeros(win_len, dtype=np.uint32)
                    for pos in range(start_win, end_win):
                        if pos in split_reads[sample][chrName][split_direction].keys():
                            split_reads_array[sample][split_direction][pos - start_win] = \
                                split_reads[sample][chrName][split_direction][pos]

            # gc_array = bw_gc.values('chr' + chrName, start_win, end_win)
            # assert len(gc_array) == win_len
            mappability_array = bw_map.values(chrName, start_win, end_win)
            assert len(mappability_array) == win_len

            # Fill the numpy vstack
            vstack_list = []
            for sample in sample_list:
                # logging.info("Considering sample %s" % sample)

                vstack_list.append(coverage_array[sample])

                vstack_list.append(read_quality_array[sample])

                for clipped_arrangement in ['left', 'right', 'D_left', 'D_right', 'I']:
                    vstack_list.append(clipped_reads_array[sample][clipped_arrangement])

                for mate_position in ['before', 'after']:
                    vstack_list.append(clipped_reads_inversion_array[sample][mate_position])
                for mate_position in ['before', 'after']:
                    vstack_list.append(clipped_reads_duplication_array[sample][mate_position])
                for orientation in ['opposite', 'same']:
                    vstack_list.append(clipped_reads_translocation_array[sample][orientation])

                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['left', 'right', 'all']:
                        # vstack_list.append(
                        #     clipped_read_distance_array[sample][direction][clipped_arrangement])
                        # vstack_list.append(
                        #     clipped_read_distance_num[sample][direction][clipped_arrangement])

                        vstack_list.append(
                            clipped_read_distance_median[sample][direction][clipped_arrangement])

                        # vstack_list.append(
                        #     clipped_read_distance_outlier[sample][direction][clipped_arrangement])

                for direction in ['left', 'right']:
                    vstack_list.append(split_reads_array[sample][direction])
                for direction in ['left', 'right']:
                    # vstack_list.append(split_read_distance_array[sample][direction])
                    # vstack_list.append(split_read_distance_num[sample][direction])
                    vstack_list.append(split_read_distance_median[sample][direction])

            # vstack_list.append(gc_array)
            vstack_list.append(mappability_array)

            # append one hot encoded sequence for the genomic region
            for nuc in ['A', 'T', 'C', 'G', 'N']:
                one_hot_n = get_one_hot_sequence(chrName, start_win, end_win, nuc, HPC_MODE)
                assert len(one_hot_n) == win_len
                vstack_list.append(one_hot_n)

            # logging.info("Shape of channel matrix: %s" % str(ch_vstack.shape))
            ch_vstack = np.vstack(vstack_list)
            ch_list.append(ch_vstack)

    outDir = os.path.dirname(outFile)
    labelDir = outDir + '/label/'
    create_dir(labelDir)

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    # Write labels for noSV category
    # if trainingMode and sampleName == 'noSV':
    if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
        label = ['noSV'] * len(ch_list)
        with gzip.GzipFile(labelDir + sampleName + '_' + chrName + '_label.npy.gz', "w") as f:
            np.save(file=f, arr=label)
        f.close()

    logging.info('Number of windows: %d' % len(ch_list))


def main():
    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    # Path on the HPC for the test BAM file
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + 'T0_dedup.bam'

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='2',
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
        level=logging.INFO)

    t0 = time()

    channel_maker(ibam=args.bam, chrom=args.chr, sampleName=args.sample, outFile=args.out)

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    print('Elapsed time channel_maker_real = %f' % (time() - t0))


if __name__ == '__main__':
    main()
