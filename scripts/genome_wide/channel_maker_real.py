import pysam
import numpy as np
import argparse
import bz2file
import gzip
import os
import pickle
from time import time
import logging

HPC_MODE = False


def count_clipped_read_positions(cpos_cnt):
    for i in range(2, 5):
        print('Number of positions with at least %d clipped reads: %d' %
              (i + 1, sum([v for k, v in cpos_cnt.items() if k > i])))


def channel_maker(ibam, chrName, outFile):

    ch_list = []

    assert os.path.isfile(ibam)
    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    win_hlen = 100
    win_len = win_hlen * 2

    if not HPC_MODE:

        sample_list = ['GS']

        clipped_read_pos_file = 'clipped_read_pos.pbz2'
        clipped_read_distance_file = 'clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads.pbz2'
        coverage_file = 'coverage.pbz2'
        split_read_distance_file = 'split_read_distance.pbz2'

        assert os.path.isfile(clipped_read_pos_file)
        assert os.path.isfile(clipped_read_distance_file)
        assert os.path.isfile(clipped_reads_file)
        assert os.path.isfile(coverage_file)
        assert os.path.isfile(split_read_distance_file)

    else:

        sample_list = ['T', 'N']
        clipped_read_pos_file = 'T/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'

        clipped_read_distance_file = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
        coverage_file = 'coverage/' + chrName + '_coverage.pbz2'
        split_read_distance_file = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

        assert os.path.isfile(clipped_read_pos_file)

        for sample in sample_list:
            assert os.path.isfile(sample + '/' + clipped_read_distance_file)
            assert os.path.isfile(sample + '/' + clipped_reads_file)
            assert os.path.isfile(sample + '/' + coverage_file)
            assert os.path.isfile(sample + '/' + split_read_distance_file)

    logging.info('Chromosome %s' % chrName)

    logging.info('Reading clipped read positions')
    with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
        clipped_pos_cnt, clipped_read_1, clipped_read_2 = pickle.load(f)
    logging.info('End of reading')

    clipped_read_distance = dict()
    clipped_reads = dict()
    coverage = dict()
    split_read_distance = dict()

    for sample in sample_list:
        logging.info('Considering %s' % sample)
        logging.info('Reading clipped read distances')
        with bz2file.BZ2File(clipped_read_distance_file, 'rb') as f:
            clipped_read_distance[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading clipped reads')
        with bz2file.BZ2File(clipped_reads_file, 'rb') as f:
            clipped_reads[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading coverage')
        with bz2file.BZ2File(coverage_file, 'rb') as f:
            coverage[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading split read distances')
        with bz2file.BZ2File(split_read_distance_file, 'rb') as f:
            split_read_distance[sample] = pickle.load(f)
        logging.info('End of reading')

    # cpos_cnt = Counter(clipped_pos_cnt.values())
    # print(clipped_pos_cnt)
    clipped_pos = [k for k, v in clipped_pos_cnt.items() if v > 2]
    # print(clipped_pos)

    clipped_read_distance_array = dict()
    clipped_reads_array = dict()
    coverage_array = dict()
    split_read_distance_array = dict()

    for center_pos in clipped_pos:
        if win_hlen <= center_pos <= (chrLen - win_hlen):

            start_win = center_pos - win_hlen
            end_win = center_pos + win_hlen

            for sample in sample_list:

                # clipped read distance
                clipped_read_distance_array[sample] = dict()
                for direction in ['forward', 'reverse']:
                    clipped_read_distance_array[sample][direction] = dict()
                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                        clipped_read_distance_array[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                       dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in clipped_read_distance[sample][direction][clipped_arrangement].keys():
                                clipped_read_distance_array[sample][direction][clipped_arrangement][pos - start_win] = \
                                    sum(clipped_read_distance[sample][direction][clipped_arrangement][pos])
                        # print(clipped_read_distance_array[direction][clipped_arrangement])

                # clipped reads
                clipped_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    clipped_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads[sample][split_direction].keys():
                            clipped_reads_array[sample][split_direction][pos - start_win] = \
                                clipped_reads[sample][split_direction][pos]

                # coverage
                coverage_array[sample] = np.array(coverage[sample][start_win:end_win])

                # split read distance
                split_read_distance_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    if pos in split_read_distance[sample][split_direction].keys():
                        split_read_distance_array[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][split_direction][pos])

            for sample in sample_list:

                ch_vstack = np.vstack((
                    clipped_reads_array[sample]['left'],
                    clipped_reads_array[sample]['right']))

                ch_vstack = np.vstack((ch_vstack,
                                       coverage_array[sample]))
                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                        ch_vstack = np.vstack((ch_vstack,
                                               clipped_read_distance_array[sample][direction][
                                                   clipped_arrangement]))
                for direction in ['left', 'right']:
                    ch_vstack = np.vstack((ch_vstack,
                                           split_read_distance_array[sample][direction]))

            ch_list.append(ch_vstack)

    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    logging.info('Number of windows: %d' % len(ch_list))


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output")
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
    channel_maker(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
