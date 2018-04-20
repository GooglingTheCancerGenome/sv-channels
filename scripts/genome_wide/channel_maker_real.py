import pysam
import numpy as np
import argparse
import bz2file
import gzip
import os
import pickle
from time import time
import logging
import csv
from collections import Counter


HPC_MODE = False


def count_clipped_read_positions(cpos_cnt):
    for i in range(0, 5):
        logging.info('Number of positions with at least %d clipped reads: %d' %
              (i + 1, len([k for k, v in cpos_cnt.items() if v > i])))


def read_gold_positions():

    gold_del_file = '/Users/lsantuari/Documents/Data/GiaB/GoldStandard/DEL_Mills2011.csv'
    gold_insdup_file = '/Users/lsantuari/Documents/Data/GiaB/GoldStandard/INS_DUP_Mills2011.csv'

    del_start = dict()
    del_end = dict()
    ins_start = dict()

    with open(gold_del_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quotechar='#')
        for row in reader:
            #print row
            if row[0]=='NA12878':
                chr = row[1][3:]
                # chr23 is chrX
                if chr == '23':
                    chr = 'X'
                if chr in del_start.keys():
                    del_start[chr].append(int(row[2]))
                else:
                    del_start[chr] = [int(row[2])]
                if chr in del_end.keys():
                     del_end[chr].append(int(row[2]))
                else:
                    del_end[chr] = [int(row[2])]
    with open(gold_insdup_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quotechar='#')
        for row in reader:
            #print(row[4])
            if row[0]=='NA12878' and row[4] == 'INS':
                chr = row[1]
                if chr in ins_start.keys():
                    ins_start[chr].append(int(row[2]))
                else:
                    ins_start[chr] = [int(row[2])]
    # print(del_start.keys())
    return del_start, del_end, ins_start


def read_bpi_positions():

    gold_del_file = '/Users/lsantuari/Documents/Data/HMF/HMF_COLO829_VCF/bpi/COLO829R_COLO829T_bpi_stats.tsv'

    del_start = dict()
    del_end = dict()

    with open(gold_del_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            #print(row)
            if row[1] == 'DEL':
                chrA, start = row[3].split(':')
                chrB, end = row[4].split(':')
                #print("%s %s %s %s" % (chrA, start, chrB, end))
                assert chrA == chrB
                if chrA in del_start.keys():
                    del_start[chrA].append(int(start))
                else:
                    del_start[chrA] = [int(start)]
                if chrB in del_end.keys():
                    del_end[chrB].append(int(end))
                else:
                    del_end[chrB] = [int(end)]

    # print(del_start.keys())
    return del_start, del_end


def read_BED():

    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    truth_file = wd + "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed"

    #wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_indel/"
    #truth_file = wd + "SV/chr17_INDEL.bed"

    chr_list_DEL = []
    start_SV_DEL = []
    end_SV_DEL = []

    chr_list_INS = []
    start_SV_INS = []

    with(open(truth_file, 'r')) as i:
        for line in i:
            line = line.rstrip()
            columns = line.split("\t")
            if columns[4] == "DEL":
                start_SV_DEL.append(int(columns[1]))
                end_SV_DEL.append(int(columns[3]))
                chr_list_DEL.append(int(columns[0]))
            elif str(columns[4]) == "INS":
                start_SV_INS.append(int(columns[1]))
                chr_list_INS.append(int(columns[0]))

    return start_SV_DEL, end_SV_DEL, start_SV_INS


def compare_positions():

    min_cr_support = 3
    confint = 100

    no_clipped_read_pos = dict()

    # del_start, del_end, ins_start = read_gold_positions()
    del_start, del_end = read_bpi_positions()
    #print(sum(map(len, list(del_start.values()))))
    #print(del_start)
    #print(del_end)

    sum_cr = 0

    for chr in del_start.keys():
        print('Considering chromosome %s' % chr)
        # clipped_read_pos_file = '/Users/lsantuari/Documents/Data/GiaB/ChannelMaker_results/T/clipped_read_pos/' \
        #                        + chr + '_clipped_read_pos.pbz2'

        clipped_read_pos_file = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/clipped_read_pos/' \
                                + chr + '_clipped_read_pos.pbz2'

        logging.info('Reading clipped read positions')
        with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
            clipped_pos_cnt = pickle.load(f)
        logging.info('End of reading')

        # Select clipped read positions with minimum support of min_cr_support clipped reads
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

        sum_cr += len(clipped_pos)
        print(len(clipped_pos))

        if chr in del_start.keys():
            start_bpj_lst = []
            for pos in set(del_start[chr]):
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                if len(int_set) > 0:
                    start_bpj_lst.extend(int_set)
                else:
                    #print(pos)
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        if chr in del_end.keys():
            end_bpj_lst = []
            for pos in set(del_end[chr]):
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                if len(int_set) > 0:
                    end_bpj_lst.extend(int_set)
                else:
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        # print(set(start_bpj_lst) & set(end_bpj_lst))
        # assert len(set(start_bpj_lst) & set(end_bpj_lst)) == 0
        labels = []
        for pos in clipped_pos:
            if pos in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_start')
            elif pos in start_bpj_lst and pos not in end_bpj_lst:
                labels.append('del_start')
            elif pos not in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_end')
            else:
                labels.append('no_del')

        assert len(clipped_pos) == len(labels)

        outFile = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                                    + chr + '_labels.npy.gz'
        with gzip.GzipFile(outFile, "w") as f:
            np.save(file=f, arr=labels)
        f.close()


    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                                + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "w") as f:
        pickle.dump(no_clipped_read_pos, f)
    f.close()

    print("Number of clipped positions: %d" % sum_cr)
    print("No clipped positions: %d" % len(list(no_clipped_read_pos.values())))

            # print(set(clipped_pos) & set(map(lambda x:x, del_end[chr])))
        # print(ins_start.keys())

        # if chr in ins_start.keys():
        #    print('ins_start:')
        #    print(set(clipped_pos) & set(map(lambda x:x-1, ins_start[chr])))


def load_NoCR_positions():

    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "r") as f:
        no_clipped_read_pos = pickle.load(f)
    f.close()
    print(list(no_clipped_read_pos))
    #print("No clipped positions: %d" % len(list(no_clipped_read_pos.values())))


def channel_maker(ibam, chrName, sampleName, trainingMode, outFile):

    ch_list = []

    assert os.path.isfile(ibam)
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromossome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    win_hlen = 25
    win_len = win_hlen * 2
    min_cr_support = 3

    # Check for file existence
    if not HPC_MODE:

        sample_list = ['GS']

        clipped_read_pos_file = 'clipped_read_pos.pbz2'
        clipped_read_distance_file = 'clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads.pbz2'
        coverage_file = 'coverage.npy.gz'
        split_read_distance_file = 'split_read_distance.pbz2'

        assert os.path.isfile(clipped_read_pos_file)
        assert os.path.isfile(clipped_read_distance_file)
        assert os.path.isfile(clipped_reads_file)
        assert os.path.isfile(coverage_file)
        assert os.path.isfile(split_read_distance_file)

    else:

        # sample_list = ['T', 'N']

        # For HMF sample
        sample_list = ['Tumor', 'Normal']
        # For NoSV sample
        # sample_list = ['NoSV/Tumor', 'NoSV/Normal']

        clipped_read_pos_file = sampleName + '/' + sample_list[0] + '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'

        clipped_read_distance_file = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
        coverage_file = 'coverage/' + chrName + '_coverage.npy.gz'
        split_read_distance_file = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

        assert os.path.isfile(clipped_read_pos_file)

        for sample in sample_list:
            assert os.path.isfile(sampleName + '/' + sample + '/' + clipped_read_distance_file)
            assert os.path.isfile(sampleName + '/' + sample + '/' + clipped_reads_file)
            assert os.path.isfile(sampleName + '/' + sample + '/' + coverage_file)
            assert os.path.isfile(sampleName + '/' + sample + '/' + split_read_distance_file)


    logging.info('Chromosome %s' % chrName)

    logging.info('Reading clipped read positions')
    with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
        #clipped_pos_cnt, clipped_read_1, clipped_read_2 = pickle.load(f)
        clipped_pos_cnt = pickle.load(f)
    logging.info('End of reading')

    count_clipped_read_positions(clipped_pos_cnt)

    # Load channel data
    clipped_read_distance = dict()
    clipped_reads = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    for sample in sample_list:

        prefix = sampleName + '/' + sample + '/' if HPC_MODE else ''

        logging.info('Considering %s' % sample)
        logging.info('Reading clipped read distances')
        with bz2file.BZ2File(prefix + clipped_read_distance_file, 'rb') as f:
            clipped_read_distance[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading clipped reads')
        with bz2file.BZ2File(prefix + clipped_reads_file, 'rb') as f:
            clipped_reads[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading coverage')
        with bz2file.BZ2File(prefix + coverage_file, 'rb') as f:
            coverage[sample] = np.load(file=f)
        logging.info('End of reading, coverage length: %d out of %d' % (len(coverage[sample]), chrLen))

        logging.info('Reading split read distances')
        with bz2file.BZ2File(prefix + split_read_distance_file, 'rb') as f:
            split_read_distance[sample], split_reads[sample] = pickle.load(f)
        logging.info('End of reading')

    # cpos_cnt = Counter(clipped_pos_cnt.values())
    # print(clipped_pos_cnt)

    if trainingMode and sampleName == 'noSV':
        clipped_pos = [k for k, v in clipped_pos_cnt.items()]
    else:
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

    # print(clipped_pos)

    clipped_read_distance_array = dict()
    clipped_reads_array = dict()
    coverage_array = dict()
    split_read_distance_array = dict()
    split_reads_array = dict()

    n_r = 10 ** 3
    # print(n_r)
    last_t = time()

    for i, center_pos in enumerate(clipped_pos, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d clipped positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Check if center_pos is within the chromosome boundaries
        if win_hlen <= center_pos <= (chrLen - win_hlen):

            start_win = center_pos - win_hlen
            end_win = center_pos + win_hlen

            # Build arrays for the numpy vstack
            for sample in sample_list:

                # clipped read distance
                clipped_read_distance_array[sample] = dict()
                for direction in ['forward', 'reverse']:
                    clipped_read_distance_array[sample][direction] = dict()
                for direction in ['forward', 'reverse']:
                    # for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                    for clipped_arrangement in ['left', 'right']:
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
                # coverage_array[sample] = np.zeros(win_len, dtype=int)
                # for pos in range(start_win, end_win):
                #    if pos in coverage[sample].keys():
                #        coverage_array[sample][pos - start_win] = coverage[sample][pos]
                coverage_array[sample] = coverage[sample][start_win:end_win]
                assert len(coverage_array[sample]) == win_len

                # split read distance
                split_read_distance_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    if pos in split_read_distance[sample][split_direction].keys():
                        split_read_distance_array[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][split_direction][pos])

                # split reads
                split_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in split_reads[sample][split_direction].keys():
                            split_reads_array[sample][split_direction][pos - start_win] = \
                                split_reads[sample][split_direction][pos]

            # Fill the numpy vstack
            for sample in sample_list:
                # logging.info("Considering sample %s" % sample)

                if sample == sample_list[0]:
                    ch_vstack = np.vstack((
                        coverage_array[sample],
                        clipped_reads_array[sample]['left'],
                        clipped_reads_array[sample]['right']))
                else:
                    ch_vstack = np.vstack((ch_vstack,
                        coverage_array[sample],
                        clipped_reads_array[sample]['left'],
                        clipped_reads_array[sample]['right']))

                for direction in ['forward', 'reverse']:
                    # for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                    for clipped_arrangement in ['left', 'right']:
                        ch_vstack = np.vstack((ch_vstack,
                                               clipped_read_distance_array[sample][direction][
                                                   clipped_arrangement]))
                for direction in ['left', 'right']:
                    ch_vstack = np.vstack((ch_vstack,
                                           split_reads_array[sample][direction]))
                for direction in ['left', 'right']:
                    ch_vstack = np.vstack((ch_vstack,
                                           split_read_distance_array[sample][direction]))

            # logging.info("Shape of channel matrix: %s" % str(ch_vstack.shape))
            ch_list.append(ch_vstack)

    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    logging.info('Number of windows: %d' % len(ch_list))


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    # wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_indel/"
    # inputBAM = wd + "BAM/S1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output")
    parser.add_argument('-s', '--sample', type=str, default='somatic',
                        help="Specify sample")
    parser.add_argument('-t', '--train', type=bool, default=True,
                        help="Specify if training mode is active")
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
    channel_maker(ibam=args.bam, chrName=args.chr, sampleName=args.sample,
                  trainingMode=args.train, outFile=args.out)

    # Compare GiaB positions with the clipped read positions
    # compare_positions()
    # load_NoCR_positions()
    # start_SV_DEL, end_SV_DEL, start_SV_INS = read_BED()

    # print(start_SV_DEL)
    # print(end_SV_DEL)
    # print(start_SV_INS)

    print('Elapsed time = %f' % (time() - t0))


if __name__ == '__main__':
    main()
