# Imports

import pysam
from pysam import VariantFile
import numpy as np
import argparse
import bz2file
import gzip
import os
import pickle
from time import time
import logging
import csv

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = True


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


def read_gold_positions():
    '''
    Deprecated. Returns the positions of the deletions and insertions of the GiaB dataset NA12878 from the paper
    Mills et al., 2011
    :return: del_start: dictionary of start positions of deletions per chromosome
             del_end:  dictionary of end positions of deletions per chromosome. Corresponds to the
                       start positions in del_start
             ins_start: dictionary of start positions of insertions per chromosome
    '''
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
                # chr23 is chrX, rename it properly
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
    '''
    Reads Manta/BPI TSV file for the COLO829 HMF dataset
    :return: Two dictionaries: del_start and del_end, with start and end positions of deletions per chromosome
    '''
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


def generate_labels():

    '''
    For each chromosome, load the clipped read positions file (clipped_read_pos) and label each position according
    to the
    :return:
    '''

    # Only clipped read positions supported by at least min_cr_support clipped reads are considered
    min_cr_support = 3
    # Half of the window length centered on SV positions, where to consider clipped read positions
    confint = 25

    # Keep track of the windows where no clipped reads can be found
    # How many SV positions reported by the SV caller do not contain any clipped reads?
    no_clipped_read_pos = dict()

    # Deletions and insertions of the GiaB dataset
    # del_start, del_end, ins_start = read_gold_positions()

    # consider HMF COLO829 Manta/BPI positions
    # del_start, del_end = read_bpi_positions()

    # read Lumpy positions
    del_start, del_end = get_Lumpy_positions()

    sum_cr = 0

    # Iterate over chromosomes
    for chr in del_start.keys():
        print('Considering chromosome %s' % chr)

        # Local path to clipped read positions
        clipped_read_pos_file = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/clipped_read_pos/' \
                                + chr + '_clipped_read_pos.pbz2'

        logging.info('Reading clipped read positions')
        with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
            clipped_pos_cnt = pickle.load(f)
        logging.info('End of reading')

        # Select clipped read positions with minimum support of min_cr_support clipped reads
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

        # Consider positions referring to start of deletions
        if chr in del_start.keys():
            start_bpj_lst = []
            for pos in set(del_start[chr]):
                # Take the intersection between the position in the [pos-confint,pos+confint] window positions and the
                # clipped read positions
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                # If some clipped read positions are included in the window, add them to the set of candidate breakpoint
                # junctions
                if len(int_set) > 0:
                    start_bpj_lst.extend(int_set)
                else:
                    # If there are no clipped read positions, record the SV position in the no_clipped_read_pos
                    # dictionary
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        # Consider positions referring to end of deletions
        if chr in del_end.keys():
            end_bpj_lst = []
            for pos in set(del_end[chr]):
                # Take the i§ntersection between the position in the [pos-confint,pos+confint] window positions and the
                # clipped read positions
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                # If some clipped read positions are included in the window, add them to the set of candidate breakpoint
                # junctions
                if len(int_set) > 0:
                    end_bpj_lst.extend(int_set)
                else:
                    # If there are no clipped read positions, record the SV position in the no_clipped_read_pos
                    # dictionary
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        # print(set(start_bpj_lst) & set(end_bpj_lst))
        # assert len(set(start_bpj_lst) & set(end_bpj_lst)) == 0
        labels = []
        for pos in clipped_pos:
            # Is a clipped read position both in the list of start and end candidate positions?
            if pos in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_start_and_end')
            # The clipped read position belongs only to the set of start breakpoint candidates
            elif pos in start_bpj_lst and pos not in end_bpj_lst:
                labels.append('del_start')
            # The clipped read position belongs only to the set of end breakpoint candidates
            elif pos not in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_end')
            # The clipped read position does not belong to any deletion endpoint positions
            else:
                labels.append('no_del')

        # A label should be assigned to each clipped read position
        assert len(clipped_pos) == len(labels)

        # (Locally) save the labels per chromosome
        outFile = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/labels/' \
                                    + chr + '_labels.npy.gz'
        with gzip.GzipFile(outFile, "w") as f:
            np.save(file=f, arr=labels)
        f.close()

    # File with SV positions without clipped reads
    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/labels/' \
                                + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "w") as f:
        pickle.dump(no_clipped_read_pos, f)
    f.close()

    print("Number of clipped positions: %d" % sum_cr)
    print("No clipped positions: %d" % len(list(no_clipped_read_pos.values())))


def get_Lumpy_positions():
    '''
    This function returns two dictionaries: del_start and del_end. They contain respectively start and end positions
    of the deletions (SVTYPE:DEL) present in the Lumpy output file.
    :return: del_start and del_end dictionaries
    '''

    del_start = dict()
    del_end = dict()

    # Local Lumpy VCF output file for the Z424 dataset
    vcf_file = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/VCF/Z424/Lumpy.vcf'
    vcf_in = VariantFile(vcf_file)

    for rec in vcf_in.fetch():
        if rec.info['SVTYPE']=='DEL' and rec.chrom not in ['Y', 'MT']:
            #print (str(rec.chrom)+' '+str(rec.pos)+' '+str(rec.stop))

            if rec.chrom in del_start.keys():
                del_start[rec.chrom].append(rec.pos)
            else:
                del_start[rec.chrom] = [rec.pos]

            if rec.chrom in del_end.keys():
                del_end[rec.chrom].append(rec.stop)
            else:
                del_end[rec.chrom] = [rec.stop]
    #print(del_start.keys())
    return  del_start, del_end


def load_NoCR_positions():

    '''
    This function provides an overview of SV positions without clipped read support that are stored in the
    no_clipped_read_positions file.
    :return: None
    '''

    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "r") as f:
        no_clipped_read_pos = pickle.load(f)
    f.close()
    print(list(no_clipped_read_pos))


def channel_maker(ibam, chrName, sampleName, trainingMode, outFile):
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
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Use a 50 bp window length
    # TODO: the window length should be an input parameter
    win_hlen = 25
    win_len = win_hlen * 2
    # Minimum clipped read support to consider
    min_cr_support = 3

    # Set the correct prefix for the path
    if trainingMode and sampleName == 'noSV':
        prefix_train = 'Training/'
    else:
        prefix_train = ''

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

        # Consider a Tumor/Normal pair
        sample_list = ['Tumor', 'Normal']

        # File with clipped read positions, output of the clipped_read_pos script
        clipped_read_pos_file = prefix_train + sampleName + '/' + sample_list[0] + '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
        # File with the clipped read distances, output of the clipped_read_distance script
        clipped_read_distance_file = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
        # File with the clipped reads, output of the clipped_reads script
        clipped_reads_file = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
        # File with the coverage array, output of the coverage script
        coverage_file = 'coverage/' + chrName + '_coverage.npy.bz2'
        # File with the split reads and split read distance, output of the split_read_distance script
        split_read_distance_file = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

        # Check file existence
        print('Checking file: %s' % clipped_read_pos_file)
        assert os.path.isfile(clipped_read_pos_file)

        for sample in sample_list:
            assert os.path.isfile(prefix_train + sampleName + '/' + sample + '/' + clipped_read_distance_file)
            assert os.path.isfile(prefix_train + sampleName + '/' + sample + '/' + clipped_reads_file)
            assert os.path.isfile(prefix_train + sampleName + '/' + sample + '/' + coverage_file)
            assert os.path.isfile(prefix_train + sampleName + '/' + sample + '/' + split_read_distance_file)

    logging.info('Chromosome %s' % chrName)

    logging.info('Reading clipped read positions')
    with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
        clipped_pos_cnt = pickle.load(f)
    logging.info('End of reading')

    # Count the number of clipped read positions with a certain minimum number of clipped reads
    count_clipped_read_positions(clipped_pos_cnt)

    # Load channel data
    # Dictionaries where to load the channel data
    clipped_read_distance = dict()
    clipped_reads = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    for sample in sample_list:

        prefix = prefix_train + sampleName + '/' + sample + '/' if HPC_MODE else ''

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

    # If in 'NoSV' mode, consider all the clipped read positions (minimum clipped read support equal to 0)
    if trainingMode and sampleName == 'noSV':
        clipped_pos = [k for k, v in clipped_pos_cnt.items()]
    else:
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

    # print(clipped_pos)

    # Dictionaries where to store the channel arrays as generated from the dictionaries
    clipped_read_distance_array = dict()
    clipped_reads_array = dict()
    coverage_array = dict()
    split_read_distance_array = dict()
    split_reads_array = dict()

    # Log info every n_r times
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
            # TODO: avoid for loop for filling the vstack. Use a list instead.
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

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    # Write labels for noSV category
    if trainingMode and sampleName == 'noSV':
        label = ['noSV']*len(ch_list)
        with gzip.GzipFile(sampleName + '_label.npy.gz', "w") as f:
            np.save(file=f, arr=label)
        f.close()

    logging.info('Number of windows: %d' % len(ch_list))


def main():

    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    # Local path
    #wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_indel/"
    #inputBAM = wd + "BAM/S1_dedup.bam"

    # Path on the HPC for the test BAM file
    wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    inputBAM = wd + 'T0_dedup.bam'

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

    # Generate labels for the datasets of real data (HMF or GiaB)
    # generate_labels()

    print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))


if __name__ == '__main__':
    main()
