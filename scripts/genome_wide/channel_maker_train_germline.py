import os, errno
import pysam
import numpy as np
import argparse
import bz2file
import gzip
import pickle
from time import time
import logging
import statistics
from collections import Counter

HPC_MODE = True


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


def read_BED(SVmode, chrName):
    '''
    
    :return: For INDELs: 
             start_SV_DEL: list of start positions for deletions.
             end_SV_DEL: list of end positions for deletions. It corresponds to the start positions.
             start_SV_INS: list of start positions for insertions.
             
             For INV and DUP:
             start_SV: list of start positions for inversion/duplication.
             end_SV: list of end positions for inversion/duplication.
             
    '''''

    if SVmode == 'INDEL':

        # Path on the HPC
        if HPC_MODE:
            wd = "/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/WG/run_"+ SVmode \
                 + "_500K/genomes/"
            truth_file = wd + "SV/GRCh37_"+SVmode+".bed"
        else:
            wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/"
            truth_file = wd + "SV/chr17_"+SVmode+".bed"

        assert os.path.isfile(truth_file)

        chr_list_DEL = []
        start_SV_DEL = []
        end_SV_DEL = []

        chr_list_INS = []
        start_SV_INS = []

        with(open(truth_file, 'r')) as i:
            for line in i:
                line = line.rstrip()
                columns = line.split("\t")
                chrom = str(columns[0])
                if columns[4] == "DEL" and chrom == chrName:
                    start_SV_DEL.append(int(columns[1]))
                    end_SV_DEL.append(int(columns[3]))
                    chr_list_DEL.append(chrom)
                elif str(columns[4]) == "INS" and chrom == chrName:
                    start_SV_INS.append(int(columns[1]))
                    chr_list_INS.append(chrom)

        return start_SV_DEL, end_SV_DEL, start_SV_INS

    elif SVmode == 'INV' or SVmode == 'DUP':

        # Path on the HPC
        if HPC_MODE:
            wd = "/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_"+SVmode+"/genomes/SV/"
            truth_file = wd + "chr17_"+SVmode+".bed"

        chr_list = []
        start_SV = []
        end_SV = []

        with(open(truth_file, 'r')) as i:
            for line in i:
                line = line.rstrip()
                columns = line.split("\t")
                chrom = str(columns[0])
                if ( columns[4] == "INV" or columns[4] == "DUP" ) and chrom == chrName:
                    start_SV.append(int(columns[1]))
                    end_SV.append(int(columns[3]))
                    chr_list.append(chrom)

        return start_SV, end_SV


def channel_maker(ibam, chrName, sampleName, trainingMode, SVmode, outFile):
    # Prefix for the relative path
    workdir = 'Training/'
    if not HPC_MODE:
        workdir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/' + workdir

    # List used to store the channel vstacks
    ch_list = []

    # Check for BAM file existence
    assert os.path.isfile(ibam)
    # Open BAM file
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

    # Consider a single sample
    sample_list = [sampleName]

    # File with clipped read positions, output of the clipped_read_pos script
    clipped_read_pos_file = workdir + sampleName + \
        '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
    # File with the clipped read distances, output of the clipped_read_distance script
    clipped_read_distance_file = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
    # File with the clipped reads, output of the clipped_reads script
    clipped_reads_file = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
    # File with the coverage array, output of the coverage script
    coverage_file = 'coverage/' + chrName + '_coverage.npy.bz2'
    # File with the split reads and split read distance, output of the split_read_distance script
    split_read_distance_file = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

    # Check existence of files
    assert os.path.isfile(clipped_read_pos_file)

    for sample in sample_list:
        assert os.path.isfile(workdir + sample + '/' + clipped_read_distance_file)
        assert os.path.isfile(workdir + sample + '/' + clipped_reads_file)
        assert os.path.isfile(workdir + sample + '/' + coverage_file)
        assert os.path.isfile(workdir + sample + '/' + split_read_distance_file)

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
    clipped_reads_inversion = dict()
    clipped_reads_duplication = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    for sample in sample_list:
        prefix = workdir + sample + '/'

        logging.info('Considering %s' % sample)
        logging.info('Reading clipped read distances')
        with bz2file.BZ2File(prefix + clipped_read_distance_file, 'rb') as f:
            clipped_read_distance[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading clipped reads')
        with bz2file.BZ2File(prefix + clipped_reads_file, 'rb') as f:
            clipped_reads[sample], clipped_reads_inversion[sample], clipped_reads_duplication[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading coverage')
        with bz2file.BZ2File(prefix + coverage_file, 'rb') as f:
            coverage[sample] = np.load(file=f)
        logging.info('End of reading, coverage length: %d out of %d' % (len(coverage[sample]), chrLen))

        logging.info('Reading split read distances')
        with bz2file.BZ2File(prefix + split_read_distance_file, 'rb') as f:
            split_read_distance[sample], split_reads[sample] = pickle.load(f)
        logging.info('End of reading')

    clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]
    print('Number of clipped read positions:%d' % len(clipped_pos))

    if SVmode == 'INDEL':

        # load the position for the artifically generated deletions (DEL, start and end) and insertion (only start)
        start_SV_DEL, end_SV_DEL, start_SV_INS = read_BED(SVmode, chrName)
        pos_list = start_SV_DEL + end_SV_DEL + start_SV_INS

        logging.info('Number of BPJ positions: %d' % len(pos_list))

        label_list = ['DEL_start'] * len(start_SV_DEL) + \
                     ['DEL_end'] * len(end_SV_DEL) + \
                     ['INS_pos'] * len(start_SV_INS)

    elif SVmode == 'INV' or SVmode == 'DUP':

        # load the position for the artifically generated deletions (DEL, start and end) and insertion (only start)
        start_SV, end_SV = read_BED(SVmode, chrName)
        pos_list = start_SV + end_SV

        logging.info('Number of BPJ positions: %d' % len(pos_list))

        label_list = [SVmode + '_start'] * len(start_SV) + \
                     [SVmode + '_end'] * len(end_SV)

    label = []
    label_BPJ = []
    distance_BPJ = []

    # Dictionaries where to store the channel arrays as generated from the dictionaries
    clipped_read_distance_array = dict()
    clipped_read_distance_num = dict()
    clipped_read_distance_median = dict()

    clipped_reads_array = dict()
    clipped_reads_inversion_array = dict()
    clipped_reads_duplication_array = dict()

    coverage_array = dict()

    split_read_distance_array = dict()
    split_read_distance_num = dict()
    split_read_distance_median = dict()

    split_reads_array = dict()

    # Log info every n_r times
    n_r = 10 ** 3
    # print(n_r)
    last_t = time()

    # Iterate over the SV positions
    for i, outzipped in enumerate(zip(pos_list, label_list), start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d candidate positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        center_pos = outzipped[0]
        lab = outzipped[1]

        # Check if center_pos is within the chromosome boundaries
        if win_hlen <= center_pos <= (chrLen - win_hlen):

            # Consider all the clipped read positions that are within win_hlen from the breakpoint junction (center_pos)
            window_cr_pos = [p for p in clipped_pos if (center_pos - win_hlen) <= p <= (center_pos + win_hlen)]

            # print('Number of clipped read positions in window %d: %d' % (center_pos, len(window_cr_pos)))

            for pos in window_cr_pos:

                start_win = pos - win_hlen
                end_win = pos + win_hlen

                # Build arrays for the numpy vstack
                for sample in sample_list:

                    # clipped read distance
                    clipped_read_distance_array[sample] = dict()
                    clipped_read_distance_num[sample] = dict()
                    clipped_read_distance_median[sample] = dict()

                    for direction in ['forward', 'reverse']:
                        clipped_read_distance_array[sample][direction] = dict()
                        clipped_read_distance_num[sample][direction] = dict()
                        clipped_read_distance_median[sample][direction] = dict()

                    for direction in ['forward', 'reverse']:
                        # for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                        for clipped_arrangement in ['left', 'right']:
                            clipped_read_distance_array[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                           dtype=int)
                            clipped_read_distance_num[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                         dtype=int)
                            clipped_read_distance_median[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                            dtype=int)
                            for pos in range(start_win, end_win):
                                if pos in clipped_read_distance[sample][direction][clipped_arrangement].keys():
                                    clipped_read_distance_array[sample][direction][clipped_arrangement][
                                        pos - start_win] = \
                                        sum(clipped_read_distance[sample][direction][clipped_arrangement][pos])
                                    clipped_read_distance_num[sample][direction][clipped_arrangement][
                                        pos - start_win] = \
                                        len(clipped_read_distance[sample][direction][clipped_arrangement][pos])
                                    clipped_read_distance_median[sample][direction][clipped_arrangement][
                                        pos - start_win] = \
                                        statistics.median(
                                            clipped_read_distance[sample][direction][clipped_arrangement][pos])
                            # print(clipped_read_distance_array[direction][clipped_arrangement])

                    # clipped reads
                    clipped_reads_array[sample] = dict()
                    for split_direction in ['left', 'right']:
                        clipped_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in clipped_reads[sample][split_direction].keys():
                                clipped_reads_array[sample][split_direction][pos - start_win] = \
                                    clipped_reads[sample][split_direction][pos]

                    # clipped reads inversions
                    clipped_reads_inversion_array[sample] = dict()
                    for mate_position in ['before', 'after']:
                        clipped_reads_inversion_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in clipped_reads_inversion[sample][mate_position].keys():
                                clipped_reads_inversion_array[sample][mate_position][pos - start_win] = \
                                    clipped_reads_inversion[sample][mate_position][pos]

                    # clipped reads duplication
                    clipped_reads_duplication_array[sample] = dict()
                    for mate_position in ['before', 'after']:
                        clipped_reads_duplication_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in clipped_reads_duplication[sample][mate_position].keys():
                                clipped_reads_duplication_array[sample][mate_position][pos - start_win] = \
                                    clipped_reads_duplication[sample][mate_position][pos]

                    # coverage
                    coverage_array[sample] = coverage[sample][start_win:end_win]
                    assert len(coverage_array[sample]) == win_len

                    # split read distance
                    split_read_distance_array[sample] = dict()
                    split_read_distance_num[sample] = dict()
                    split_read_distance_median[sample] = dict()

                    for split_direction in ['left', 'right']:
                        split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                        split_read_distance_num[sample][split_direction] = np.zeros(win_len, dtype=int)
                        split_read_distance_median[sample][split_direction] = np.zeros(win_len, dtype=int)

                        if pos in split_read_distance[sample][split_direction].keys():
                            split_read_distance_array[sample][split_direction][pos - start_win] = \
                                sum(split_read_distance[sample][split_direction][pos])
                            split_read_distance_array[sample][split_direction][pos - start_win] = \
                                len(split_read_distance[sample][split_direction][pos])
                            split_read_distance_array[sample][split_direction][pos - start_win] = \
                                statistics.median(split_read_distance[sample][split_direction][pos])

                    # split reads
                    split_reads_array[sample] = dict()
                    for split_direction in ['left', 'right']:
                        split_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in split_reads[sample][split_direction].keys():
                                split_reads_array[sample][split_direction][pos - start_win] = \
                                    split_reads[sample][split_direction][pos]

                # Fill the numpy vstack

                vstack_list = []
                for sample in sample_list:
                    # logging.info("Considering sample %s" % sample)

                    vstack_list.append(coverage_array[sample])

                    for clipped_arrangement in ['left', 'right']:
                        vstack_list.append(clipped_reads_array[sample][clipped_arrangement])

                    for mate_position in ['before', 'after']:
                        vstack_list.append(clipped_reads_inversion_array[sample][mate_position])
                    for mate_position in ['before', 'after']:
                        vstack_list.append(clipped_reads_duplication_array[sample][mate_position])

                    for direction in ['forward', 'reverse']:
                        for clipped_arrangement in ['left', 'right']:
                            vstack_list.append(clipped_read_distance_array[sample][direction][clipped_arrangement])
                            vstack_list.append(clipped_read_distance_num[sample][direction][clipped_arrangement])
                            vstack_list.append(clipped_read_distance_median[sample][direction][clipped_arrangement])
                    for direction in ['left', 'right']:
                        vstack_list.append(split_reads_array[sample][direction])
                    for direction in ['left', 'right']:
                        vstack_list.append(split_read_distance_array[sample][direction])
                        vstack_list.append(split_read_distance_num[sample][direction])
                        vstack_list.append(split_read_distance_median[sample][direction])

                # logging.info("Shape of channel matrix: %s" % str(ch_vstack.shape))
                ch_vstack = np.vstack(vstack_list)
                ch_list.append(ch_vstack)
                # print(len(ch_list))

                # print('Vstack:')
                # for d in np.arange(ch_vstack.shape[0]):
                #    print(ch_vstack[d])

                # Append the label to the list of labels
                label.append(lab)
                # Is the position a breakpoint junction?
                if pos == center_pos:
                    # print('Center!')
                    label_BPJ.append(True)
                else:
                    label_BPJ.append(False)
                # Record the distance of the position from the breakpoint junction for postprocessing
                distance_BPJ.append(abs(pos - center_pos))

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    outDir = os.path.dirname(outFile)
    labelDir = outDir +'/label/'
    create_dir(labelDir)

    # print(Counter(label))
    # Save the list of labels
    with gzip.GzipFile(labelDir + sampleName + '_' + chrName + '_label.npy.gz', "w") as f:
        np.save(file=f, arr=label)
    f.close()
    # Save the list of flags for breakpoint junctions (True/False)
    with gzip.GzipFile(labelDir + sampleName + '_' + chrName + '_label_BPJ.npy.gz', "w") as f:
        np.save(file=f, arr=label_BPJ)
    f.close()
    # Save the list of distances from positions to breakpoint junctions
    with gzip.GzipFile(labelDir + sampleName + '_' + chrName + '_distance_BPJ.npy.gz', "w") as f:
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

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='germline.npy.gz',
                        help="Specify output")
    parser.add_argument('-s', '--sample', type=str, default='germline',
                        help="Specify sample")
    parser.add_argument('-m', '--svmode', type=str, default='INDEL',
                        help="Specify SV type")
    parser.add_argument('-t', '--train', type=bool, default=True,
                        help="Specify if training mode is active")
    parser.add_argument('-l', '--logfile', default='channel_maker_train.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    channel_maker(ibam=args.bam, chrName=args.chr, sampleName=args.sample, SVmode=args.svmode,
                  trainingMode=args.train, outFile=args.out)

    print('Elapsed time channel_maker_train on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))


if __name__ == '__main__':
    main()
