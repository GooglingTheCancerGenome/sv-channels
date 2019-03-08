# Imports
import argparse
import logging
import os
import pickle
import bz2file
from collections import Counter
from time import time
import pysam
import functions as fun


def get_clipped_read_positions(ibam, chrName, outFile):
    '''
    
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the dictionary of clipped read positions
    :return: None. Outputs a dictionary with the positions of clipped read positions as keys and
    the number of clipped reads per position as values
    '''''

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    # Minimum read mapping quality to consider
    minMAPQ = 30

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # List to store the clipped read positions
    clipped_pos = []

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # Pysam iterator to fetch the reads
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    # Print every n_r alignments processed
    n_r = 10 ** 6
    # Record the current time
    last_t = time()

    for i, read in enumerate(iter, start=1):

        # Every n_r alignments, write log informations
        if not i % n_r:
            # Record the current time
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Both read and mate should be mapped, read should have a minimum mapping quality
        # if (not read.is_unmapped) and (not read.mate_is_unmapped) and read.mapping_quality >= minMAPQ:
        if (not read.is_unmapped) and read.mapping_quality >= minMAPQ:

            # if fun.has_indels(read):
            #     # print(read)
            #     dels_start, dels_end, ins = fun.get_indels(read)
            #     dels = dels_start + dels_end + ins
            #     clipped_pos.extend(dels)

            if fun.is_left_clipped(read):
                # read.reference_start is the 1-based start position of the read mapped on the reference genome
                cpos = read.reference_start
                clipped_pos.append(cpos)
            if fun.is_right_clipped(read):
                # read.reference_end is the 0-based end position of the read mapped on the reference genome
                cpos = read.reference_end + 1
                clipped_pos.append(cpos)

    # Close the BAM file
    bamfile.close()

    # Count the number of clipped reads per position
    clipped_pos_cnt = Counter(clipped_pos)

    logging.info('Number of unique positions: %d' % len(clipped_pos_cnt))

    # Write the output in pickle format
    with bz2file.BZ2File(outFile, 'wb') as f:
        pickle.dump(clipped_pos_cnt, f)


def main():

    # Default BAM file for testing
    # On the HPC
    #wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    #inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"
    #wd = '/Users/lsantuari/Documents/mount_points/hpc_mnt/Datasets/CretuStancu2017/Patient1/'
    #inputBAM = wd + 'Patient1.bam'

    # Default chromosome is 17 for the artificial data

    # Parse the arguments of the script
    parser = argparse.ArgumentParser(description='Get clipped reads positions')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='clipped_read_pos.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='clipped_read_pos.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    get_clipped_read_positions(ibam=args.bam, chrName=args.chr, outFile=args.out)
    logging.info('Time: clipped read positions on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
