# Imports
import argparse
import pysam
import bz2file
import pickle
from time import time
import logging
from functions import *
from collections import defaultdict


def get_split_read_distance(ibam, chrName, outFile):

    # minimum mapping quality to consider
    minMAPQ = 30

    # Dictionary to store left/right split read distance
    split_read_distance = dict()
    for split_direction in ['left', 'right']:
        split_read_distance[split_direction] = defaultdict(list)

    # Dictionary to store number of left/right split reads
    split_reads = dict()
    for split_direction in ['left', 'right']:
        split_reads[split_direction] = defaultdict(int)

    # Load BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Get chromosome length from BAM file
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    # print(chrLen)

    # Fetch reads mapped on the chromosome
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    # Log information every n_r reads
    n_r = 10 ** 6
    # print(n_r)
    last_t = time()
    # print(type(last_t))
    for i, read in enumerate(iter, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Read and mate should be mapped, read should have a minimum mapping quality
        if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ:
            # Check if the read has a supplementary alignment: is the read a split read?
            if read.has_tag('SA'):
                # The read is left clipped
                # if is_left_clipped(read):
                #     # print('Left clipped')
                #     # print(read)
                #     chr, pos, strand = get_suppl_aln(read)
                #     # print('%s:%d-%s' % (chr, pos, strand))
                #     # The read and the supplementary alignment are on the same chromosome
                #     if chr == read.reference_name:
                #         # print('Left split')
                #         # print(str(read))
                #         refpos = read.reference_start
                #         #if pos not in split_read_distance['left'].keys():
                #         #    split_read_distance['left'][pos] = [abs(refpos - pos)]
                #         #else:
                #         split_read_distance['left'][pos].append(abs(refpos - pos))
                #         #if pos not in split_reads['left'].keys():
                #         #    split_reads['left'][pos] = 1
                #         #else:
                #         # print('Adding for left at %d position' % refpos)
                #         # print('Adding for right at %d position' % pos)
                #         split_reads['left'][refpos] += 1
                #         # split_reads['right'][pos] += 1

                # The read is right clipped
                if is_right_clipped(read):
                    # print('Right clipped')
                    # print(read)
                    chr, pos, strand = get_suppl_aln(read)
                    # print('%s:%d-%s' % (chr, pos, strand))
                    # The read and the supplementary alignment are on the same chromosome
                    if chr == read.reference_name:
                        print('Right split')
                        print(str(read))
                        refpos = read.reference_end
                        #if pos not in split_read_distance['right'].keys():
                        #    split_read_distance['right'][pos] = [abs(pos - refpos)]
                        #else:
                        split_read_distance['right'][pos].append(abs(pos - refpos))
                        #if pos not in split_reads['right'].keys():
                        #    split_reads['right'][pos] = 1
                        #else:
                        print('Adding for left at %d position' % pos)
                        print('Adding for right at %d position' % refpos)
                        split_reads['right'][refpos] += 1
                        split_reads['left'][pos] += 1


    # Save two dictionaries: split_read_distance and split_reads
    with bz2file.BZ2File(outFile, 'w') as f:
        pickle.dump((split_read_distance, split_reads), f)


def main():

    # Default BAM file for testing
    # On the HPC
    #wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    #inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create channels with split read distance for left/right split reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='split_read_distance.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='split_read_distance.log',
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
    get_split_read_distance(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()