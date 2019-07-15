
# Imports

import argparse
import pysam
import json
import gzip
from time import time
from functions import *
import logging
from collections import defaultdict
import os


def get_clipped_read_distance(ibam, chrName, outFile):

    '''

    :param ibam: BAM file in input
    :param chrName: chromosome to consider
    :param outFile: output file where to store the clipped_read_distance dictionary
    :return:
    '''

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    # minimum read mapping quality to consider
    config = get_config_file()
    minMAPQ = config["DEFAULT"]["MIN_MAPQ"]

    # open BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # get chromosome length from BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Dictionary with positions as keys and list of clipped read distances as values
    clipped_read_distance = dict()
    # For clipped reads mapped in the forward and reverse orientation
    for direction in ['forward', 'reverse']:
        clipped_read_distance[direction] = dict()
        # For left- and right-clipped reads
        for clipped_arrangement in ['left', 'right', 'all']:
            clipped_read_distance[direction][clipped_arrangement] = defaultdict(list)

    def set_distance(direction, read, dist):
        '''

        :param direction: forward/reverse read direction
        :param read: read object of the class pysam.AlignedSegment
        :param dist: read to mate distance
        :return: None. Adds dist to the list of distances at a clipped read position for a certain read direction
        '''

        if direction == 'forward':
            pos = read.reference_end + 1
        elif direction == 'reverse':
            pos = read.reference_start
        clipped_read_distance[direction]['all'][pos].append(dist)

        if is_left_clipped(read):
            pos = read.reference_start
            #if pos not in clipped_read_distance[direction]['left'].keys():
            #    clipped_read_distance[direction]['left'][pos] = [dist]
            #else:
            clipped_read_distance[direction]['left'][pos].append(dist)
        elif is_right_clipped(read):
            pos = read.reference_end + 1
            #if pos not in clipped_read_distance[direction]['right'].keys():
            #    clipped_read_distance[direction]['right'][pos] = [dist]
            #else:
            clipped_read_distance[direction]['right'][pos].append(dist)

    # Consider all the chromosome: interval [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # Fetch the reads mapped on the chromosome
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

        # Both read and mate should be mapped
        if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ:
            # Read and mate should be mapped on the same chromosome
            if read.reference_name == read.next_reference_name:
                # Calculate absolute read to mate distance
                dist = abs(read.reference_start - read.next_reference_start)
                # Read is mapped in forward orientation, mate is in reverse orientation, read is mapped before mate
                if not read.is_reverse and read.mate_is_reverse and read.reference_start <= read.next_reference_start:
                    set_distance('forward', read, dist)
                # Read is mapped in reverse orientation, mate is in forward orientation, read is mapped after mate
                elif read.is_reverse and not read.mate_is_reverse and read.reference_start > read.next_reference_start:
                    set_distance('reverse', read, dist)

    # Write clipped read distance dictionaries
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(clipped_read_distance).encode('utf-8'))

    # to load it:
    # with gzip.GzipFile(outFile, 'r') as fin:
    #     clipped_read_distance = json.loads(fin.read().decode('utf-8'))


def main():

    # Default BAM file for testing
    # On the HPC
    #wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    #inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create channels with distance between clipped/non-clipped reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='clipped_read_distance.json.gz',
                        help="Specify output")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-l', '--logfile', default='clipped_read_distance.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    cmd_name = 'clipped_read_distance'
    output_dir = os.path.join(args.outputpath, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, '_'.join((args.chr, args.logfile)))
    output_file = os.path.join(output_dir, '_'.join((args.chr, args.out)))

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    get_clipped_read_distance(ibam=args.bam, chrName=args.chr, outFile=output_file)
    logging.info('Time: clipped read distance on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
