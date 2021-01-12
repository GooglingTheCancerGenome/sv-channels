import argparse
import gzip
import json
import logging
import os
from collections import defaultdict
from time import time

import pysam
from functions import *


def get_clipped_read_distance(ibam, chrName, min_mapq, outFile):
    '''
    :param ibam: BAM file in input
    :param chrName: chromosome to consider
    :param outFile: output file where to store the clipped_read_distance dictionary
    :return:
    '''
    bamfile = pysam.AlignmentFile(ibam, "rb")
    bam_mean, bam_stddev = get_insert_size(ibam, bamfile, min_mapq)
    # get chromosome length
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Dictionary with positions as keys and list of clipped read distances as values
    clipped_read_distance = dict()
    # For clipped reads mapped in the forward and reverse orientation
    for direction in ['forward', 'reverse']:
        clipped_read_distance[direction] = dict()
        # For left- and right-clipped reads
        for clipped_arrangement in ['left', 'right', 'all']:
            clipped_read_distance[direction][
                clipped_arrangement] = defaultdict(list)

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
            # if pos not in clipped_read_distance[direction]['left'].keys():
            #    clipped_read_distance[direction]['left'][pos] = [dist]
            # else:
            clipped_read_distance[direction]['left'][pos].append(dist)
        elif is_right_clipped(read):
            pos = read.reference_end + 1
            # if pos not in clipped_read_distance[direction]['right'].keys():
            #    clipped_read_distance[direction]['right'][pos] = [dist]
            # else:
            clipped_read_distance[direction]['right'][pos].append(dist)

    # Consider all the chromosome: interval [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # Log information every n_r reads
    n_r = 10 ** 6
    last_t = time()
    for i, read in enumerate(bamfile.fetch(chrName, start_pos, stop_pos), start=1):
        if not i % n_r:
            now_t = time()
            logging.info("%d alignments processed (%f alignments / s)" %
                         (i, n_r / (now_t - last_t)))
            last_t = time()

        # Both read and mate should be mapped
        if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= min_mapq:
            # Read and mate should be mapped on the same chromosome
            if read.reference_name == read.next_reference_name:
                # Calculate absolute read to mate distance
                dist = abs(read.reference_start - read.next_reference_start)
                dist = (dist - bam_mean) / bam_stddev
                # Read is mapped in forward orientation, mate is in reverse orientation, read is mapped before mate
                if not read.is_reverse and read.mate_is_reverse and read.reference_start <= read.next_reference_start:
                    set_distance('forward', read, dist)
                # Read is mapped in reverse orientation, mate is in forward orientation, read is mapped after mate
                elif read.is_reverse and not read.mate_is_reverse and read.reference_start > read.next_reference_start:
                    set_distance('reverse', read, dist)

    # Write clipped read distance dictionaries
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(clipped_read_distance).encode('utf-8'))


def main():
    parser = argparse.ArgumentParser(
        description='Create channels with distance between clipped/non-clipped reads')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chr',
                        type=str,
                        default='12',
                        help="Specify chromosome")
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        default='clipped_read_distance.json.gz',
                        help="Specify output")
    parser.add_argument('-p',
                        '--outputpath',
                        type=str,
                        default='.',
                        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='clipped_read_distance.log',
                        help='File in which to write logs.')
    parser.add_argument('-m',
                        '--min_mapq',
                        type=int,
                        default=10,
                        help='Minimum read mapping quality')
    args = parser.parse_args()
    cmd_name = 'clipped_read_distance'
    output_dir = os.path.join(args.outputpath, cmd_name)
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, '_'.join((args.chr, args.logfile)))
    output_file = os.path.join(output_dir, '_'.join((args.chr, args.out)))
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    get_clipped_read_distance(ibam=args.bam,
                              chrName=args.chr,
                              min_mapq=args.min_mapq,
                              outFile=output_file)
    logging.info('Time: clipped read distance on BAM %s and Chr %s: %f' %
                 (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
