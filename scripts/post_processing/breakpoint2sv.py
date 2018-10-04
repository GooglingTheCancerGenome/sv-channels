'''
This script takes in input a BED file with candidate breakpoint positions and
returns a VCF file of SVs where the breakpoints are connected using paired end information
and information on clipped read positions.
'''
import itertools
import pysam
import logging
import os
from pathlib import Path
import argparse as ap
from collections import defaultdict, Counter
from intervaltree import Interval, IntervalTree
import numpy as np
from multiprocessing import Process, Pool



__authors__ = ["Luca Santuari", "Tilman Schäfers"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.0.1"
__status__ = "alpha"

###########
parser = ap.ArgumentParser(description='Provide breakpoint2sv arguments.')
parser.add_argument('--BAM', type=str, nargs='?', dest='bam_file', default='/Users/tschafers/Test_data/CNN/BAM/G1_dedup.bam')
parser.add_argument('--BED', type=str, nargs='?', dest='bed_file', default='/Users/tschafers/Test_data/CNN/SV/chr17B_T.proper.bed')
parser.add_argument('--OUT_DIR', type=str, nargs='?', dest='out_dir', default='/Users/tschafers/Test_data/CNN/Results/')
parser.add_argument('--MAX_READ_COUNT', type=int, nargs='?', dest='max_read_count', default=5000)
parser.add_argument('--MIN_MAPQ', type=int, nargs='?', dest='min_mapq', default=20)
parser.add_argument('--WIN_H_LEN', type=int, nargs='?', dest='win_h_len', default=250)

##################################
args = parser.parse_args()
# Window half length
win_hlen = args.win_h_len
# Window size
win_len = win_hlen * 2
# Minimum read mapping quality

'''
Generic functions used in the channel scripts
'''


# Return if a read is clipped on the left
def is_left_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the left, False otherwise
    '''
    if read.cigartuples is not None:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the right, False otherwise
    '''
    if read.cigartuples is not None:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right or on the left
def is_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the left or on the right, False otherwise
    '''
    if read.cigartuples is not None:
        if is_left_clipped(read) or is_right_clipped(read):
            return True
    return False


def read_breakpoints(bed_file):
    print('Reading BED file for breakpoints')
    assert os.path.isfile(bed_file)
    breakpoints = defaultdict(list)
    with(open(bed_file, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom = str(columns[0])
            if columns[3][:3] == "DEL":
                breakpoints[chrom].append(int(columns[1]))
                breakpoints[chrom].append(int(columns[2]))
    for chr in breakpoints.keys():
        breakpoints[chr] = sorted(breakpoints[chr])
    # print(breakpoints)
    return breakpoints

##Accepts a list of arguments(breakpoints,chr)##
def breakpoint_to_sv(pargs):
    ###Extract arguments
    breakpoints = pargs[1]
    chr = pargs[0]
    ##Logging
    basename = os.path.splitext(os.path.basename(args.bed_file))[0]
    logging.basicConfig(filename=args.out_dir+basename+'_'+chr+'.log',level=logging.DEBUG, 
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    

    #open BAM file
    logging.info('Reading bam file: '+ args.bam_file)
    assert os.path.isfile(args.bam_file)
    aln = pysam.AlignmentFile(args.bam_file, "rb")
    # Check if the BAM file in input exists
    logging.info('Createing IntervalTree...')
    chr_tree = defaultdict(IntervalTree)
    # Create interval windows around candidate breakpoint positions
    for pos in breakpoints[chr]:
        chr_tree[chr][pos - win_hlen: pos + win_hlen + 1] = int(pos)
    logging.info('IntervalTree created.')
    links = []
    no_cr_pos = []
    npos = 1
    scanned_reads = set()
    logging.info('Processing Chr %s' % str(chr))
    for pos in breakpoints[chr]:
        if not npos % 1000:
            logging.info('Processed %d positions...' % npos)
        # print('Pos: %d' % pos)
        start = pos - win_hlen
        end = pos + win_hlen + 1
        right_clipped_array = np.zeros(win_len)
        left_clipped_array = np.zeros(win_len)
        count_reads = aln.count(chr, start, end)
        # Fetch the reads mapped on the chromosome
        if count_reads <= args.max_read_count:
            logging.info('Fetching %d reads in region %s:%d-%d' % (count_reads, chr, start, end))
            for read in aln.fetch(chr, start, end, multiple_iterators=True):
                # Both read and mate should be mapped
                if not read.is_unmapped and not read.mate_is_unmapped and \
                        read.mapping_quality >= args.min_mapq and \
                        read.reference_name == read.next_reference_name:
                    # Filled vectors with counts of clipped reads at clipped read positions
                    if is_right_clipped(read):
                        cpos = read.reference_end + 1
                        if start <= cpos <= end:
                            right_clipped_array[cpos - start - 2] += 1
                    if is_left_clipped(read):
                        cpos = read.reference_start
                        if start <= cpos <= end:
                            left_clipped_array[cpos - start - 1] += 1
                    if read.query_name not in scanned_reads:
                        match = chr_tree[read.next_reference_name][read.next_reference_start]
                        if chr != read.next_reference_name or start > read.next_reference_start or \
                                end < read.next_reference_start:
                            if match:
                                for m in match:
                                    int_start, int_end, int_data = m
                                    # print('%s -> %s' % (pos, int_data))
                                    # Pay attention of double insertions. The same read pair will be added
                                    # from both intervals, leading to double count.
                                    links.append(
                                        frozenset({chr + '_' + str(pos),
                                                    read.next_reference_name + '_' + str(int_data)}))
                                    scanned_reads.add(read.query_name)

        # print('Right clipped:\n%s' % right_clipped_array)
        # print('Left clipped:\n%s' % left_clipped_array)
        if sum(right_clipped_array) + sum(left_clipped_array) == 0:
            # print('Pos %d has no CR pos' % pos)
            no_cr_pos.append(chr + '_' + str(pos))
        else:
            if max(right_clipped_array) > max(left_clipped_array):
                max_i = np.where(right_clipped_array == max(right_clipped_array))
                # print('Right: %d -> %s' % (pos, max_i[0]))
            else:
                max_i = np.where(left_clipped_array == max(left_clipped_array))
                # print('Left: %d -> %s' % (pos, max_i[0]))

        npos += 1

    links_counts = Counter(links)
    logging.info('Set size: %d' % len(links_counts))
    logging.info('No CR pos: %d' % len(no_cr_pos))
    logging.info('Connections with min 3 read pairs: %d' % len([v for l, v in links_counts.items() if v > 2]))

    i = 0
    while len([v for l, v in links_counts.items() if v > i]) > 5000:
        i += 1
    logging.info('%d connections with min %d RP' % (len([v for l, v in links_counts.items() if v > i]), i))
    # Return link positions, and counts
    vcf_out = args.out_dir+basename+'_'+chr+'.vcf'
    linksToVcf(links_counts, vcf_out)
    logging.info('Finished breakoint assembly for chr%s ' % (chr))

def linksToVcf(links_counts, filename):
    cols = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    # Write VCF Header
    with open(filename, 'w') as sv_calls:
        sv_calls.write('##fileformat=VCFv4.2\n')
        sv_calls.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        sv_calls.write('##fileDate=20180726\n')
        sv_calls.write('##reference=GATK-GRCh-hg19')
        sv_calls.write('##ALT=<ID=DEL,Description="Deletion">\n')
        sv_calls.write('##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
        sv_calls.write('##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n')
        sv_calls.write('##FORMAT=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        sv_calls.write(cols)
        for l, v in links_counts.items():
            interval = list(l)
            s1 = interval[0].split('_')
            s2 = interval[1].split('_')
            chr = s1[0]

            s1[1] = int(s1[1])
            s2[1] = int(s2[1])

            if s1[1] < s2[1]:
                # print('%d < %d' % (s1[1], s2[1]))
                start = s1[1]
                stop = s2[1]
            else:
                # print('%d > %d' % (s1[1], s2[1]))
                start = s2[1]
                stop = s1[1]
            f_line = "SVTYPE:PE:END"
            s_line = 'SVTYPE=%s;PE=%s;END=%s' % ('DEL', v, stop)
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chr, start, '.', 'N','<DEL>', '.', 'PASS', '.', f_line, s_line)
            # print(line)
            sv_calls.write(line + '\n')
        print("VCF file written!")

def main():
    breakpoints = read_breakpoints(args.bed_file)
    ##Spawn processes for each chromosome
    P = Pool(processes=len(breakpoints.keys()))
    pargs = zip(breakpoints.keys(), itertools.repeat(breakpoints))
    P.map(breakpoint_to_sv, pargs)
    P.close()
    P.join()

   




       
    




if __name__ == '__main__':
    main()
