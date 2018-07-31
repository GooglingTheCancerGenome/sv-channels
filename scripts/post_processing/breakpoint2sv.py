'''
This script takes in input a BED file with candidate breakpoint positions and
returns a VCF file of SVs where the breakpoints are connected using paired end information
and information on clipped read positions.
'''

import pysam
import os
from collections import defaultdict, Counter
from intervaltree import Interval, IntervalTree
import numpy as np

__authors__ = ["Luca Santuari"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.0.1"
__status__ = "alpha"

# parameters

# Locally:
work_dir = '/Users/tschafers/CNN/scripts/post_processing/Test/'
bed_file = work_dir + 'genomes/SV/chr17B_T.proper_small.bed'
bam_file = work_dir + 'samples/G1/BAM/G1/mapping/G1_dedup.bam'

# On HPC:
#work_dir = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/'
#bed_file = 'genomes/SV/chr17B_T.proper.bed'
#bam_file = 'samples/G1/BAM/G1/mapping/G1_dedup.bam'

# Window half length
win_hlen = 250
# Window size
win_len = win_hlen * 2

# Minimum read mapping quality
min_mapq = 20


class Location:

    def __init__(self, chrA, posA_start, posA_end,
                 chrB, posB_start, posB_end):

        if chrA == chrB and posA_start <= posB_start or chrA < chrB:

            self.chr1 = chrA
            self.pos1_start = posA_start
            self.pos1_end = posA_end

            self.chr2 = chrB
            self.pos2_start = posB_start
            self.pos2_end = posB_end

        elif chrA > chrB:

            self.chr1 = chrB
            self.pos1_start = posB_start
            self.pos1_end = posB_end

            self.chr2 = chrA
            self.pos2_start = posA_start
            self.pos2_end = posA_end


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


def breakpoint_to_sv():
    breakpoints = read_breakpoints(bed_file)
    # Check if the BAM file in input exists
    assert os.path.isfile(bam_file)
    # open the BAM file
    aln = pysam.AlignmentFile(bam_file, "rb")

    chr_tree = defaultdict(IntervalTree)
    for chr in breakpoints.keys():
        #Create interval windows around candidate breakpoint positions
        for pos in breakpoints[chr]:
            chr_tree[chr][pos - win_hlen: pos + win_hlen + 1] = pos
    links = []
    no_cr_pos = []
    npos = 1
    scanned_reads = set()
    for chr in breakpoints.keys():
        print('Chr %s' % str(chr))
        for pos in breakpoints[chr]:
            if not npos % 1000:
                print('Processed %d positions...' % npos)
            # print('Pos: %d' % pos)
            start = pos - win_hlen
            end = pos + win_hlen + 1
            right_clipped_array = np.zeros(win_len)
            left_clipped_array = np.zeros(win_len)

            # Fetch the reads mapped on the chromosome
            for read in aln.fetch(chr, start, end):
                # Both read and mate should be mapped
                if not read.is_unmapped and not read.mate_is_unmapped and \
                        read.mapping_quality >= min_mapq:

                    # Filled vectors with counts of clipped reads at clipped read positions
                    if is_right_clipped(read):
                        cpos = read.reference_end + 1
                        if start <= cpos <= end:
                            right_clipped_array[cpos-start-2] += 1
                    if is_left_clipped(read):
                        cpos = read.reference_start
                        if start <= cpos <= end:
                            left_clipped_array[cpos-start-1] += 1

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


            #print('Right clipped:\n%s' % right_clipped_array)
            #print('Left clipped:\n%s' % left_clipped_array)
            if sum (right_clipped_array) + sum(left_clipped_array) == 0:
                # print('Pos %d has no CR pos' % pos)
                no_cr_pos.append(chr + '_' + str(pos))
            else:
                if max(right_clipped_array) > max(left_clipped_array):
                    max_i = np.where(right_clipped_array == max(right_clipped_array))
                    #print('Right: %d -> %s' % (pos, max_i[0]))
                else:
                    max_i = np.where(left_clipped_array == max(left_clipped_array))
                    #print('Left: %d -> %s' % (pos, max_i[0]))

            npos += 1

    links_counts = Counter(links)
    print('Set size: %d' % len(links_counts))
    print('No CR pos: %d' % len(no_cr_pos))
    print('Connections with min 3 read pairs: %d' % len([v for l, v in links_counts.items() if v > 2]))

    i = 0
    while len([v for l, v in links_counts.items() if v > i]) > 5000:
        i += 1
    print('%d connections with min %d RP' % (len([v for l, v in links_counts.items() if v > i]), i))
    #Return link positions, and counts
    return links_counts

def linksToVcf(links_counts):
    filename = 'sv_calls2.vcf'
    cols = '#CHROM\tPOS\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    #Write VCF Header
    with open(filename, 'a') as sv_calls:
        sv_calls.write('##fileformat=VCFv4.2\n')
        sv_calls.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        sv_calls.write('##fileDate=20180726\n')
        sv_calls.write('##ALT=<ID=DEL,Description="Deletion">\n')
        sv_calls.write('##INFO=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
        sv_calls.write('##INFO=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n')
        sv_calls.write('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')
        sv_calls.write(cols)
        for l,v in links_counts.items():
            interval = list(l)
            s1 = interval[0].split('_')
            s2 = interval[1].split('_')
            chr = s1[0]
            start = s1[1]
            stop = s2[1]
            f_line = 'SVTYPE=%s;PE=%s;END=%s' % ('DEL',v,stop)
            line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chr,start,'.','.','<DEL>','q30','PASS',f_line)
            print(line)
            sv_calls.write(line+'\n')




def main():
    links_counts = breakpoint_to_sv()
    linksToVcf(links_counts)


if __name__ == '__main__':
    main()