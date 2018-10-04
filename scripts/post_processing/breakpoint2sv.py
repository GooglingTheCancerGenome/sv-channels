'''
This script takes in input a BED file with candidate breakpoint positions and
returns a VCF file of SVs where the breakpoints are connected using paired end information
and information on clipped read positions.
'''
import itertools
import pysam
import os
from collections import defaultdict, Counter
from intervaltree import Interval, IntervalTree
import numpy as np
from multiprocessing import Process, Pool

__authors__ = ["Luca Santuari", "Tilman Schäfers"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.0.1"
__status__ = "alpha"

# parameters

HPC_MODE = False

if not HPC_MODE:

    # Locally:
    #work_dir = '/Users/tschafers/CNN/scripts/post_processing/Test/'
    #bed_file = work_dir + 'genomes/SV/chr17B_T.proper_small.bed'
    #bam_file = work_dir + 'samples/G1/BAM/G1/mapping/G1_dedup.bam'

    #work_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL'
    #bed_file = os.path.join(work_dir, 'SV/chr17B_T.proper.bed')
    #bam_file = os.path.join(work_dir, 'BAM/G1_dedup.bam')
    #vcf_output = os.path.join(work_dir, 'VCF/chr17B_T.vcf')
    work_dir = '/Users/tschafers/Test_data/CNN/'
    bed_file = os.path.join(work_dir, 'SV/chr17B_T.proper.bed')
    bam_file = os.path.join(work_dir, 'BAM/G1_dedup.bam')
    vcf_output = os.path.join(work_dir, 'VCF/chr17B_T_m.vcf')


else:

    # On HPC:
    patient_number = str(1)
    bed_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/060818/TestData_060818/PATIENT' + \
               patient_number + '_DEL.sorted.bed'
    bam_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient' + patient_number + '/Patient' + \
               patient_number + '.bam'
    vcf_output = '/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Test/060818/TestData_060818/PATIENT' + \
                 patient_number + '_DEL.vcf'

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
def breakpoint_to_sv(args):
    print("Starting")
    breakpoints = args[1]
    chr = args[0]
    assert os.path.isfile(bam_file)
    # open the BAM file
    aln = pysam.AlignmentFile(bam_file, "rb")
    # Check if the BAM file in input exists
    #print('Create IntervalTree...')
    chr_tree = defaultdict(IntervalTree)
    # Create interval windows around candidate breakpoint positions
    for pos in breakpoints[chr]:
        chr_tree[chr][pos - win_hlen: pos + win_hlen + 1] = int(pos)
    #print('IntervalTree created.')
    links = []
    no_cr_pos = []
    npos = 1
    scanned_reads = set()
    print('Processing Chr %s' % str(chr))
    for pos in breakpoints[chr]:
        if not npos % 1000:
            print('Processed %d positions...' % npos)
        # print('Pos: %d' % pos)
        start = pos - win_hlen
        end = pos + win_hlen + 1
        right_clipped_array = np.zeros(win_len)
        left_clipped_array = np.zeros(win_len)
        count_reads = aln.count(chr, start, end)
        # Fetch the reads mapped on the chromosome
        if count_reads <= 10000:
            # print('Fetching %d reads in region %s:%d-%d' % (count_reads, chr, start, end))
            for read in aln.fetch(chr, start, end, multiple_iterators=True):
                # Both read and mate should be mapped
                if not read.is_unmapped and not read.mate_is_unmapped and \
                        read.mapping_quality >= min_mapq and \
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
    print('Set size: %d' % len(links_counts))
    print('No CR pos: %d' % len(no_cr_pos))
    print('Connections with min 3 read pairs: %d' % len([v for l, v in links_counts.items() if v > 2]))

    i = 0
    while len([v for l, v in links_counts.items() if v > i]) > 5000:
        i += 1
    print('%d connections with min %d RP' % (len([v for l, v in links_counts.items() if v > i]), i))
    # Return link positions, and counts
    linksToVcf(links_counts)

def linksToVcf(links_counts):
    filename = vcf_output
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
    breakpoints = read_breakpoints(bed_file)
    ##Spawn processes for each chromosome
    P = Pool(processes=len(breakpoints.keys()))
    args = zip(breakpoints.keys(), itertools.repeat(breakpoints))
    P.map(breakpoint_to_sv, args)
    P.close()
    P.join()

   




       
    




if __name__ == '__main__':
    main()
