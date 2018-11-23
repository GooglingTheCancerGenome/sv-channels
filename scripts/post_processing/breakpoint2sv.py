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
import multiprocessing
import datetime
import time



__authors__ = ["Luca Santuari", "Tilman Schäfers"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.0.1"
__status__ = "alpha"

###########
parser = ap.ArgumentParser(description='Provide breakpoint2sv arguments.')
parser.add_argument('--BAM', type=str, nargs='?', dest='bam_file',
                    #default='/Users/tschafers/Test_data/CNN/BAM/G1_dedup.bam')
                    default='/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/'+\
                            'run_test_INDEL/BAM/G1_dedup.bam')
parser.add_argument('--BED', type=str, nargs='?', dest='bed_file',
                    #default='/Users/tschafers/Test_data/CNN/SV/chr17B_T.proper.bed')
                    default='/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/'+\
                            'run_test_INDEL/SV/chr17B_T.proper.bed')
parser.add_argument('--OUT_DIR', type=str, nargs='?', dest='out_dir',
                    #default='/Users/tschafers/Test_data/CNN/Results/')
                    default='/Users/lsantuari/Documents/Processed/Breakpoint2SV/Results/')
parser.add_argument('--MAX_READ_COUNT', type=int, nargs='?', dest='max_read_count', default=5000)
parser.add_argument('--MIN_MAPQ', type=int, nargs='?', dest='min_mapq', default=20)
parser.add_argument('--WIN_H_LEN', type=int, nargs='?', dest='win_h_len', default=250)
parser.add_argument('--VCF_OUT', type=str, nargs='?', dest='vcf_out', default='G1_deepsv_indels.vcf')

##################################
args = parser.parse_args()
# Window half length
win_hlen = args.win_h_len
# Window size
win_len = win_hlen * 2

strand = {False:'+', True:'-'}

bp_pos_dict = defaultdict(dict)
bp_counter_sum = []

class ExceptionWrapper(object):

    def __init__(self, ee):
        self.ee = ee
        __,  __, self.tb = sys.exc_info()

    def re_raise(self):
        raise self.ee.with_traceback(self.tb)
        # for Python 2 replace the previous line by:
        # raise self.ee, None, self.tb

'''
Generic functions used in the channel scripts
'''

def get_chr_len_dict(ibam):

    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header

    chr_len_dict = {}
    # print(header_dict['SQ'])
    for d in header_dict['SQ']:
        chr_len_dict[d['SN']] = d['LN']
    # print(chr_len_dict)
    return chr_len_dict


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


def read_breakpoints_single_position(bed_file):
    print('Reading BED file for breakpoints')
    assert os.path.isfile(bed_file)
    breakpoints = defaultdict(list)
    with(open(bed_file, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom = str(columns[0])
            breakpoints[chrom].append(int(columns[1]))
    for chr in breakpoints.keys():
        breakpoints[chr] = sorted(breakpoints[chr])
    # print(breakpoints)
    return breakpoints

##Accepts a list of arguments(breakpoints,chr)##
def breakpoint_to_sv(chr,breakpoints):
    # Read BAM file
    assert os.path.isfile(args.bam_file)
    aln = pysam.AlignmentFile(args.bam_file, "rb")
    ##Logging
    basename = os.path.splitext(os.path.basename(args.bed_file))[0]
    logging.basicConfig(filename=args.out_dir+basename+'_'+chr+'.log',level=logging.DEBUG, filemode='w',
                        format='%(asctime)s %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    #open BAM file
    logging.info('Reading bam file: '+ args.bam_file)
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
                                        frozenset({chr + '_' + str(pos) + '_' + strand[read.is_reverse],
                                                   read.next_reference_name + '_' +
                                                   str(int_data) + '_' + strand[read.mate_is_reverse]}))
                                    scanned_reads.add(read.query_name)

        # print('Right clipped:\n%s' % right_clipped_array)
        # print('Left clipped:\n%s' % left_clipped_array)
        if sum(right_clipped_array) + sum(left_clipped_array) == 0:
            # print('Pos %d has no CR pos' % pos)
            no_cr_pos.append(chr + '_' + str(pos))
        else:
            if max(right_clipped_array) > max(left_clipped_array):
                bp_pos_dict[chr][pos] = start + np.where(right_clipped_array == max(right_clipped_array))[0][0] + 1
                # print('Right: %d -> %s' % (pos, max_i[0]))
            else:
                bp_pos_dict[chr][pos] = start + np.where(left_clipped_array == max(left_clipped_array))[0][0] + 1
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
    logging.info('Finished breakpoint assembly for chr%s ' % (chr))
    return(links_counts)

def linksToVcf(links_counts, filename, ibam):
    cols = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    # Write VCF Header
    with open(filename, 'w') as sv_calls:

        now = datetime.datetime.now()

        sv_calls.write('##fileformat=VCFv4.2\n')
        sv_calls.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        sv_calls.write('##fileDate='+now.strftime("%Y-%m-%d %H:%M")+'\n')
        sv_calls.write('##fileTime=' + now.strftime("%H:%M") + '\n')
        sv_calls.write('##reference=GATK-GRCh-hg19')
        sv_calls.write('##ALT=<ID=DEL,Description="Deletion">\n')
        sv_calls.write('##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n')
        sv_calls.write('##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n')
        sv_calls.write('##FORMAT=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n')

        # writing contig info
        chr_dict = get_chr_len_dict(ibam)
        for k in chr_dict.keys():
            sv_calls.write('##contig=<ID='+k+',length='+str(chr_dict[k])+'>\n')

        sv_calls.write(cols)

        sv_evaluated = []

        for l, v in links_counts.items():

            # logging.info('Candidate link => %s' % (l))

            interval = list(l)
            s1 = interval[0].split('_')
            s2 = interval[1].split('_')

            chrA = s1[0]
            chrB = s2[0]

            s1[1] = int(s1[1])
            s2[1] = int(s2[1])

            if s1[1] in bp_pos_dict[chrA] and s2[1] in bp_pos_dict[chrB] and v > 2:

                posA = bp_pos_dict[chrA][s1[1]]
                posB = bp_pos_dict[chrB][s2[1]]

                strandA = s1[2]
                strandB = s2[2]

                fs = frozenset({str(chrA) + '_' + str(posA) + '_' + strandA,
                               str(chrB) + '_' + str(posB) + '_' + strandB})

                if fs not in sv_evaluated:
                    # logging.info('Writing link => %s:%d-%s:%d' % (chrA, posA, chrB, posB))

                    sv_evaluated.append(fs)

                    if posA <= posB:
                        # print('%d < %d' % (s1[1], s2[1]))
                        start = posA
                        stop = posB
                    else:
                        # print('%d > %d' % (s1[1], s2[1]))
                        start = posB
                        stop = posA

                    if chrA == chrB:
                        svtype = 'DEL'
                    else:
                        svtype = 'BND'

                    f_line = "SVTYPE:PE:END"
                    s_line = 'SVTYPE=%s;PE=%s;END=%s;STRAND=%s' % ('DEL', v, stop, strandA+strandB)
                    line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (chrA, start, '.', 'N',
                                                                       '<'+svtype+'>', '.', 'PASS',
                                                                       '.', f_line, s_line)
                    # print(line)
                    sv_calls.write(line + '\n')
                # else:
                    # logging.info('Link %s:%d-%s:%d considered already' % (chrA, posA, chrB, posB))
        print("VCF file written!")


def on_return(retval):
    bp_counter_sum.extend(retval)

def main():
    breakpoints = read_breakpoints(args.bed_file)
    ##Spawn processes for each chromosome
    print('Found chromosomes:')
    for k in breakpoints.keys(): print (k)
    print('Starting assembly')
    ###### Parallel execution ########
    start_time = time.time()  
    P = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    for chr in breakpoints.keys(): 
        P.apply_async(breakpoint_to_sv, args=(chr,breakpoints), callback=on_return)
    P.get()
    P.close()
    P.join()
    print('Finished breakpoint assembly')
    print("Writing intervals to VCF")
    temp = sum(bp_counter_sum,Counter())
    linksToVcf(temp, args.vcf_out, ibam = args.bam_file)
    print('Finished breakpoint assembly')
    print("--- %s seconds ---" % (time.time() - start_time))
    # mychr = '17'
    # breakpoint_to_sv([mychr, breakpoints])

    


if __name__ == '__main__':
    main()
