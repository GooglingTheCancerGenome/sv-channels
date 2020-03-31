'''
This script takes in input a BED file with candidate breakpoint positions and
returns a VCF file of SVs where the breakpoints are connected using paired end information
and information on clipped read positions.
'''

import argparse as ap
import datetime
import itertools
import logging
import multiprocessing
import os
import time
from collections import Counter, defaultdict
from pathlib import Path
from subprocess import call

import numpy as np
import pysam
import twobitreader as twobit
from intervaltree import Interval, IntervalTree

from aux_functions import *

__authors__ = ["Luca Santuari", "Tilman Schäfers"]
__license__ = "Apache License, Version 2.0"
__version__ = "0.0.1"
__status__ = "alpha"

###########
parser = ap.ArgumentParser(description='Provide breakpoint2sv arguments.')
parser.add_argument('--BAM', type=str, nargs='?', dest='bam_file',
                    # default='/Users/tschafers/Test_data/CNN/BAM/G1_dedup.bam')
                    default='/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/' + \
                            'run_test_INDEL/BAM/G1_dedup.bam')
parser.add_argument('--BED', type=str, nargs='?', dest='bed_file',
                    # default='/Users/tschafers/Test_data/CNN/SV/bed/Patient1_94.bed')
                    default='/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/' + \
                            'run_test_INDEL/SV/chr17B_T.proper.bed')
parser.add_argument(
    '--OUT_DIR',
    type=str,
    nargs='?',
    dest='out_dir',
    #default=os.getcwd())
    default='/Users/lsantuari/Documents/Processed/Breakpoint2SV/Results/')
parser.add_argument('--MAX_READ_COUNT',
                    type=int,
                    nargs='?',
                    dest='max_read_count',
                    default=5000)
parser.add_argument('--MIN_MAPQ',
                    type=int,
                    nargs='?',
                    dest='min_mapq',
                    default=20)
parser.add_argument('--WIN_H_LEN',
                    type=int,
                    nargs='?',
                    dest='win_h_len',
                    default=250)
parser.add_argument(
    '--VCF_OUT',
    type=str,
    nargs='?',
    dest='vcf_out',
    # default='/Users/tschafers/Test_data/CNN/Results/G1_deepsv_indels.vcf')
    default=
    '/Users/lsantuari/Documents/Processed/Breakpoint2SV/Results/G1_DELs.vcf')

##################################
args = parser.parse_args()
# Window half length
win_hlen = args.win_h_len
# Window size
win_len = win_hlen * 2

strand = {False: '+', True: '-'}

bp_counter_sum = []


class Read:
    def __init__(self, location, strand):
        self.location = location
        self.strand = strand


class Location:
    def __init__(self, chrom, position):
        self.chrom = str(chrom)
        self.position = position

    # Return True if a read is clipped on the right or on the left
    def is_clipped_at(self, read):
        '''

        :param read: read object of the class pysam.AlignedSegment
        :return: True if the read is soft (4) or hard (5) clipped on the left or on the right, False otherwise
        '''
        if read.cigartuples is not None and read.reference_name == self.chrom:
            if is_left_clipped(read):
                return read.reference_start == self.position
            elif is_right_clipped(read):
                return read.reference_end == self.position

        return False


class Breakpoint:
    def __init__(self, chrom, position):

        self.location = Location(chrom, position)

        self.support = dict()
        for e in ['CR_R', 'CR_L', 'SR_L', 'SR_R', 'D', 'I', 'PE_L', 'PE_R']:
            self.support[e] = []
        self.mate_breakpoints = []
        self.strand = ''

    def get_distance_to_breakpoint(self, read):

        if read.is_reverse:
            return read.reference_start - self.location.position
        else:
            return self.location.position - read.reference_start

    def print(self):

        logging.info('Breakpoint => %s:%d' %
                     (self.location.chrom, self.location.position))

        for k in self.support.keys():
            if k in ['SR_L', 'SR_R', 'CR_L', 'CR_R']:
                clipped_read_count = Counter([
                    ':'.join(
                        [r.location.chrom,
                         str(r.location.position), r.strand])
                    for r in self.support[k]
                ])
                logging.info('Support %s => %s' %
                             (k, clipped_read_count.most_common(1)))
            else:
                logging.info('Support %s => %s' % (k, self.support[k]))

        if len(self.mate_breakpoints) > 0:
            mate_read_count = Counter([
                ':'.join(
                    [r.location.chrom,
                     str(r.location.position), r.strand])
                for r in self.mate_breakpoints
            ])
            logging.info('Support mate_breakpoints => %s' %
                         mate_read_count.most_common(1))

        logging.info('---------')


def read_bam():
    # Read BAM file
    logging.info('Reading bam file: ' + args.bam_file)
    # Check if the BAM file in input exists
    assert os.path.isfile(args.bam_file)
    aln = pysam.AlignmentFile(args.bam_file, "rb")
    return aln


def breakpoint_to_sv_v2(chr, breakpoints):
    # Read BAM
    aln = read_bam()

    ##Logging
    basename = os.path.splitext(os.path.basename(args.bed_file))[0]
    log_filename = args.out_dir + basename + '_' + str(chr) + '.log'

    # print('Writing %s' % log_filename)

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(filename=log_filename,
                        level=logging.INFO,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    logger = logging.getLogger()  # get the root logger
    logger.warning('This should go in the file.')

    logging.info('Processing Chr %s' % str(chr))
    logging.info('Processing %d candidate breakpoints' % len(breakpoints[chr]))

    breakpoint_list = []

    pos_counter = 1

    for pos in breakpoints[chr]:
        if not pos_counter % 1000:
            logging.info('Processed %d positions...' % pos_counter)

        current_brkpnt = Breakpoint(chr, pos)
        breakpoint_list.append(current_brkpnt)
        pos_counter = pos_counter + 1

        # print('Pos: %d' % pos)
        start = current_brkpnt.location.position - win_hlen
        end = current_brkpnt.location.position + win_hlen + 1
        count_reads = aln.count(current_brkpnt.location.chrom, start, end)
        # Fetch the reads mapped on the chromosome
        if count_reads <= args.max_read_count:
            logging.info(
                'Fetching %d reads in region %s:%d-%d' %
                (count_reads, current_brkpnt.location.chrom, start, end))
            for read in aln.fetch(current_brkpnt.location.chrom,
                                  start,
                                  end,
                                  multiple_iterators=True):
                # Both read and mate should be mapped
                if not read.is_unmapped and \
                        read.mapping_quality >= args.min_mapq:

                    if is_right_clipped(
                            read
                    ) and read.reference_end == current_brkpnt.location.position:
                        current_brkpnt.support['CR_R'].append(
                            Read(
                                Location(read.reference_name,
                                         read.reference_end),
                                strand[read.is_reverse]))
                    elif is_left_clipped(
                            read
                    ) and read.reference_start == current_brkpnt.location.position:
                        current_brkpnt.support['CR_L'].append(
                            Read(
                                Location(read.reference_name,
                                         read.reference_start),
                                strand[read.is_reverse]))

                    # INDEL case
                    if read.reference_start <= current_brkpnt.location.position <= read.reference_end:
                        if has_indels(read):

                            dels, ins = get_indels(read)

                            for i in ins:
                                if i[1] <= current_brkpnt.location.position <= i[
                                        2]:
                                    current_brkpnt.support['I'].append(i)
                                    current_brkpnt.mate_breakpoints.append(
                                        Breakpoint(read.reference_name, i[0]))

                            for d in dels:
                                if d[1] <= current_brkpnt.location.position <= d[
                                        2]:
                                    current_brkpnt.support['D'].append(d)
                                    current_brkpnt.mate_breakpoints.append(
                                        Breakpoint(read.reference_name, d[1]))

                    # PE case
                    if is_supporting_pe_read(read, current_brkpnt):

                        if read.is_reverse:
                            current_brkpnt.support['PE_L'].append(
                                current_brkpnt.get_distance_to_breakpoint(
                                    read))
                        else:
                            current_brkpnt.support['PE_R'].append(
                                current_brkpnt.get_distance_to_breakpoint(
                                    read))
                        current_brkpnt.mate_breakpoints.append(
                            Read(
                                Location(read.next_reference_name,
                                         read.next_reference_start),
                                strand_type(strand[read.is_reverse],
                                            strand[read.mate_is_reverse])))
                    # Supplementary Alignment (SA) case
                    if current_brkpnt.location.is_clipped_at(
                            read) and has_suppl_aln(read):
                        sa_chrom, sa_pos, sa_strand = get_suppl_aln(read)

                        if sa_chrom == current_brkpnt.location.chrom and sa_strand == strand[
                                read.is_reverse]:

                            if is_right_clipped(
                                    read) and sa_pos - read.reference_end > 0:
                                current_brkpnt.support['SR_R'].append(
                                    Read(
                                        Location(sa_chrom, sa_pos),
                                        strand_type(strand[read.is_reverse],
                                                    sa_strand)))
                            elif is_left_clipped(
                                    read) and sa_pos - read.reference_end < 0:
                                current_brkpnt.support['SR_L'].append(
                                    Read(
                                        Location(sa_chrom, sa_pos),
                                        strand_type(strand[read.is_reverse],
                                                    sa_strand)))
                            current_brkpnt.mate_breakpoints.append(
                                Read(Location(sa_chrom, sa_pos), sa_strand))

    for bp in breakpoint_list:
        bp.print()

    return breakpoint_list


##Accepts a list of arguments(breakpoints,chr)##
def breakpoint_to_sv(chr, breakpoints):
    aln = read_bam()
    ##Logging
    basename = os.path.splitext(os.path.basename(args.bed_file))[0]
    logging.basicConfig(filename=args.out_dir + basename + '_' + str(chr) +
                        '.log',
                        level=logging.INFO,
                        filemode='w',
                        format='%(asctime)s %(message)s',
                        datefmt='%m/%d/%Y %I:%M:%S %p')

    logging.info('Creating IntervalTree...')
    chr_tree = defaultdict(IntervalTree)
    bp_pos_dict = defaultdict(dict)
    # Create interval windows around candidate breakpoint positions
    for pos in breakpoints[chr]:
        chr_tree[chr][pos - win_hlen:pos + win_hlen + 1] = int(pos)
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
            logging.info('Fetching %d reads in region %s:%d-%d' %
                         (count_reads, chr, start, end))
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
                        match = chr_tree[read.next_reference_name][
                            read.next_reference_start]
                        if chr != read.next_reference_name or start > read.next_reference_start or \
                                end < read.next_reference_start:
                            if match:
                                for m in match:
                                    int_start, int_end, int_data = m
                                    # print('%s -> %s' % (pos, int_data))
                                    # Pay attention of double insertions. The same read pair will be added
                                    # from both intervals, leading to double count.
                                    links.append(
                                        frozenset({
                                            chr + '_' + str(pos) + '_' +
                                            strand[read.is_reverse],
                                            read.next_reference_name + '_' +
                                            str(int_data) + '_' +
                                            strand[read.mate_is_reverse]
                                        }))
                                    scanned_reads.add(read.query_name)

        # print('Right clipped:\n%s' % right_clipped_array)
        # print('Left clipped:\n%s' % left_clipped_array)
        if sum(right_clipped_array) + sum(left_clipped_array) == 0:
            # print('Pos %d has no CR pos' % pos)
            no_cr_pos.append(chr + '_' + str(pos))
        else:
            if max(right_clipped_array) > max(left_clipped_array):
                bp_pos_dict[chr][pos] = start + np.where(
                    right_clipped_array == max(right_clipped_array))[0][0] + 1
                # print('Right: %d -> %s' % (pos, max_i[0]))
            else:
                bp_pos_dict[chr][pos] = start + np.where(
                    left_clipped_array == max(left_clipped_array))[0][0] + 1
                # print('Left: %d -> %s' % (pos, max_i[0]))

        npos += 1

    links_counts = Counter(links)
    logging.info('Set size: %d' % len(links_counts))
    logging.info('No CR pos: %d' % len(no_cr_pos))
    logging.info('Connections with min 3 read pairs: %d' %
                 len([v for l, v in links_counts.items() if v > 2]))

    i = 0
    while len([v for l, v in links_counts.items() if v > i]) > 5000:
        i += 1
    logging.info('%d connections with min %d RP' %
                 (len([v for l, v in links_counts.items() if v > i]), i))
    # Return link positions, and counts
    logging.info('Finished breakpoint assembly for chr%s ' % (chr))
    ### Create result dict
    res_dict = defaultdict(dict)
    res_dict['links'] = links_counts
    res_dict['pos'] = bp_pos_dict
    return (res_dict)


def linksToVcf(links_counts, filename, ibam):
    cols = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n'
    # Write VCF Header
    with open(filename, 'w') as sv_calls:

        now = datetime.datetime.now()

        sv_calls.write('##fileformat=VCFv4.2\n')
        sv_calls.write('##FILTER=<ID=PASS,Description="All filters passed">\n')
        sv_calls.write('##fileDate=' + now.strftime("%Y-%m-%d %H:%M") + '\n')
        sv_calls.write('##fileTime=' + now.strftime("%H:%M") + '\n')
        sv_calls.write('##reference=GATK-GRCh-hg19')
        sv_calls.write('##ALT=<ID=DEL,Description="Deletion">\n')
        sv_calls.write(
            '##FORMAT=<ID=END,Number=1,Type=Integer,Description="End position of the structural variant">\n'
        )
        sv_calls.write(
            '##FORMAT=<ID=PE,Number=1,Type=Integer,Description="Paired-end support of the structural variant">\n'
        )
        sv_calls.write(
            '##FORMAT=<ID=SVTYPE,Number=1,Type=String,Description="Type of structural variant">\n'
        )

        # writing contig info
        chr_dict = get_chr_len_dict(ibam)
        for k in chr_dict.keys():
            sv_calls.write('##contig=<ID=' + k + ',length=' +
                           str(chr_dict[k]) + '>\n')

        sv_calls.write(cols)

        sv_evaluated = []

        for res in links_counts:
            links = res.get('links')
            position = res.get('pos')
            for l, v in Counter(links).items():
                logging.info('Candidate link => %s' % (l))
                interval = list(l)
                s1 = interval[0].split('_')
                s2 = interval[1].split('_')

                chrA = s1[0]
                chrB = s2[0]

                s1[1] = int(s1[1])
                s2[1] = int(s2[1])
                # Crashes if differnet from v >= 1
                if s1[1] in position[chrA] and s2[1] in position[
                        chrB] and v >= 1:
                    posA = position[chrA][s1[1]]
                    posB = position[chrB][s2[1]]

                    strandA = s1[2]
                    strandB = s2[2]

                    fs = frozenset({
                        str(chrA) + '_' + str(posA) + '_' + strandA,
                        str(chrB) + '_' + str(posB) + '_' + strandB
                    })

                    if fs not in sv_evaluated:
                        logging.info('Writing link => %s:%d-%s:%d' %
                                     (chrA, posA, chrB, posB))

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

                        # Get reference base
                        ref_base = get_ref_sequence(chrA, start)
                        sv_id = 'DEEPSV_' + chrA + '_' + str(start)

                        f_line = "SVTYPE:PE:END"
                        s_line = 'SVTYPE=%s;PE=%s;END=%s;STRAND=%s' % (
                            'DEL', v, stop, strandA + strandB)
                        line = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (
                            chrA, start, sv_id, ref_base, '<' + svtype + '>',
                            '.', 'PASS', s_line, '.', '.')
                        sv_calls.write(line + '\n')
                else:
                    logging.info('Link %s:%d-%s:%d considered already' %
                                 (chrA, posA, chrB, posB))
        print("VCF file written!")


def on_return(retval):
    bp_counter_sum.append(retval)


def main():

    breakpoints = read_breakpoints_single_position(args.bed_file)
    ##Spawn processes for each chromosome
    print('Found chromosomes:')
    for k in breakpoints.keys():
        print(k)
    print('Starting assembly')
    ###### Parallel execution ########
    start_time = time.time()

    # P = multiprocessing.Pool(processes=multiprocessing.cpu_count())
    # for chr in breakpoints.keys():
    #     P.apply_async(breakpoint_to_sv_v2, args=(chr,breakpoints), callback=on_return)
    # P.close()
    # P.join()

    #print("Writing intervals to VCF")
    #linksToVcf(bp_counter_sum, args.vcf_out, ibam = args.bam_file)

    for chr in breakpoints.keys():
        breakpoint_to_sv_v2(chr, breakpoints)

    print('Finished breakpoint assembly')
    # print("Sorting VCF file")
    # sort_command = "bcftools sort %s -o %s" % (args.vcf_out, args.vcf_out+"sorted.vcf")
    # call(sort_command)
    print("--- %s seconds ---" % (time.time() - start_time))

    # mychr = '17'
    # breakpoint_to_sv([mychr, breakpoints])


if __name__ == '__main__':
    main()
