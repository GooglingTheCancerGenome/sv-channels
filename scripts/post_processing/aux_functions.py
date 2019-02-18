'''
Generic functions used in the channel scripts
'''

import pysam
import os
from collections import defaultdict
import twobitreader as twobit

del_min_size = 30

def get_ref_sequence(chrname, pos):

    # Path on the HPC of the 2bit version of the human reference genome (hg19)
    #genome = twobit.TwoBitFile('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes/hg19.2bit')
    # Local path
    genome = twobit.TwoBitFile('/Users/lsantuari/Documents/Data/GiaB/reference/hg19.2bit')

    ref_pos = int(pos)-1
    if chrname[:2] == 'chr' or chrname[:2] == 'Chr':
        ref_base = genome['chr'+chrname[3:]][ref_pos].upper()
    else:
        ref_base = genome['chr' + chrname][ref_pos].upper()
    return ref_base


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


# Return start and end position of deletions and insertions
def get_indels(read):

    dels = []
    ins = []
    pos = read.reference_start
    if read.cigartuples is not None:
        for ct in read.cigartuples:
            # D is 2, I is 1
            if ct[0] == 2 and ct[1] >= del_min_size:
                dels.append(('D', pos, pos+ct[0]))
            if ct[0] == 1:
                ins.append(('I', pos, pos+ct[0]))
            pos = pos + ct[1]

    return dels, ins


def has_indels(read):

    if read.cigartuples is not None:
        cigar_set = set([ct[0] for ct in read.cigartuples])
        # D is 2, I is 1
        if len(set([1,2]) & cigar_set) > 0:
            return True
        else:
            return False
    else:
        return False


def is_supporting_pe_read(read, breakpoint):

    read_len = read.infer_query_length()
    if read_len is None:
        read_len = 0

    if not read.mate_is_unmapped:
        if read.next_reference_name == breakpoint.location.chrom:
            if (not read.is_reverse and
                    read.reference_end <= breakpoint.location.position <= read.next_reference_start ) or \
                (read.is_reverse and
                  read.next_reference_start + read_len <= breakpoint.location.position <= read.reference_start):
                return True
        else:
            if (not read.is_reverse and
                    read.reference_end <= breakpoint.location.position) or \
                (read.is_reverse and
                    breakpoint.location.position <= read.reference_start):
                return True
    return False


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


def has_suppl_aln(read):
    return read.has_tag('SA')


def strand_type(strand, mate_strand):
    if strand == mate_strand:
        return 'same'
    else:
        return 'opposite'

# Return chromosome and starting position of a supplementary alignment (split reads)
def get_suppl_aln(read):
    '''
    This function returns the chromosome and start position of the first supplementary alignment ('SA' tag) for a read.
    :param read: read object of the class pysam.AlignedSegment
    :return: a tuple with chromosome and start position of the first supplementary alignment. None if there are no
    supplementary alignments.
    '''
    if len(read.get_tag('SA')) > 0:
        # print(read.get_tag('SA'))
        supp_aln = read.get_tag('SA').split(';')[0]
        sa_info = supp_aln.split(',')
        # print(supp_aln)
        # print(sa_info)
        chr_sa = sa_info[0]
        start_sa = int(sa_info[1])
        strand_sa = sa_info[2]
        return chr_sa, start_sa, strand_sa
    else:
        return None


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
            if columns[3][:3] == "DEL":
                chrom = str(columns[0])
                breakpoints[chrom].append(int(columns[1]))
    for chr in breakpoints.keys():
        breakpoints[chr] = sorted(breakpoints[chr])
    # print(breakpoints)
    return breakpoints


def print_read(read, breakpoint):

    # DEBUG
    print('---\nSupporting read:')
    breakpoint.print()
    if read.is_reverse:
        print('<=')
    else:
        print('=>')
    print(read)