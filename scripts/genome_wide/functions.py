import gzip
import json
import logging
import os
from collections import Counter
from itertools import groupby

import numpy as np
import pysam
import twobitreader as twobit

del_min_size = 50
ins_min_size = 50
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


def has_suppl_aln(read):
    return read.has_tag('SA')


# Return chromosome and starting position of a supplementary alignment (split reads)
def get_suppl_aln(read):
    '''
    This function returns the chromosome and start position of the first supplementary alignment ('SA' tag) for a read.
    :param read: read object of the class pysam.AlignedSegment
    :return: a tuple with chromosome and start position of the first supplementary alignment. None if there are no
    supplementary alignments.
    '''
    def query_len(cigar_string):
        """
        Given a CIGAR string, return the number of bases consumed from the
        query sequence.
        """
        read_consuming_ops = ("M", "D", "N", "=", "X")
        result = 0
        cig_iter = groupby(cigar_string, lambda chr: chr.isdigit())
        for _, length_digits in cig_iter:
            length = int(''.join(length_digits))
            op = next(next(cig_iter)[1])
            if op in read_consuming_ops:
                result += length
        return result

    if len(read.get_tag('SA')) > 0:

        # print(read.get_tag('SA'))
        # get first supplemental alignment
        supp_aln = read.get_tag('SA').split(';')[0]
        sa_info = supp_aln.split(',')
        # print(supp_aln)
        # print(sa_info)
        chr_sa = sa_info[0]
        start_sa = int(sa_info[1])
        strand_sa = sa_info[2]
        cigar_sa = sa_info[3]
        mapq_sa = sa_info[4]
        nm_sa = sa_info[5]

        if int(mapq_sa) < 10:
            return None

        # print('{} {} {}'.format(chr_sa, start_sa, strand_sa))
        start_sa -= 1

        return chr_sa, start_sa, strand_sa, cigar_sa
    else:
        return None


# Return start and end position of deletions and insertions
def get_indels(read):
    dels_start = []
    dels_end = []
    ins = []
    pos = read.reference_start
    if read.cigartuples is not None:
        for ct in read.cigartuples:
            # D is 2, I is 1
            if ct[0] == 2 and ct[1] >= del_min_size:
                # dels.append(('D', pos, pos+ct[0]))
                dels_start.append(pos)
                dels_end.append(pos + ct[1])
            if ct[0] == 1 and ct[1] >= ins_min_size:
                # ins.append(('I', pos, pos+ct[0]))
                ins.append(pos)
            pos = pos + ct[1]

    return dels_start, dels_end, ins


def has_indels(read):
    if read.cigartuples is not None:
        cigar_set = set([ct[0] for ct in read.cigartuples])
        # D is 2, I is 1
        if len(set([1, 2]) & cigar_set) > 0:
            return True
        else:
            return False
    else:
        return False


# Return the mate of a read. Get read mate from BAM file
def get_read_mate(read, bamfile):
    '''
    This function was used in a previous version of the code when we needed to check if the mate was clipped.
    Retrieving the mate for each read is time consuming and should be avoided when not strictly necessary.

    :param read: object of the class pysam.AlignedSegment whose mate needs to be retrieved
    :param bamfile: BAM file containing both the read and its mate
    :return: mate, object of the class pysam.AlignedSegment, if a mate for the read is found. Return None otherwise.
    '''
    # print(read)
    # The mate is located at:
    # chromosome: read.next_reference_name
    # positions: [read.next_reference_start, read.next_reference_start+1]
    # Fetch all the reads in that location and retrieve the mate
    iter = bamfile.fetch(read.next_reference_name,
                         read.next_reference_start,
                         read.next_reference_start + 1,
                         multiple_iterators=True)
    for mate in iter:
        # A read and its mate have the same query_name
        if mate.query_name == read.query_name:
            # Check if read is first in pair (read1) and mate is second in pair (read2) or viceversa
            if (read.is_read1 and mate.is_read2) or (read.is_read2
                                                     and mate.is_read1):
                # print('Mate is: ' + str(mate))
                return mate
    return None


def get_reference_sequence(HPC_MODE, REF_GENOME):
    if HPC_MODE:
        # Path on the HPC of the 2bit version of the human reference genome
        genome = twobit.TwoBitFile(
            os.path.join(
                '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes',
                REF_GENOME + '.2bit'))
    else:
        # Path on the local machine of the 2bit version of the human reference genome
        genome = twobit.TwoBitFile(
            os.path.join('/Users/lsantuari/Documents/Data/GiaB/reference',
                         REF_GENOME + '.2bit'))

    return genome


def is_flanked_by_n(chrname, pos, HPC_MODE, REF_GENOME):
    genome = get_reference_sequence(HPC_MODE, REF_GENOME)

    if "N" in genome['chr' + chrname][pos - 1:pos + 1].upper():
        return True
    else:
        return False


# Return a one-hot encoding for the chromosome region chr:start-stop
# with Ns encoded as 1 and other chromosomes encoded as 0
def get_one_hot_sequence(chrname, start, stop, nuc, HPC_MODE, REF_GENOME):
    genome = get_reference_sequence(HPC_MODE, REF_GENOME)

    # ltrdict = {'a': 1, 'c': 2, 'g': 3, 't': 4, 'n': 0}

    # N one-hot
    # ltrdict = {'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 1}
    # return np.array([ltrdict[x.lower()] for x in genome['chr'+chrname][start:stop]])

    if chrname == 'MT':
        chrname = 'M'

    chrname = chrname if REF_GENOME == 'GRCh38' else 'chr' + chrname

    return np.array([
        1 if x.lower() == nuc.lower() else 0
        for x in genome[chrname][start:stop]
    ],
                    dtype=np.uint8)


def get_one_hot_sequence_by_list(twobitfile, chrname, positions):
    genome = twobit.TwoBitFile(twobitfile)
    whole_chrom = str(genome[chrname])
    nuc_list = ['A', 'T', 'C', 'G', 'N']
    res = np.zeros(shape=(len(positions), len(nuc_list)), dtype=np.uint32)

    for i, nuc in enumerate(nuc_list, start=0):
        res[:, i] = np.array([
            1 if whole_chrom[pos].lower() == nuc.lower() else 0
            for pos in positions
        ])

    return res


# From https://github.com/joferkington/oost_paper_code/blob/master/utilities.py
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.
    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.
    Returns:
    --------
        mask : A numobservations-length boolean array.
    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:, None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh


def get_config_file():
    with open(os.path.join(os.path.dirname(__file__), 'parameters.json'),
              'r') as f:
        config = json.load(f)
    return config


def get_chr_list():

    chrlist = [str(c) for c in list(np.arange(1, 23))]
    chrlist.extend(['X', 'Y'])

    return chrlist


# def get_chr_list():
#
#     chrlist = ['chr' + str(c) for c in list(np.arange(1, 23))]
#     chrlist.extend(['chrX', 'chrY'])
#
#     return chrlist


def get_chr_len_dict(ibam):
    chr_list = get_chr_list()
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chr_dict = {i['SN']: i['LN'] for i in header_dict['SQ']}
    chr_dict = {k: v for k, v in chr_dict.items() if k in chr_list}

    return chr_dict


def load_clipped_read_positions_by_chr(sampleName, chrName, chr_dict, win_hlen,
                                       channel_dir):
    def get_filepath(vec_type):
        fn = os.path.join(channel_dir, sampleName, vec_type,
                          vec_type + '.json.gz')
        return fn

    logging.info('Loading SR positions for Chr%s' % chrName)

    with gzip.GzipFile(get_filepath('split_read_pos'), 'rb') as fin:
        positions, locations = json.loads(fin.read().decode('utf-8'))

    with gzip.GzipFile(get_filepath('clipped_read_pos'), 'rb') as fin:
        positions_cr = json.loads(fin.read().decode('utf-8'))

    # print(locations)
    locations = [(chr1, pos1, chr2, pos2)
                 for chr1, pos1, chr2, pos2 in locations
                 if chr1 in chr_dict.keys() and chr2 in chr_dict.keys()
                 and win_hlen <= pos1 <= (chr_dict[chr1] - win_hlen)
                 and win_hlen <= pos2 <= (chr_dict[chr2] - win_hlen)]

    positions_cr_l = set(
        [int(k) + 1 for k, v in positions_cr.items() if v >= min_CR_support])
    positions_cr_r = set(
        [int(k) - 1 for k, v in positions_cr.items() if v >= min_CR_support])
    positions_cr = positions_cr_l | positions_cr_r

    # for pos in positions_cr:
    #     print('{}:{}'.format(chrName, pos))

    # print(positions_cr)
    locations = [(chr1, pos1, chr2, pos2)
                 for chr1, pos1, chr2, pos2 in locations
                 if (chr1 == chrName and pos1 in positions_cr) or (
                     chr2 == chrName and pos2 in positions_cr)]

    logging.info('{} positions'.format(len(locations)))

    return locations


def load_all_clipped_read_positions_by_chr(sampleName, win_hlen, chr_dict,
                                           output_dir):

    cr_pos_file = os.path.join(
        output_dir, sampleName,
        'candidate_positions_' + sampleName + '.json.gz')

    if os.path.exists(cr_pos_file):

        logging.info('Loading existing candidate positions file...')

        with gzip.GzipFile(cr_pos_file, 'rb') as fin:
            cpos_list = json.loads(fin.read().decode('utf-8'))
        fin.close()

        return cpos_list

    else:

        cpos_list = []

        chrlist = get_chr_list()
        chr_list = chrlist if sampleName != 'T1' else ['17']

        for chrName in chr_list:
            logging.info(
                'Loading candidate positions for Chr{}'.format(chrName))
            cpos = load_clipped_read_positions(sampleName, chrName, chr_dict,
                                               win_hlen, output_dir)
            cpos_list.extend(cpos)
            logging.info('Candidate positions for Chr{}: {}'.format(
                chrName, len(cpos)))

        logging.info('Writing candidate positions file {}'.format(cr_pos_file))

        with gzip.GzipFile(cr_pos_file, 'wb') as f:
            f.write(json.dumps(cpos_list).encode('utf-8'))
        f.close()

        return cpos_list


def load_all_clipped_read_positions(win_hlen,
                                    svtype,
                                    chr_dict,
                                    output_dir,
                                    clipped_type="SR"):

    config = get_config_file()
    min_CR_support = config["DEFAULT"]["MIN_CR_SUPPORT"]

    # cr_pos_file = os.path.join(output_dir, svtype + '_candidate_positions_' + clipped_type + '.json.gz')
    #
    # if os.path.exists(cr_pos_file):
    #
    #     logging.info('Loading existing candidate positions file...')
    #
    #     with gzip.GzipFile(cr_pos_file, 'rb') as fin:
    #         cpos_list = json.loads(fin.read().decode('utf-8'))
    #     fin.close()
    #
    #     return cpos_list
    #
    # else:

    def get_filepath(vec_type):
        fn = os.path.join(output_dir, vec_type, vec_type + '.json.gz')
        return fn

    logging.info('Loading SR positions')

    chrlist = get_chr_list()
    chr_list = chrlist  # if sampleName != 'T1' else ['17']

    with gzip.GzipFile(get_filepath('split_reads'), 'rb') as fin:
        positions_with_min_support_ls, positions_with_min_support_rs, total_reads_coord_min_support_json, \
        split_reads, split_read_distance = json.loads(fin.read().decode('utf-8'))

    with gzip.GzipFile(get_filepath('clipped_read_pos'), 'rb') as fin:
        left_clipped_pos_cnt, right_clipped_pos_cnt = json.loads(
            fin.read().decode('utf-8'))

    if svtype == 'DEL':
        total_reads_coord_min_support = total_reads_coord_min_support_json['DEL'] + \
            total_reads_coord_min_support_json['INDEL_DEL']
    elif svtype == 'INS':
        total_reads_coord_min_support = total_reads_coord_min_support_json['INS'] + \
            total_reads_coord_min_support_json['INDEL_INS']
    else:
        total_reads_coord_min_support = total_reads_coord_min_support_json[
            svtype]

    locations_sr = dict()
    locations_cr_r = dict()
    locations_cr_l = dict()

    positions_cr = dict()

    for chrom in chr_list:

        if clipped_type == 'SR':

            locations_sr[chrom] = [
                (chr1, pos1, chr2, pos2)
                for chr1, pos1, chr2, pos2 in total_reads_coord_min_support
                if chr1 in chr_dict.keys() and chr2 in chr_dict.keys() and chr1
                == chrom and win_hlen <= pos1 <= (chr_dict[chr1] - win_hlen)
                and win_hlen <= pos2 <= (chr_dict[chr2] - win_hlen)
            ]

            if svtype in ['DEL', 'INV', 'DUP', 'TRA']:

                if chrom in left_clipped_pos_cnt.keys():
                    positions_cr_l = set([
                        int(k) for k, v in left_clipped_pos_cnt[chrom].items()
                        if v >= min_CR_support
                    ])
                else:
                    positions_cr_l = set()
                if chrom in right_clipped_pos_cnt.keys():
                    positions_cr_r = set([
                        int(k)
                        for k, v in right_clipped_pos_cnt[chrom].items()
                        if v >= min_CR_support
                    ])
                else:
                    positions_cr_r = set()

                positions_cr[chrom] = positions_cr_l | positions_cr_r

                # for pos in positions_cr:
                #     print('{}:{}'.format(chrName, pos))

                # print(positions_cr)
                locations_sr[chrom] = [
                    (chr1, pos1, chr2, pos2)
                    for chr1, pos1, chr2, pos2 in locations_sr[chrom]
                    if (chr1 == chrom and pos1 in positions_cr[chr1]) or (
                        chr2 == chrom and pos2 in positions_cr[chr2])
                ]

                logging.info('Chr{}: {} positions'.format(
                    chrom, len(locations_sr[chrom])))

        elif clipped_type == 'CR':

            if chrom in left_clipped_pos_cnt.keys():
                positions_cr_l = set([
                    int(k) for k, v in left_clipped_pos_cnt[chrom].items()
                    if v >= min_CR_support
                ])
            else:
                positions_cr_l = set()

            if chrom in right_clipped_pos_cnt.keys():
                positions_cr_r = set([
                    int(k) for k, v in right_clipped_pos_cnt[chrom].items()
                    if v >= min_CR_support
                ])
            else:
                positions_cr_r = set()

            if len(positions_cr_r) > 0:
                locations_cr_r[chrom] = [
                    (chrom, pos) for pos in sorted(list(positions_cr_r))
                ]
            if len(positions_cr_l) > 0:
                locations_cr_l[chrom] = [
                    (chrom, pos) for pos in sorted(list(positions_cr_l))
                ]

    if clipped_type == 'SR':

        cpos_list = []
        for chrom in chr_list:
            if chrom in locations_sr.keys():
                cpos_list.extend(locations_sr[chrom])

        logging.info('{} candidate positions'.format(len(cpos_list)))

        # logging.info('Writing candidate positions file {}'.format(cr_pos_file))

        # with gzip.GzipFile(cr_pos_file, 'wb') as f:
        #     f.write(json.dumps(cpos_list).encode('utf-8'))
        # f.close()

        return cpos_list

    elif clipped_type == 'CR':

        cpos_list_right = []
        cpos_list_left = []

        for chrom in chr_list:
            if chrom in locations_cr_r.keys():
                cpos_list_right.extend(locations_cr_r[chrom])

        for chrom in chr_list:
            if chrom in locations_cr_l.keys():
                cpos_list_left.extend(locations_cr_l[chrom])

        logging.info('Right-clipped: {} candidate positions'.format(
            len(cpos_list_right)))
        logging.info('Left-clipped: {} candidate positions'.format(
            len(cpos_list_left)))

        logging.info('Writing candidate positions file {}'.format(cr_pos_file))

        # with gzip.GzipFile(cr_pos_file, 'wb') as f:
        #     f.write(json.dumps((cpos_list_right, cpos_list_left)).encode('utf-8'))
        # f.close()

        return cpos_list_right, cpos_list_left
