# Imports
import re

import pysam
from pysam import VariantFile

from collections import Counter
from intervaltree import IntervalTree
from collections import defaultdict
import numpy as np
import gzip
import bz2file
import os, errno
import pickle
import logging
import gzip

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = False

# Window half length
win_hlen = 100
# Window size
win_len = win_hlen * 2

__bpRE__ = None
__symbolicRE__ = None


class Breakpoint:

    def __init__(self, chr, pos, strand):
        self.chr = chr
        self.pos = pos
        self.strand = strand

    def id(self):
        return '_'.join([self.chr, str(self.pos), self.strand])


class StructuralVariant:

    def __init__(self, bp1, bp2):

        if bp1.chr == bp2.chr:
            if bp1.pos <= bp2.pos:
                self.tuple = (bp1, bp2)
            else:
                self.tuple = (bp2, bp1)
        elif bp1.chr < bp2.chr:
            self.tuple = (bp1, bp2)
        else:
            self.tuple = (bp2, bp1)

    def id(self):
        bp1, bp2 = self.tuple
        return '<=>'.join([bp1.id(), bp2.id()])


# Classes
class SVRecord_generic:

    def __init__(self, record, sv_caller):

        if type(record) != pysam.VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')
        # print(record)

        if record.info['SVTYPE'] == 'BND':
            ct, chr2, pos2, indellen = self.get_bnd_info(record)
        else:
            # print(record.info['SVTYPE'])
            ct = None
            chr2 = record.chrom
            pos2 = record.stop
            if 'SVLEN' in record.info.keys():
                indellen = record.info['SVLEN']
            else:
                indellen = abs(record.stop - record.pos)

        # print(record.info.keys())

        self.id = record.id
        self.chrom = record.chrom
        self.start = record.pos
        self.chrom2 = chr2
        self.end = pos2
        self.alt = record.alts[0]

        # CIPOS
        if 'CIPOS' in record.info.keys():
            if 'CIPOS95' in record.info.keys():
                self.cipos = record.info['CIPOS95']
            else:
                self.cipos = record.info['CIPOS']
        else:
            self.cipos = (0, 0)

        # CIEND
        if 'CIEND' in record.info.keys():
            if 'CIEND95' in record.info.keys():
                self.ciend = record.info['CIEND95']
            else:
                self.ciend = record.info['CIEND']
        elif 'CIRPOS' in record.info.keys():
            self.ciend = record.info['CIRPOS']
        else:
            self.ciend = (0, 0)

        self.filter = record.filter

        # Deletions are defined by 3to5 connection, same chromosome for start and end, start before end
        if self.chrom == self.chrom2:
            if ct == '3to5':
                if self.start < self.end:
                    self.svtype = 'DEL'
                elif self.start == self.end - 1:
                    self.svtype = 'INS'
                else:
                    self.svtype = record.info['SVTYPE']
            elif ct == '5to5' or ct == '3to3':
                self.svtype = 'INV'
            elif ct == '5to3':
                self.svtype = 'DUP'
            else:
                self.svtype = record.info['SVTYPE']
        else:
            self.svtype = 'BND'

    @staticmethod
    def stdchrom(chrom):

        if chrom[0] == 'c':
            return chrom[3:]
        else:
            return chrom

        # Modified from the function ctAndLocFromBkpt in mergevcf

    def locFromBkpt(self, ref, pre, delim1, pair, delim2, post):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''

        chpos = pair.split(':')
        # print(chpos[0])
        chr2 = self.stdchrom(chpos[0])
        pos2 = int(chpos[1])
        assert delim1 == delim2  # '['..'[' or ']'...']'
        joinedAfter = True
        extendRight = True
        connectSeq = ""

        if len(pre) > 0:
            connectSeq = pre
            joinedAfter = True
            assert len(post) == 0
        elif len(post) > 0:
            connectSeq = post
            joinedAfter = False

        if delim1 == "]":
            extendRight = False
        else:
            extendRight = True

        indellen = len(connectSeq) - len(ref)

        if joinedAfter:
            if extendRight:
                ct = '3to5'
            else:
                ct = '3to3'
        else:
            if extendRight:
                ct = '5to5'
            else:
                ct = '5to3'

        return ct, chr2, pos2, indellen

    def get_bnd_info(self, record):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''
        setupREs()

        altstr = str(record.alts[0])
        resultBP = re.match(__bpRE__, altstr)

        if resultBP:
            ct, chr2, pos2, indellen = self.locFromBkpt(str(record.ref), resultBP.group(1),
                                                        resultBP.group(2), resultBP.group(3), resultBP.group(4),
                                                        resultBP.group(5))
        return (ct, chr2, pos2, indellen)


class SVRecord_SUR:

    def __init__(self, record):
        if type(record) != pysam.VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')

        self.chrom = record.chrom
        self.chrom2 = record.info['CHR2']
        self.start = record.pos
        self.end = record.stop
        self.supp_vec = record.info['SUPP_VEC']
        self.svtype = record.info['SVTYPE']
        self.samples = record.samples


class SVRecord_nanosv:

    def __init__(self, record, sv_caller):

        if type(record) != pysam.VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')
        # print(record)

        ct, chr2, pos2, indellen = self.get_bnd_info(record)

        # print(record.info.keys())

        self.id = record.id
        self.chrom = record.chrom
        self.start = record.pos
        self.chrom2 = chr2
        self.end = pos2
        self.alt = record.alts[0]
        self.cipos = record.info['CIPOS']
        self.ciend = record.info['CIEND']
        self.filter = record.filter

        # Deletions are defined by 3to5 connection, same chromosome for start and end, start before end
        if ct == '3to5' and self.chrom == self.chrom2 and self.start <= self.end:
            self.svtype = 'DEL'
        else:
            self.svtype = record.info['SVTYPE']

    @staticmethod
    def stdchrom(chrom):

        if chrom[0] == 'c':
            return chrom[3:]
        else:
            return chrom

        # Modified from the function ctAndLocFromBkpt in mergevcf

    def locFromBkpt(self, ref, pre, delim1, pair, delim2, post):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''

        chpos = pair.split(':')
        # print(chpos[0])
        chr2 = self.stdchrom(chpos[0])
        pos2 = int(chpos[1])
        assert delim1 == delim2  # '['..'[' or ']'...']'
        joinedAfter = True
        extendRight = True
        connectSeq = ""

        if len(pre) > 0:
            connectSeq = pre
            joinedAfter = True
            assert len(post) == 0
        elif len(post) > 0:
            connectSeq = post
            joinedAfter = False

        if delim1 == "]":
            extendRight = False
        else:
            extendRight = True

        indellen = len(connectSeq) - len(ref)

        if joinedAfter:
            if extendRight:
                ct = '3to5'
            else:
                ct = '3to3'
        else:
            if extendRight:
                ct = '5to5'
            else:
                ct = '5to3'

        return ct, chr2, pos2, indellen

    def get_bnd_info(self, record):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''
        setupREs()

        altstr = str(record.alts[0])
        resultBP = re.match(__bpRE__, altstr)

        if resultBP:
            ct, chr2, pos2, indellen = self.locFromBkpt(str(record.ref), resultBP.group(1),
                                                        resultBP.group(2), resultBP.group(3), resultBP.group(4),
                                                        resultBP.group(5))
        return (ct, chr2, pos2, indellen)


def setupREs():
    '''
    Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
    URL: https://github.com/ljdursi/mergevcf
    '''
    global __symbolicRE__
    global __bpRE__
    if __symbolicRE__ is None or __bpRE__ is None:
        __symbolicRE__ = re.compile(r'.*<([A-Z:]+)>.*')
        __bpRE__ = re.compile(r'([ACGTNactgn\.]*)([\[\]])([a-zA-Z0-9\.]+:\d+)([\[\]])([ACGTNacgtn\.]*)')


def initialize_vcf_paths(sampleName):
    vcf_files = dict()

    if HPC_MODE:

        base_dir = '/hpc/cog_bioinf/kloosterman/shared/nanosv_comparison'

        if sampleName[:7] == 'NA12878':

            vcf_dir = os.path.join(base_dir, sampleName[:7])

            for mapper in ['bwa', 'minimap2', 'ngmlr', 'last']:

                vcf_files[mapper] = dict()
                vcf_files[mapper]['nanosv'] = os.path.join(vcf_dir, mapper, mapper + '_nanosv.sorted.vcf')
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

                vcf_files[mapper]['nanosv_sniffles_settings'] = os.path.join(vcf_dir, mapper,
                                                                             mapper +
                                                                             '_nanosv_with_sniffles_' +
                                                                             'settings.sorted.vcf')
                assert os.path.isfile(vcf_files[mapper]['nanosv_sniffles_settings'])

                if mapper in ['bwa', 'ngmlr']:
                    vcf_files[mapper]['sniffles'] = os.path.join(vcf_dir, mapper, mapper + '_sniffles.sorted.vcf')
                    assert os.path.isfile(vcf_files[mapper]['sniffles'])

                    vcf_files[mapper]['sniffles_nanosv_settings'] = os.path.join(vcf_dir, mapper,
                                                                                 mapper + '_sniffles_with_nanosv_' +
                                                                                 'settings.sorted.vcf')
                    assert os.path.isfile(vcf_files[mapper]['sniffles_nanosv_settings'])

    else:

        base_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth_no_INDELS/'

        if sampleName[:7] == 'NA12878':

            vcf_dir = os.path.join(base_dir, sampleName[:7], 'SV')

            for mapper in ['bwa', 'last']:

                vcf_files[mapper] = dict()
                vcf_files[mapper]['nanosv'] = os.path.join(vcf_dir, mapper, mapper + '_nanosv_pysam.sorted.vcf')
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

                if mapper in ['bwa']:
                    vcf_files[mapper]['sniffles'] = os.path.join(vcf_dir, mapper, mapper + '_sniffles.sorted.vcf')
                    assert os.path.isfile(vcf_files[mapper]['sniffles'])

        elif sampleName == 'Patient1' or sampleName == 'Patient2':

            vcf_dir = os.path.join(base_dir, sampleName, 'SV')

            for mapper in ['last']:
                vcf_files[mapper] = dict()
                vcf_files[mapper]['nanosv'] = os.path.join(vcf_dir, mapper, mapper + '_nanosv.sorted.vcf')
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

    return vcf_files


def read_nanosv_vcf(sampleName):
    '''
    This function parses the entries of the nanosv VCF file into SVRecord_nanosv objects and
    returns them in a list
    :param sampleName: str, name of the sample to consider
    :return: list, list of SVRecord_nanosv objects with the SV entries from the nanosv VCF file
    '''

    # Initialize regular expressions
    setupREs()

    # Setup locations of VCF files
    vcf_files = initialize_vcf_paths(sampleName)

    if sampleName[:7] == 'NA12878' or sampleName == 'Patient1' or sampleName == 'Patient2':

        # Reading the Last-mapped NanoSV VCF file
        filename = vcf_files['last']['nanosv']
        vcf_in = VariantFile(filename, 'r')

        sv = []

        # create sv list with SVRecord objects
        for rec in vcf_in.fetch():
            resultBP = re.match(__bpRE__, rec.alts[0])
            if resultBP:
                svrec = SVRecord_nanosv(rec, 'nanosv')
                sv.append(svrec)

        # Select good quality (no LowQual, only 'PASS') deletions (DEL)
        sv = [svrec for svrec in sv if svrec.svtype == 'DEL'
              # if 'LowQual' not in list(svrec.filter)]
              # if 'PASS' in list(svrec.filter)]
              if 'PASS' in list(svrec.filter) or \
              'CIPOS' in list(svrec.filter) or \
              'CIEND' in list(svrec.filter)]

        return sv
    else:
        return None


def load_bam(ibam):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    return pysam.AlignmentFile(ibam, "rb")


def get_chr_len_dict(ibam):
    bamfile = load_bam(ibam)
    # Extract chromosome length from the BAM header
    header_dict = bamfile.header

    chrLen = {i['SN']: i['LN'] for i in header_dict['SQ']}
    return chrLen


def read_bed_sv(inbed):
    # Check file existence
    assert os.path.isfile(inbed)
    # Dictionary with chromosome keys to store SVs
    sv_dict = defaultdict(list)

    with(open(inbed, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom = str(columns[0])
            if columns[3][:3] == "DEL":
                sv_dict[chrom].append((int(columns[1]), int(columns[2]), columns[3]))

    # print(sv_dict)
    return sv_dict


def read_bedpe_sv(inbed):
    #
    # Check file existence
    assert os.path.isfile(inbed)
    # Dictionary with chromosome keys to store SVs
    sv_dict = defaultdict(list)

    with(open(inbed, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            if columns[0] == columns[3]:
                chrom = str(columns[0])
                chrom = chrom.replace('chr', '')
                sv_dict[chrom].append(
                    (int(columns[1]), int(columns[2]), 'DEL_start',
                     int(columns[4]), int(columns[5]), 'DEL_end')
                )

    # print(sv_dict)
    return sv_dict


def load_candidate_pairs(sample, ibam):
    logging.info('Loading data for candidate_pairs:')

    chr_len = get_chr_len_dict(ibam)

    prefix = '' if HPC_MODE else '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth/'
    ch = 'candidate_pairs'

    candidate_pairs_file = os.path.join(prefix, sample, ch + '.pgz')

    if os.path.exists(candidate_pairs_file):

        # Save candidate_pairs
        with gzip.GzipFile(candidate_pairs_file, 'rb') as f:
            candidate_pairs = pickle.load(f)
        f.close()

    else:

        candidate_pairs = dict()

        for chrom in chr_len.keys():

            logging.info('Loading candidate_pairs for Chr%s' % chrom)
            suffix = '.npy.bz2' if ch == 'coverage' else '.pbz2'
            if HPC_MODE:
                filename = os.path.join(prefix, sample, ch, '_'.join([chrom, ch + suffix]))
            else:
                filename = os.path.join(prefix, sample, ch, chrom + '_' + ch + suffix)

            assert os.path.isfile(filename)

            # logging.info('Reading %s for Chr%s' % (ch, chrom))
            with bz2file.BZ2File(filename, 'rb') as f:
                if suffix == '.npy.bz2':
                    candidate_pairs[chrom] = np.load(f)
                else:
                    candidate_pairs[chrom] = pickle.load(f)
            # logging.info('End of reading')

        # Save candidate_pairs
        with gzip.GzipFile(candidate_pairs_file, 'wb') as f:
            pickle.dump(candidate_pairs, f)
        f.close()

    return candidate_pairs


def get_pairs_with_overlap(sv_list, pairs, mode):
    '''

    :param sv_list: list, list of SVs
    :param candidate_pairs: list, list of StructuralVariant
    :return: list, list of cStructuralVariant whose Breakpoint windows completely overlap the CIPOS interval
    and the CIEND interval
    '''

    def get_tree(pairs):
        # Tree with windows for candidate pairs
        print('Building tree with %d pairs' % len(pairs))

        tree = dict()
        tree['pair_bp1'] = IntervalTree()
        tree['pair_bp2'] = IntervalTree()
        # Populate tree
        for p in pairs:
            bp1, bp2 = p.tuple
            tree['pair_bp1'][bp1.pos - win_hlen:bp1.pos + win_hlen + 1] = p.id()
            tree['pair_bp2'][bp2.pos - win_hlen:bp2.pos + win_hlen + 1] = p.id()
        return tree

    def search_tree_with_sv(sv_list, tree):

        overlap = dict()

        if mode == 'VCF':

            overlap['cipos'] = [tree['pair_bp1'][var.start + var.cipos[0]: var.start + var.cipos[1] + 1]
                                for var in sv_list]
            overlap['ciend'] = [tree['pair_bp2'][var.end + var.ciend[0]: var.end + var.ciend[1] + 1]
                                for var in sv_list]
        elif mode == 'BED':

            overlap['cipos'] = [tree['pair_bp1'][bp1_start: bp1_end + 1]
                                for bp1_start, bp1_end, bp1_lab, bp2_start, bp2_end, bp2_lab in sv_list]
            overlap['ciend'] = [tree['pair_bp2'][bp2_start: bp2_end + 1]
                                for bp1_start, bp1_end, bp1_lab, bp2_start, bp2_end, bp2_lab in sv_list]

        return overlap

    def get_overlap(tree, sv_list):

        full_overlap_ids = []
        partial_overlap_ids = []

        rg_overlap = search_tree_with_sv(sv_list, tree)

        # print('Length rg_overlap for cipos: %d' % len(rg_overlap['cipos']))
        # print('Length rg_overlap for ciend: %d' % len(rg_overlap['ciend']))

        for rg_bp1, rg_bp2, sv in zip(rg_overlap['cipos'], rg_overlap['ciend'], sv_list):

            bp1_start, bp1_end, bp1_lab, bp2_start, bp2_end, bp2_lab = sv

            rg_bp1_id_set = set([i for s, e, i in rg_bp1])
            rg_bp2_id_set = set([i for s, e, i in rg_bp2])

            common_ids = rg_bp1_id_set & rg_bp2_id_set

            if len(common_ids) > 0:

                # are SV CIPOS start and end positions fully included in the pair start endpoints?
                set_bp1 = set()
                for r in rg_bp1:
                    s, e, i = r
                    if i in common_ids and s <= bp1_start and e >= bp1_end:
                        set_bp1.add(i)
                    else:
                        partial_overlap_ids.append(i)

                # are SV CIEND start and end positions fully included in the pair end endpoints?
                set_bp2 = set()
                for r in rg_bp2:
                    s, e, i = r
                    if i in common_ids and s <= bp2_start and e >= bp2_end:
                        set_bp2.add(i)
                    else:
                        partial_overlap_ids.append(i)

                full_overlap_ids.extend(list(set_bp1 & set_bp2))

        partial_overlap_ids = list(set(partial_overlap_ids))

        print('Total candidate pairs with partial CI overlap: %d' % len(partial_overlap_ids))
        print('Total candidate pairs with full CI overlap: %d' % len(full_overlap_ids))

        return partial_overlap_ids, full_overlap_ids

    t = get_tree(pairs)
    partial, full = get_overlap(t, sv_list)

    return partial, full


def create_dir(directory):
    '''
    Create a directory if it does not exist. Raises an exception if the directory exists.
    :param directory: directory to create
    :return: None
    '''
    try:
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
