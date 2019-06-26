# Imports
import argparse
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
from time import time
import pandas as pd
import json

# from plotnine import *
import pprint

import logging
import csv
import statistics

# Default BAM file for testing
# On the HPC
# wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
# inputBAM = wd + "T0_dedup.bam"
# Locally
# wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
# inputBAM = wd + "T1_dedup.bam"

ci_slop = 0

# Chromosome lengths for reference genome hg19/GRCh37
chrom_lengths = {'1': 249250621, '2': 243199373, '3': 198022430, '4': 191154276, '5': 180915260, '6': 171115067, \
                 '7': 159138663, '8': 146364022, '9': 141213431, '10': 135534747, '11': 135006516, '12': 133851895, \
                 '13': 115169878, '14': 107349540, '15': 102531392, '16': 90354753, '17': 81195210, '18': 78077248, \
                 '19': 59128983, '20': 63025520, '21': 48129895, '22': 51304566, 'X': 155270560, \
                 'Y': 59373566, 'MT': 16569}

__bpRE__ = None
__symbolicRE__ = None

with open('./genome_wide/parameters.json', 'r') as f:
    config = json.load(f)

HPC_MODE = config["DEFAULT"]["HPC_MODE"]
CANDIDATE_POSITIONS = config["DEFAULT"]["CANDIDATE_POSITIONS"]
# Window size
win_hlen = config["DEFAULT"]["WIN_HLEN"]
win_len = config["DEFAULT"]["WIN_HLEN"] * 2
# Only clipped read positions supported by at least min_cr_support clipped reads are considered
min_cr_support = config["DEFAULT"]["MIN_CR_SUPPORT"]
min_sr_support = config["DEFAULT"]["MIN_SR_SUPPORT"]


# Classes

class SVRecord_generic:

    def __init__(self, record):

        # For CHM[1|13] SVs
        svtype_dict = {'deletion': 'DEL', 'insertion': 'INS', 'inversion': 'INV'}

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
        self.chrom = record.chrom.replace('chr', '')
        self.start = record.pos
        self.chrom2 = chr2.replace('chr', '')
        self.end = pos2
        self.alt = record.alts[0]

        # CIPOS
        if 'CIPOS' in record.info.keys():
            if 'CIPOS95' in record.info.keys():
                self.cipos = record.info['CIPOS95']
            else:
                self.cipos = record.info['CIPOS']
        else:
            self.cipos = (-ci_slop, ci_slop)

        # CIEND
        if 'CIEND' in record.info.keys():
            if 'CIEND95' in record.info.keys():
                self.ciend = record.info['CIEND95']
            else:
                self.ciend = record.info['CIEND']
        elif 'CIRPOS' in record.info.keys():
            self.ciend = record.info['CIRPOS']
        else:
            self.ciend = (-ci_slop, ci_slop)

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

        if self.svtype in svtype_dict.keys():
            self.svtype = svtype_dict[self.svtype]

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


# class Label:
#
#     def __init__(self, chr, position, id, label_dict):
#         self.chr = chr
#         self.position = position
#         self.id = chr + '_' + position
#         self.label_dict = label_dict


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


def get_chr_len_by_chr(ibam, chrName):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    return chrLen


def get_chr_len_dict(ibam):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chr_dict = {i['SN']: i['LN'] for i in header_dict['SQ']}

    return chr_dict


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


def load_clipped_read_positions(sampleName, chrName):
    def get_filepath(vec_type):

        if HPC_MODE:
            channel_dir = ''
            fn = os.path.join(channel_dir, sampleName, vec_type, chrName + '_' + vec_type + '.pbz2')
        else:
            channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
            fn = '/'.join((channel_dir, sampleName, vec_type, chrName + '_' + vec_type + '.pbz2'))

        return fn

    # vec_type = 'clipped_read_pos' if CANDIDATE_POSITIONS == "CR" else 'split_read_pos'

    print('Loading CR positions for Chr%s' % chrName)
    # Load files

    with bz2file.BZ2File(get_filepath('clipped_read_pos'), 'rb') as f:
        cpos = pickle.load(f)

    print('Loading SR positions for Chr%s' % chrName)

    with bz2file.BZ2File(get_filepath('split_read_pos'), 'rb') as f:
        positions, locations = pickle.load(f)
        spos = positions

    # Filter by minimum support
    # cr_pos = [elem for elem, cnt in cpos.items() if cnt >= min_cr_support]

    print('Calculating CR positions for min support {}, length={}'.format(min_cr_support, len(cpos)))
    cr_pos = [k for k, v in cpos.items() if v >= min_cr_support]
    print('Calculating SR positions for min support {}, length={}'.format(min_sr_support, len(spos)))
    sr_pos = [k for k, v in spos.items() if v >= min_sr_support]
    print('Calculating CR positions in SR positions')
    cr_pos = sorted(list(set(cr_pos) & set(sr_pos)))
    print('Final CR positions={}'.format(len(cr_pos)))

    # Remove positions with windows falling off chromosome boundaries
    # print(f'win_hlen = {win_hlen}, chrom_lengths[{chrName}] = {chrom_lengths[chrName]}')
    cr_pos = [pos for pos in cr_pos if win_hlen <= pos <= (chrom_lengths[chrName] - win_hlen)]

    return cr_pos


def load_all_clipped_read_positions(sampleName):
    cr_pos_dict = {}
    for chrName in chrom_lengths.keys():
        # for chrName in ['22']:
        cr_pos_dict[chrName] = load_clipped_read_positions(sampleName, chrName)
    return cr_pos_dict


def initialize_nanosv_vcf_paths(sampleName):
    vcf_files = dict()

    if HPC_MODE:

        if sampleName[:7] == 'NA12878':

            # vcf_dir = '/hpc/cog_bioinf/kloosterman/shared/nanosv_comparison/' + sampleName[:7]
            vcf_dir_last = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Data_for_labels/',
                                        sampleName, 'VCF')

            # for mapper in ['bwa', 'minimap2', 'ngmlr', 'last']:
            for mapper in ['last']:

                vcf_files[mapper] = dict()

                if mapper == 'last':
                    vcf_files[mapper]['nanosv'] = os.path.join(vcf_dir_last, mapper + '_nanosv.sorted.vcf')
                # else:
                #     vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv.sorted.vcf'
                # assert os.path.isfile(vcf_files[mapper]['nanosv'])
                #
                # vcf_files[mapper]['nanosv_sniffles_settings'] = vcf_dir + '/' + mapper + '/' + \
                #                                                 mapper + '_nanosv_with_sniffles_settings.sorted.vcf'
                # assert os.path.isfile(vcf_files[mapper]['nanosv_sniffles_settings'])
                #
                # if mapper in ['bwa', 'ngmlr']:
                #     vcf_files[mapper]['sniffles'] = vcf_dir + '/' + mapper + '/' + mapper + '_sniffles.sorted.vcf'
                #     assert os.path.isfile(vcf_files[mapper]['sniffles'])
                #
                #     vcf_files[mapper]['sniffles_nanosv_settings'] = vcf_dir + '/' + mapper + '/' + \
                #                                                     mapper + '_sniffles_with_nanosv_settings.sorted.vcf'
                #     assert os.path.isfile(vcf_files[mapper]['sniffles_nanosv_settings'])

    else:

        if sampleName[:7] == 'NA12878':

            vcf_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth_no_INDELS/' + sampleName[:7] + '/SV'
            vcf_files = dict()

            for mapper in ['bwa', 'last']:

                vcf_files[mapper] = dict()

                vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv_pysam.sorted.vcf'
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

                if mapper in ['bwa']:
                    vcf_files[mapper]['sniffles'] = vcf_dir + '/' + mapper + '/' + mapper + '_sniffles.sorted.vcf'
                    assert os.path.isfile(vcf_files[mapper]['sniffles'])

        elif sampleName == 'Patient1' or sampleName == 'Patient2':

            vcf_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth_no_INDELS/' + sampleName + '/SV'
            vcf_files = dict()

            for mapper in ['last']:
                vcf_files[mapper] = dict()

                vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv.sorted.vcf'
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
    vcf_files = initialize_nanosv_vcf_paths(sampleName)

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

        # othersv = [svrec for svrec in sv if svrec.svtype != 'DEL']

        # How many distinct FILTERs?
        # filter_set = set([f for svrec in sv for f in svrec.filter])
        # print(filter_set)

        # How many VCF record with a specific FILTER?
        # filter_list = sorted(filter_set)
        # s = pd.Series([sum(list(map(lambda x: int(f in x),
        #                             [list(svrec.filter) for svrec in sv])))
        #                for f in filter_list],
        #               index=filter_list)
        # s = s.append(pd.Series([len(sv)], index=['Total']))
        # s = s.sort_values()
        # Print pd.Series with stats on FILTERs
        # print(s)

        return sv


def read_vcf(sampleName, sv_caller):
    '''
    This function parses the entries of the nanosv VCF file into SVRecord_nanosv objects and
    returns them in a list
    :param sampleName: str, name of the sample to consider
    :return: list, list of SVRecord_nanosv objects with the SV entries from the nanosv VCF file
    '''

    # Initialize regular expressions
    setupREs()

    # Setup locations of VCF files

    if HPC_MODE:

        if sampleName == 'NA24385':
            filename = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/',
                                    'Datasets/GiaB/HG002_NA24385_son/NIST_SVs_Integration_v0.6/',
                                    'processed/HG002_SVs_Tier1_v0.6.PASS.vcf')

        elif sampleName in ['CHM1', 'CHM13']:
            filename = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/Huddleston2016/',
                                    'structural_variants/',sampleName+'_SVs.annotated.vcf.gz')
        else:
            filename = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Data_for_labels',
                                    sampleName, 'VCF', sv_caller + '.sym.vcf')
    else:

        if sampleName[:7] == 'NA12878':
            filename = '/Users/lsantuari/Documents/Data/germline/trio/' + \
                       sampleName[:7] + '/SV/Filtered/' + sv_caller + '.sym.vcf'
        elif sampleName == 'Patient1' or sampleName == 'Patient2':
            filename = '/Users/lsantuari/Documents/Data/germline/patients/' + \
                       sampleName + '/SV/Filtered/' + sv_caller + '.sym.vcf'

    print('Reading VCF file %s\nfor SV caller %s' % (filename, sv_caller))
    vcf_in = VariantFile(filename, 'r')

    sv = []

    # create sv list with SVRecord objects
    for rec in vcf_in.fetch():
        svrec = SVRecord_generic(rec)
        sv.append(svrec)

    print('SVs read: %d' % len(sv))
    print(Counter([svrec.svtype for svrec in sv]))

    # Select good quality (no LowQual, only 'PASS') deletions (DEL)
    sv = [svrec for svrec in sv if svrec.svtype == 'DEL'
          # if 'LowQual' not in list(svrec.filter)]
          # if 'PASS' in list(svrec.filter)
          ]

    return sv


def get_labels_from_nanosv_vcf(sampleName):
    '''
    This function writes the label files based on the nanosv VCF file information

    :param sampleName: str, name of the sample considered
    :param chrName: str, chromosome name
    :param ibam: str, path of the BAM file in input
    :return: None
    '''

    # Lines to write in the BED file
    lines = []

    # def closest_loc(pos, pos_list):
    #     '''
    #
    #     :param pos: reference position
    #     :param pos_list: list of clipped read positions
    #     :return: tuple of (position, distance) for the closest position to the reference position
    #     '''
    #     pos_array = np.asarray(pos_list)
    #     deltas = np.abs(pos_array - pos)
    #     idx = np.argmin(deltas)
    #     return (pos_list[idx], deltas[idx])

    # Load SV list
    sv_list = read_nanosv_vcf(sampleName)

    # Select deletions (DELs)
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL']

    # list of chromosomes
    chr_list = set([var.chrom for var in sv_list])

    # print('Plotting CI distribution')
    plot_ci_dist(sv_list, sampleName)

    # print(chr_list)
    print('Total # of DELs: %d' % len(sv_list))

    cnt = Counter([sv.chrom for sv in sv_list])
    chr_series = pd.Series([v for v in cnt.values()], index=cnt.keys())
    # print('# SVs per chromosome:')
    # print(chr_series)

    assert sum(chr_series) == len(sv_list)

    confint = 100

    labels_list = defaultdict(list)

    for chrName in chr_list:

        sv_list_chr = [var for var in sv_list if var.chrom == chrName]

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # print(sorted(cr_pos))

        # Using IntervalTree for interval search
        t = IntervalTree()

        # print('# SVs in Chr: %d' % len(sv_list_chr))

        for var in sv_list_chr:
            # cipos[0] and ciend[0] are negative in the VCF file
            id_start = var.svtype + '_start'
            id_end = var.svtype + '_end'

            # id_start = '_'.join((var.chrom, str(var.start+var.cipos[0]),  str(var.start+var.cipos[1])))
            # id_end = '_'.join((var.chrom, str(var.end + var.ciend[0]), str(var.end+var.ciend[1])))
            assert var.start <= var.end, "Start: " + str(var.start) + " End: " + str(var.end)

            # print('var start -> %s:%d CIPOS: (%d, %d)' % (chrName, var.start, var.cipos[0], var.cipos[1]))
            # print('var end -> %s:%d CIEND: (%d, %d)' % (chrName, var.end, var.ciend[0], var.ciend[1]))

            t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
            t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

            # t[var.start - confint:var.start + confint + 1] = var.svtype + '_start'
            # t[var.end - confint:var.end + confint + 1] = var.svtype + '_end'

        label_search = [sorted(t[p - win_hlen: p + win_hlen + 1]) for p in cr_pos]

        crpos_full_ci, crpos_partial_ci = get_crpos_win_with_ci_overlap(sv_list_chr, cr_pos)

        # print('Clipped read positions with complete CI overlap: %s' % crpos_full_ci)
        # print('Clipped read positions with partial CI overlap: %s' % crpos_partial_ci)
        # crpos_ci_isec = set(crpos_full_ci) & set(crpos_partial_ci)
        # print('Intersection: %s' % crpos_ci_isec)

        print('# CRPOS in CI: %d' % len([l for l in label_search if len(l) != 0]))

        count_zero_hits = 0
        count_multiple_hits = 0

        label_ci_full_overlap = []

        for elem, pos in zip(label_search, cr_pos):
            if len(elem) == 1:
                # print(elem)
                if pos in crpos_full_ci:
                    label_ci_full_overlap.append(elem[0].data)

                    lines.append(bytes(chrName + '\t' + str(elem[0].begin) + '\t' \
                                       + str(elem[0].end) + '\t' + elem[0].data + '\n', 'utf-8'))

                elif pos in crpos_partial_ci:
                    label_ci_full_overlap.append('UK')
                else:
                    label_ci_full_overlap.append('noSV')
            elif len(elem) == 0:
                count_zero_hits += 1
                label_ci_full_overlap.append('noSV')
            elif len(elem) > 1:
                count_multiple_hits += 1
                label_ci_full_overlap.append('UK')
                # if pos in crpos_full_ci:
                #     label_ci_full_overlap.append('Multiple_Full')
                #     #print('Multiple full: %s -> %s' % ( [d for s,e,d in elem], set([d for s,e,d in elem]) ) )
                #     #for s, e, d in elem:
                #     #    print('%d -> %d %d %s' % (pos, s, e, d))
                # #else:
                #     #label_ci_full_overlap.append('Multiple_Partial')

        print('CR positions: %d' % len(cr_pos))
        print('Label length: %d' % len(label_search))
        assert len(label_ci_full_overlap) == len(cr_pos)

        print('Label_CI_full_overlap: %s' % Counter(label_ci_full_overlap))
        print('Zero hits:%d' % count_zero_hits)
        print('Multiple hits:%d' % count_multiple_hits)

        # Write labels for chromosomes
        if not HPC_MODE:
            channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
        else:
            channel_dir = '.'

        output_dir = '/'.join((channel_dir, sampleName, 'label'))
        create_dir(output_dir)

        # print(output_dir)

        with gzip.GzipFile('/'.join((output_dir, chrName + '_label_ci_full_overlap.npy.gz')), "w") as f:
            np.save(file=f, arr=label_ci_full_overlap)
        f.close()

    # Write BED file with labelled CI positions
    outfile = sampleName + '_nanosv_vcf_ci_labelled.bed.gz'
    f = gzip.open(outfile, 'wb')
    try:
        # use set to make lines unique
        for l in set(lines):
            f.write(l)
    finally:
        f.close()


def write_sv_without_cr(sampleName, ibam):
    '''
    Writes NanoSV SVs with no clipped read support into a BED file
    :param sampleName: name of the sample to consider
    :param ibam: path of the BAM file. Used to get chromosome lengths
    :return: None
    '''

    # Load SV list
    sv_list = read_nanosv_vcf(sampleName)
    # Select deletions
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL' if sv.chrom == sv.chrom2 if sv.start < sv.end]
    # list of chromosomes
    chr_list = set([var.chrom for var in sv_list])

    var_with_cr = 0

    bedout = open(sampleName + '_nanosv_no_cr.bed', 'w')

    for chrName in sorted(chr_list):

        sv_list_chr = [var for var in sv_list if var.chrom == chrName]

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # Using IntervalTree for interval search
        t = IntervalTree()

        for var in sv_list_chr:
            id_start = var.svtype + '_start'
            id_end = var.svtype + '_end'

            t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
            t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

        hits = [sorted(t[p]) for p in cr_pos]

        inter_list = [nested_elem for elem in hits for nested_elem in elem]

        start_var_list = []
        end_var_list = []

        for start, end, data in inter_list:
            se = data.split('_')[1]
            if se == 'start':
                start_var_list.append(start)
            elif se == 'end':
                end_var_list.append(start)

        for var in sv_list_chr:

            if var.start + var.cipos[0] in start_var_list and \
                    var.end + var.ciend[0] in end_var_list:
                var_with_cr += 1

            if var.start + var.cipos[0] not in start_var_list:
                bedout.write('\t'.join((var.chrom,
                                        str(var.start + var.cipos[0]),
                                        str(var.start + var.cipos[1]))) + '\n')
            if var.end + var.ciend[0] not in end_var_list:
                bedout.write('\t'.join((var.chrom,
                                        str(var.end + var.ciend[0]),
                                        str(var.end + var.ciend[1]))) + '\n')

    bedout.close()

    print('VCF entries with CR on both sides: %d/%d' % (var_with_cr, len(sv_list)))


def plot_ci_dist(sv_list, sampleName):
    '''
    Saves the plots of the distributions of the confidence intervals reported by NanoSV
    :param sv_list: list, a list of SVs
    :param sampleName: str, name of the sample to consider
    :return: None
    '''

    # Plot distribution of CIPOS and CIEND
    df = pd.DataFrame({
        "cipos": np.array([var.cipos[1] + abs(var.cipos[0]) for var in sv_list]),
        "ciend": np.array([var.ciend[1] + abs(var.ciend[0]) for var in sv_list])
    })

    print('Max CIPOS:%d, max CIEND:%d' % (max(df['cipos']), max(df['ciend'])))

    output_dir = '/Users/lsantuari/Documents/Data/germline/plots'
    # the histogram of the data
    p = ggplot(aes(x='cipos'), data=df) + \
        geom_histogram(binwidth=1) + ggtitle(' '.join((sampleName, 'CIPOS', 'distribution')))
    p.save(filename='_'.join((sampleName, 'CIPOS', 'distribution')), path=output_dir)

    p = ggplot(aes(x='ciend'), data=df) + \
        geom_histogram(binwidth=1) + ggtitle(' '.join((sampleName, 'CIEND', 'distribution')))
    p.save(filename='_'.join((sampleName, 'CIEND', 'distribution')), path=output_dir)


def get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName):
    '''

    :param sampleName: sample ID
    :return: list of SVRecord_nanosv that overlap Manta SVs as reported in SURVIVOR merge VCF file
    '''
    sv_sur = read_SURVIVOR_merge_VCF(sampleName)
    sv_nanosv = read_nanosv_vcf(sampleName)

    common_sv = defaultdict(list)
    for sv in sv_sur:
        # print(sv.supp_vec)
        # 0:Delly, 1:GRIDSS, 2:last_nanosv, 3:Lumpy, 4:Manta
        if sv.supp_vec[2] == '1' and sv.supp_vec[4] == '1':
            # print(sv.samples.keys())
            svtype = sv.samples.get('NanoSV').get('TY')
            coord = sv.samples.get('NanoSV').get('CO')
            if svtype == 'DEL' and coord != 'NaN':
                coord_list = re.split('-|_', coord)
                assert len(coord_list) == 4
                common_sv[coord_list[0]].append(int(coord_list[1]))
    # print(common_sv.keys())
    sv_nanosv_manta = [sv for sv in sv_nanosv
                       if sv.chrom in common_sv.keys() if sv.start in common_sv[sv.chrom]]
    # print(common_sv['1'])
    # print(len(sv_nanosv_manta))
    return sv_nanosv_manta


# END: NanoSV specific functions

# START: BED specific functions

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


def get_labels_from_bed(sampleName, ibam, inbed):
    '''

    :param sampleName: name of sample considered
    :param ibam: path to the BAM file. Needed to get chromosome length
    :param inbed: path to the BED file with SVs
    :return: dictionary with list of labels per chromosome
    '''

    print('sample = %s' % sampleName)
    print('window = %d' % win_len)

    sv_list = read_bed_sv(inbed)

    # chr_list = sv_list.keys()
    # Use Chr1 for testing
    chr_list = ['1']

    labels_list = defaultdict(list)

    for chrName in chr_list:

        sv_list_chr = sv_list[chrName]

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # print(sorted(cr_pos))

        # Using IntervalTree for interval search
        t = IntervalTree()

        # print('# Breakpoints in Chr: %d' % len(sv_list_chr))

        for start, end, lab in sv_list_chr:
            t[start:end + 1] = lab

        label = [sorted(t[p - win_hlen: p + win_hlen + 1]) for p in cr_pos]

        crpos_full_ci, crpos_partial_ci = get_crpos_win_with_bed_overlap(sv_list_chr, cr_pos)

        # print('Clipped read positions with complete CI overlap: %s' % crpos_full_ci)
        # print('Clipped read positions with partial CI overlap: %s' % crpos_partial_ci)

        crpos_ci_isec = set(crpos_full_ci) & set(crpos_partial_ci)
        # print('Intersection should be empty: %s' % crpos_ci_isec)
        assert len(crpos_ci_isec) == 0

        # print('# CRPOS in BED: %d' % len([l for l in label if len(l) != 0]))

        count_zero_hits = 0
        count_multiple_hits = 0

        label_ci_full_overlap = []

        for elem, pos in zip(label, cr_pos):
            # Single match
            if len(elem) == 1:
                # print(elem)
                if pos in crpos_full_ci:
                    label_ci_full_overlap.append(elem[0].data)
                elif pos in crpos_partial_ci:
                    label_ci_full_overlap.append('UK')
                else:
                    label_ci_full_overlap.append('noSV')
            # No match
            elif len(elem) == 0:
                count_zero_hits += 1
                label_ci_full_overlap.append('noSV')
            # Multiple match
            elif len(elem) > 1:
                count_multiple_hits += 1
                label_ci_full_overlap.append('UK')

        # print('CR positions: %d' % len(cr_pos))
        # print('Label length: %d' % len(label))
        assert len(label_ci_full_overlap) == len(cr_pos)

        # print('Label_CI_full_overlap: %s' % Counter(label_ci_full_overlap))
        # print('Zero hits:%d' % count_zero_hits)
        # print('Multiple hits:%d' % count_multiple_hits)

        # if not HPC_MODE:
        #     channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
        # else:
        #     channel_dir = ''
        #
        # output_dir = '/'.join((channel_dir, sampleName, 'label'))
        # create_dir(output_dir)

        # print(output_dir)

        # with gzip.GzipFile('/'.join((output_dir, chrName + '_label_ci_full_overlap.npy.gz')), "w") as f:
        #     np.save(file=f, arr=label_ci_full_overlap)
        # f.close()

        labels_list[chrName] = label_ci_full_overlap

    return labels_list


def get_crpos_win_with_ci_overlap(sv_list, cr_pos):
    '''

    :param sv_list: list, list of SVs
    :param cr_pos: list, list of clipped read positions
    :return: list, list of clipped read positions whose window completely overlap either the CIPOS interval
    or the CIEND interval
    '''

    def get_tree(cr_pos):
        # Tree with windows for CR positions
        tree = IntervalTree()
        # Populate tree
        for pos in cr_pos:
            tree[pos - win_hlen:pos + win_hlen + 1] = pos
        return tree

    def search_tree_with_sv(sv_list, tree, citype):

        if citype == 'CIPOS':
            return [sorted(tree[var.start + var.cipos[0]: var.start + var.cipos[1] + 1])
                    for var in sv_list]
        elif citype == 'CIEND':
            return [sorted(tree[var.end + var.ciend[0]: var.end + var.ciend[1] + 1])
                    for var in sv_list]

    def get_overlap(tree, sv_list, citype):

        rg_overlap = search_tree_with_sv(sv_list, tree, citype)

        if citype == 'CIPOS':
            start_ci = [var.start + var.cipos[0] for var in sv_list]
            end_ci = [var.start + var.cipos[1] for var in sv_list]
        elif citype == 'CIEND':
            start_ci = [var.end + var.ciend[0] for var in sv_list]
            end_ci = [var.end + var.ciend[1] for var in sv_list]

        full = []
        partial = []
        for rg, start, end in zip(rg_overlap, start_ci, end_ci):
            for elem in rg:
                elem_start, elem_end, elem_data = elem
                if start >= elem_start and end <= elem_end:
                    # print('CIPOS->Full: %s\t%d\t%d' % (elem, start, end))
                    full.append(elem_data)
                else:
                    # print('CIPOS->Partial: %s\t%d\t%d' % (elem, start, end))
                    partial.append(elem_data)

        return partial, full

    t = get_tree(cr_pos)
    partial_cipos, full_cipos = get_overlap(t, sv_list, 'CIPOS')
    partial_ciend, full_ciend = get_overlap(t, sv_list, 'CIEND')

    cr_full_overlap = sorted(full_cipos + full_ciend)
    cr_partial_overlap = sorted(partial_cipos + partial_ciend)

    return sorted(list(set(cr_full_overlap))), sorted(list(set(cr_partial_overlap) - set(cr_full_overlap)))


def get_crpos_win_with_bed_overlap(sv_list, cr_pos):
    '''
    :param sv_list: list, list of SV bed intervals
    :param cr_pos: list, list of clipped read positions
    :return: list, list of clipped read positions whose window completely overlap either the CIPOS interval
    or the CIEND interval
    '''
    # Tree with windows for CR positions
    t_cr = IntervalTree()

    for pos in cr_pos:
        t_cr[pos - win_hlen:pos + win_hlen + 1] = pos

    cr_full_overlap = []
    cr_partial_overlap = []

    rg_overlap = [sorted(t_cr[start: end + 1]) for start, end, lab in sv_list]
    # print('Range overlap: %s' % rg_overlap)

    for rg, start, end in zip(rg_overlap,
                              [start for start, end, lab in sv_list],
                              [end for start, end, lab in sv_list]):
        for elem in rg:
            elem_start, elem_end, elem_data = elem
            if start >= elem_start and end <= elem_end:
                cr_full_overlap.append(elem_data)
            else:
                cr_partial_overlap.append(elem_data)

    cr_full_overlap = sorted(cr_full_overlap)
    cr_partial_overlap = sorted(cr_partial_overlap)

    return sorted(list(set(cr_full_overlap))), sorted(list(set(cr_partial_overlap) - set(cr_full_overlap)))


# END: BED specific functions

def read_SURVIVOR_merge_VCF(sampleName):
    '''
    Reads the SURVIVOR merge VCF output and returns a list of SVs in as SVRecord_SUR objects
    :param sampleName: sample to consider
    :return: a list of SVs in as SVRecord_SUR objects
    '''

    if HPC_MODE:
        # To fill with HPC path
        survivor_vcf = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari',
                                    'Processed/Data_for_labels', sampleName, 'SURVIVOR',
                                    'survivor_merge.vcf')
    else:
        if sampleName[:7] == 'NA12878':
            context = 'trio'
            survivor_vcf = '/Users/lsantuari/Documents/Data/germline/' + context + \
                           '/' + sampleName[:7] + '/SV/Filtered/survivor_merge.vcf'
        elif sampleName == 'Patient1' or sampleName == 'Patient2':
            context = 'patients'
            survivor_vcf = '/Users/lsantuari/Documents/Data/germline/' + context + \
                           '/' + sampleName + '/SV/Filtered/survivor_merge.vcf'

    vcf_in = VariantFile(survivor_vcf)
    samples_list = list((vcf_in.header.samples))
    samples = samples_list
    # samples = [w.split('_')[0].split('/')[1] for w in samples_list]
    # print(samples)

    sv = []

    # create sv list with SVRecord_SUR objects
    for rec in vcf_in.fetch():
        # avoid SVs on chromosomes Y and MT
        if rec.chrom not in ['Y', 'MT'] and rec.info['CHR2'] not in ['Y', 'MT']:
            # print(rec)
            # print(dir(rec))
            sv.append(SVRecord_SUR(rec))

    return sv


def load_NoCR_positions():
    '''
    This function provides an overview of SV positions without clipped read support that are stored in the
    no_clipped_read_positions file.
    :return: None
    '''

    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "r") as f:
        no_clipped_read_pos = pickle.load(f)
    f.close()
    print(list(no_clipped_read_pos))


# Methods to save to BED format

def clipped_read_positions_to_bed(sampleName, ibam):
    chrlist = list(map(str, range(1, 23)))
    chrlist.extend(['X', 'Y'])
    # print(chrlist)

    lines = []
    for chrName in chrlist:
        crpos_list = load_clipped_read_positions(sampleName, chrName)
        lines.extend(
            [bytes(chrName + '\t' + str(crpos) + '\t' + str(crpos + 1) + '\n', 'utf-8') for crpos in crpos_list])

    crout = sampleName + '_clipped_read_pos.bed.gz'
    f = gzip.open(crout, 'wb')
    try:
        for l in lines:
            f.write(l)
    finally:
        f.close()


def nanosv_vcf_to_bed(sampleName):
    # Load SV list
    sv_list = read_nanosv_vcf(sampleName)

    plot_ci_dist(sv_list, sampleName)

    # nanoSV & Manta SVs
    # sv_list = get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName)

    # Select deletions
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL' if sv.chrom == sv.chrom2 if sv.start < sv.end]

    plot_ci_dist(sv_list, sampleName)

    lines = []
    for sv in sv_list:
        lines.append(bytes(sv.chrom + '\t' + str(sv.start + sv.cipos[0]) + '\t' \
                           + str(sv.start + sv.cipos[1] + 1) + '\t' + 'DEL_start:' + sv.chrom + '_' + str(sv.start) \
                           + '\n', 'utf-8'))
        lines.append(bytes(sv.chrom + '\t' + str(sv.end + sv.ciend[0]) + '\t' \
                           + str(sv.end + sv.ciend[1] + 1) + '\t' + 'DEL_end:' + sv.chrom + '_' + str(sv.end) \
                           + '\n', 'utf-8'))

    outfile = sampleName + '_nanosv_vcf_ci.bed.gz'
    # outfile = sampleName + '_manta_nanosv_vcf_ci.bed.gz'
    f = gzip.open(outfile, 'wb')
    try:
        for l in lines:
            f.write(l)
    finally:
        f.close()


# Get labels
def get_labels(sampleName):
    print(f'running {sampleName}')

    def get_win_id(chr, position):
        return {'chromosome': chr, 'position': position}

    def make_tree_from_bed(sv_list):

        # Using IntervalTree for interval search
        t = IntervalTree()
        # print('# Breakpoints in Chr: %d' % len(sv_list_chr))
        for start, end, lab in sv_list:
            t[start:end + 1] = lab
        return t

    def make_tree_from_vcf(sv_list):

        # Using IntervalTree for interval search
        t = IntervalTree()

        for var in sv_list:
            # cipos[0] and ciend[0] are negative in the VCF file
            id_start = var.svtype + '_start'
            id_end = var.svtype + '_end'

            assert var.start <= var.end, "Start: " + str(var.start) + " End: " + str(var.end)

            # print('var start -> %s:%d CIPOS: (%d, %d)' % (chrName, var.start, var.cipos[0], var.cipos[1]))
            # print('var end -> %s:%d CIEND: (%d, %d)' % (chrName, var.end, var.ciend[0], var.ciend[1]))

            t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
            t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

        return t

    def get_sv_dict():

        hpc_path = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Data_for_labels/', sampleName)
        sv_dict = dict()

        if sampleName in ['NA12878', 'PATIENT1', 'PATIENT2']:
            sv_dict['nanosv'] = read_nanosv_vcf(sampleName)
            sv_dict['nanosv_manta'] = get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName)

            for sv_caller in ['manta', 'delly', 'lumpy', 'gridss']:
                # for sv_caller in ['gridss']:
                sv_dict[sv_caller] = read_vcf(sampleName, sv_caller)

        if sampleName[:7] == 'NA12878':
            # Mills2011
            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/'
            inbed = os.path.join(inbed_path, 'lumpy-Mills2011-call-set.bed')
            sv_dict['Mills2011'] = read_bed_sv(inbed)

            # inbed_path = hpc_path if HPC_MODE else \
            #     '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
            # inbed = os.path.join(inbed_path, 'Mills2011_nanosv_full_inclusion.bed')
            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/'
            inbed = os.path.join(inbed_path, 'lumpy-Mills2011-call-set.nanosv.sorted.bed')
            sv_dict['Mills2011_nanosv'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
            inbed = os.path.join(inbed_path, 'Mills2011_nanosv_full_inclusion.unique.bed')
            sv_dict['Mills2011_nanosv_unique'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
            inbed = os.path.join(inbed_path, 'NA12878_nanosv_Mills2011.bed')
            sv_dict['nanosv_Mills2011'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                os.path.join('/Users/lsantuari/Documents/IGV/Screenshots/', sampleName[:7], 'overlaps')
            inbed = os.path.join(inbed_path, 'lumpy-Mills2011_manta_nanosv.bed')
            sv_dict['Mills2011_nanosv_manta'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                os.path.join('/Users/lsantuari/Documents/Processed/', sampleName[:7], 'Long_read_validation')
            inbed = os.path.join(inbed_path, 'lumpy-Mills2011-DEL.pacbio_moleculo.bed')
            sv_dict['Mills2011_PacBio_Moleculo'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
            inbed = os.path.join(inbed_path, 'Mills2011_pacbio_moleculo_nanosv_full_inclusion.unique.bed')
            sv_dict['Mills2011_PacBio_Moleculo_nanosv'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
            inbed = os.path.join(inbed_path, 'NA12878_nanosv_Mills2011-DEL.pacbio_moleculo.bed')
            sv_dict['nanosv_Mills2011_PacBio_Moleculo'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                os.path.join('/Users/lsantuari/Documents/Processed/', sampleName[:7], 'Long_read_validation')
            inbed = os.path.join(inbed_path, 'lumpy-Mills2011_pacbio_moleculo_manta_nanosv.bed')
            sv_dict['Mills2011_PacBio_Moleculo_nanosv_manta'] = read_bed_sv(inbed)

            inbed_path = hpc_path if HPC_MODE else \
                os.path.join('/Users/lsantuari/Documents/Processed/', sampleName[:7], 'Long_read_validation',
                             'Data_sources', 'Lumpy_paper_2014')
            inbed = os.path.join(inbed_path, 'lumpy-GASVPro-DELLY-Pindel-Mills2011_PacBio_Moleculo.bed')
            sv_dict['Mills2011_PacBio_Moleculo_Lumpy_GASVPro_DELLY_Pindel'] = read_bed_sv(inbed)

        elif sampleName == 'NA24385':
            sv_dict['sv_tier1'] = read_vcf(sampleName, '')

        elif sampleName in ['CHM1', 'CHM13']:
            sv_dict['huddleston2016'] = read_vcf(sampleName, '')

        return sv_dict

    def get_crpos_overlap_with_sv_callsets(sv_dict, cr_pos_dict):

        print(f'Creating crpos_overlap_with_sv_callsets')
        crpos_all_sv = dict()

        for chrName in chrom_lengths.keys():

            print(f'Considering Chr{chrName}')

            # Build two sets: crpos_full_all_sv and crpos_partial_all_sv with clipped read positions that
            # fully/partially overlap at least one SV callset of the caller_list_all_sv
            sv_list_all_sv = dict()
            crpos_full_all_sv_per_caller = dict()
            crpos_partial_all_sv_per_caller = dict()
            caller_list_all_sv = ['manta', 'gridss', 'lumpy', 'delly', 'nanosv']

            for caller in caller_list_all_sv:
                print(caller)
                sv_list_all_sv[caller] = [var for var in sv_dict[caller] if var.chrom == chrName]
                crpos_full_all_sv_per_caller[caller], crpos_partial_all_sv_per_caller[caller] = \
                    get_crpos_win_with_ci_overlap(sv_list_all_sv[caller], cr_pos_dict[chrName])

            crpos_full_all_sv = set()
            crpos_partial_all_sv = set()

            for caller in caller_list_all_sv:
                crpos_full_all_sv = crpos_full_all_sv.union(set(crpos_full_all_sv_per_caller[caller]))
                crpos_partial_all_sv = crpos_partial_all_sv.union(set(crpos_partial_all_sv_per_caller[caller]))

            crpos_all_sv[chrName] = crpos_full_all_sv | crpos_partial_all_sv

        print(f'Finished crpos_overlap_with_sv_callsets')

        return crpos_all_sv

    cr_pos_dict = load_all_clipped_read_positions(sampleName)

    sv_dict = get_sv_dict()

    # Get overlap of candidate positions with all SV breakpoints (all 4 SV callers)
    # crpos_all_sv = get_crpos_overlap_with_sv_callsets(sv_dict, cr_pos_dict)

    labels = dict()

    # for sv_dict_key in sv_dict.keys():
    #    print(sv_dict_key)
    #    print(sv_dict[sv_dict_key])

    for sv_dict_key in sv_dict.keys():
        # for sv_dict_key in ['Mills2011_nanosv']:
        # for sv_dict_key in ['Mills2011_PacBio_Moleculo_Lumpy_GASVPro_DELLY_Pindel']:

        print(f'running {sv_dict_key}')

        labels[sv_dict_key] = {}

        sv_list = sv_dict[sv_dict_key]

        if type(sv_list) is list:
            print('VCF mode')
            # Select deletions (DELs)
            print('%d SVs (all)' % len(sv_list))

            sv_list = [sv for sv in sv_list if sv.svtype == 'DEL']
            print('%d SVs' % len(sv_list))

            # list of chromosomes
            chr_list = set([var.chrom for var in sv_list])

        else:
            print('BED mode')
            chr_list = sv_list.keys()

        for chrName in chr_list:
            # DEBUG
            # for chrName in ['22']:

            print(f'running Chr{chrName}')

            labels[sv_dict_key][chrName] = []

            # Load CR positions, once
            cr_pos = cr_pos_dict[chrName]

            # print('CRPOS:')
            # print(cr_pos)

            # print(type(sv_list))

            # VCF file SVs
            if type(sv_list) is list:

                sv_list_chr = [var for var in sv_list if var.chrom == chrName]
                tree = make_tree_from_vcf(sv_list_chr)

                crpos_full, crpos_partial = get_crpos_win_with_ci_overlap(sv_list_chr, cr_pos)

            # BED file SVs
            else:

                sv_list_chr = sv_list[chrName]
                tree = make_tree_from_bed(sv_list_chr)

                crpos_full, crpos_partial = get_crpos_win_with_bed_overlap(sv_list_chr, cr_pos)

            # print(f'crpos_full = {crpos_full}')
            # print(f'crpos_partial = {crpos_partial}')

            label_search = [sorted(tree[p - win_hlen: p + win_hlen + 1]) for p in cr_pos]

            # print(tree)
            # print(label_search)

            count_zero_hits = 0
            count_multiple_hits = 0

            for elem, pos in zip(label_search, cr_pos):

                if len(elem) == 1:
                    # print(elem)
                    if pos in crpos_full:
                        labels[sv_dict_key][chrName].append(elem[0].data)
                    elif pos in crpos_partial: # or \
                            # pos in crpos_all_sv[chrName] / crpos_full:
                        labels[sv_dict_key][chrName].append('UK')
                    else:
                        labels[sv_dict_key][chrName].append('noSV')
                elif len(elem) == 0:
                    count_zero_hits += 1
                    labels[sv_dict_key][chrName].append('noSV')
                elif len(elem) > 1:
                    count_multiple_hits += 1
                    labels[sv_dict_key][chrName].append('UK')

            assert len(labels[sv_dict_key][chrName]) == len(cr_pos)

            print(Counter(labels[sv_dict_key][chrName]))

    # Add window ID
    labels['id'] = {}
    for chrName in chr_list:
        labels['id'][chrName] = []
        for pos in cr_pos_dict[chrName]:
            labels['id'][chrName].append(get_win_id(chrName, pos))
        assert len(labels['id'][chrName]) == len(cr_pos_dict[chrName])

    # pp = pprint.PrettyPrinter(depth=6)
    # for key in sv_dict:
    #     pp.pprint(sv_dict[key])

    if not HPC_MODE:
        channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
    else:
        channel_dir = '.'

    output_dir = '/'.join((channel_dir, sampleName, 'label_npy'+'_win'+str(win_len)))
    create_dir(output_dir)

    data_file = '/'.join((output_dir, 'labels.pickle'))
    # print(output_dir)
    pickle.dump(labels, open(data_file, "wb"))
    os.system('gzip ' + data_file)


def load_labels(sampleName):
    if not HPC_MODE:
        channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
    else:
        channel_dir = ''

    output_dir = '/'.join((channel_dir, sampleName+'_win'+win_len, 'label_npy'))

    pickle_file = '/'.join((output_dir, 'labels.pickle.gz'))
    with gzip.GzipFile(pickle_file, "rb") as f:
        labels = pickle.load(f)
    f.close()

    print(labels['id'])


def main():
    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    # parser.add_argument('-b', '--bam', type=str,
    #                     default=inputBAM,
    #                     help="Specify input file (BAM)")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output")
    parser.add_argument('-s', '--sample', type=str, default='NA12878',
                        help="Specify sample")

    args = parser.parse_args()

    # bed_dict = dict()
    # for sampleName in ['NA12878', 'Patient1', 'Patient2']:
    #     bed_dict[sampleName] = dict()
    #     bed_dict[sampleName]['NanoSV_Manta'] = '/Users/lsantuari/Documents/IGV/Screenshots/' + sampleName + '/' + \
    #                                            sampleName + '_manta_nanosv_vcf_ci.bed'
    #     if sampleName == 'NA12878':
    #         bed_dict[sampleName]['NanoSV_Manta_Mills2011'] = '/Users/lsantuari/Documents/IGV/Screenshots/' + sampleName \
    #                                                          + '/overlaps/' + 'lumpy-Mills2011_manta_nanosv.bed'

    t0 = time()

    # Iterate over samples and BED files
    # for sampleName in bed_dict.keys():
    #     print('Running get_labels for %s' % sampleName)
    #     get_labels(sampleName)
    # for bed_file in bed_dict[sampleName].keys():
    #     print('Processing BED file: %s' % bed_file)
    #     labels_list = get_labels_from_bed(sampleName=sampleName, ibam=args.bam, inbed=bed_dict[sampleName][bed_file])
    #     for c in labels_list.keys():
    #         print('%s -> %s'% (c, Counter(labels_list[c])))

    # write_sv_without_cr(sampleName=args.sample, ibam=args.bam)

    # clipped_read_positions_to_bed(sampleName=args.sample, ibam=args.bam)
    # nanosv_vcf_to_bed(sampleName=args.sample)

    # get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName=args.sample)
    #
    # for sampleName in ['NA12878']:
    #     nanosv_vcf_to_bed(sampleName)

    # for sampleName in ['NA12878']:
    # get_labels_from_nanosv_vcf(sampleName=sampleName)
    # load_labels(sampleName=sampleName)

    # for sampleName in ['NA12878', 'Patient1', 'Patient2']:
    # for sampleName in ['NA24385', 'CHM1', 'CHM13']:
    for sampleName in ['NA12878']:
        get_labels(sampleName)
        # nanosv_vcf_to_bed(sampleName)

    # crpos_giab = load_all_clipped_read_positions('NA12878')
    # crpos_ena = load_all_clipped_read_positions('NA12878_ENA')
    #
    # print(Counter(crpos_giab))
    # print(Counter(crpos_ena))

    # for chr in crpos_ena.keys():
    # print(chr)
    # print(set(crpos_giab[chr])-set(crpos_ena[chr]))

    print('Elapsed time making labels = %f' % (time() - t0))


if __name__ == '__main__':
    main()
