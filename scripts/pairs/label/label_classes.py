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

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = False

# Window half length
win_hlen = 100
# Window size
win_len = win_hlen * 2

__bpRE__ = None
__symbolicRE__ = None


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
                indellen = abs(record.stop-record.pos)

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
            self.cipos = (0,0)

        # CIEND
        if 'CIEND' in record.info.keys():
            if 'CIEND95' in record.info.keys():
                self.ciend = record.info['CIEND95']
            else:
                self.ciend = record.info['CIEND']
        elif 'CIRPOS' in record.info.keys():
            self.ciend = record.info['CIRPOS']
        else:
            self.ciend = (0,0)

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
                vcf_files[mapper]['nanosv'] = os.path.join(vcf_dir, mapper, mapper+'_nanosv.sorted.vcf')
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
              #if 'PASS' in list(svrec.filter)]
              if 'PASS' in list(svrec.filter) or \
              'CIPOS' in list(svrec.filter) or \
              'CIEND' in list(svrec.filter)]

        return sv
    else:
        return None
