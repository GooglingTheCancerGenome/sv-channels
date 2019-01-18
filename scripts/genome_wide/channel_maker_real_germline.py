# Imports

import re
import pysam
import statistics
from pysam import VariantFile
from collections import Counter
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import numpy as np
import argparse
import bz2file
import gzip
import os, errno
import pickle
from time import time
import logging
import csv
import pandas as pd
from plotnine import *
import pyBigWig
from functions import get_one_hot_sequence, is_outlier
from itertools import chain

# import matplotlib.pyplot as plt

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = True

# Only clipped read positions supported by at least min_cr_support clipped reads are considered
min_cr_support = 3
# Window half length
win_hlen = 100
# Window size
win_len = win_hlen * 2


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

    def __init__(self, record):

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


__bpRE__ = None
__symbolicRE__ = None


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


def get_chr_len(ibam, chrName):
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    return chrLen


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
    channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'

    vec_type = 'clipped_read_pos'
    print('Loading CR positions for Chr %s' % chrName)
    # Load files
    if HPC_MODE:
        fn = '/'.join((sampleName, vec_type, chrName + '_' + vec_type + '.pbz2'))
    else:
        fn = '/'.join((channel_dir, sampleName, vec_type, chrName + '_' + vec_type + '.pbz2'))
    with bz2file.BZ2File(fn, 'rb') as f:
        cpos = pickle.load(f)

    cr_pos = [elem for elem, cnt in cpos.items() if cnt >= min_cr_support]

    return cr_pos


def count_clipped_read_positions(cpos_cnt):
    '''

    :param cpos_cnt: dictionary of clipped read positions (keys) and counts of clipped reads per position (values) as
    returned by the clipped_read_pos.py script
    :return: None. Prints the number of clipped read positions with clipped read support greater than the integers
    specified in the range
    '''
    for i in range(0, 5):
        logging.info('Number of positions with at least %d clipped reads: %d' %
                     (i + 1, len([k for k, v in cpos_cnt.items() if v > i])))


def read_Mills2011_positions():
    '''
    Deprecated. Returns the positions of the deletions and insertions of the GiaB dataset NA12878 from the paper
    Mills et al., 2011
    :return: del_start: dictionary of start positions of deletions per chromosome
             del_end:  dictionary of end positions of deletions per chromosome. Corresponds to the
                       start positions in del_start
             ins_start: dictionary of start positions of insertions per chromosome
    '''
    gold_del_file = '/Users/lsantuari/Documents/Data/GiaB/GoldStandard/DEL_Mills2011.csv'
    gold_insdup_file = '/Users/lsantuari/Documents/Data/GiaB/GoldStandard/INS_DUP_Mills2011.csv'

    del_start = dict()
    del_end = dict()
    ins_start = dict()

    with open(gold_del_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quotechar='#')
        for row in reader:
            # print row
            if row[0] == 'NA12878':
                chr = row[1][3:]
                # chr23 is chrX, rename it properly
                if chr == '23':
                    chr = 'X'
                if chr in del_start.keys():
                    del_start[chr].append(int(row[2]))
                else:
                    del_start[chr] = [int(row[2])]
                if chr in del_end.keys():
                    del_end[chr].append(int(row[2]))
                else:
                    del_end[chr] = [int(row[2])]

    with open(gold_insdup_file, 'r') as csvfile:
        reader = csv.reader(csvfile, delimiter=';', quotechar='#')
        for row in reader:
            # print(row[4])
            if row[0] == 'NA12878' and row[4] == 'INS':
                chr = row[1]
                if chr in ins_start.keys():
                    ins_start[chr].append(int(row[2]))
                else:
                    ins_start[chr] = [int(row[2])]
    # print(del_start.keys())
    return del_start, del_end, ins_start


def read_bpi_positions():
    '''
    Reads Manta/BPI TSV file for the COLO829 HMF dataset
    :return: Two dictionaries: del_start and del_end, with start and end positions of deletions per chromosome
    '''
    gold_del_file = '/Users/lsantuari/Documents/Data/HMF/HMF_COLO829_VCF/bpi/COLO829R_COLO829T_bpi_stats.tsv'

    del_start = dict()
    del_end = dict()

    with open(gold_del_file, 'r') as tsvfile:
        reader = csv.reader(tsvfile, delimiter='\t')
        for row in reader:
            # print(row)
            if row[1] == 'DEL':
                chrA, start = row[3].split(':')
                chrB, end = row[4].split(':')
                # print("%s %s %s %s" % (chrA, start, chrB, end))
                assert chrA == chrB
                if chrA in del_start.keys():
                    del_start[chrA].append(int(start))
                else:
                    del_start[chrA] = [int(start)]
                if chrB in del_end.keys():
                    del_end[chrB].append(int(end))
                else:
                    del_end[chrB] = [int(end)]

    # print(del_start.keys())
    return del_start, del_end


# Generate labels for the HMF COLO829 dataset
def generate_labels():
    '''
    For each chromosome, load the clipped read positions file (clipped_read_pos) and label each position according
    to the informations in the VCF file for the HMF COLO829 dataset
    :return: None
    '''

    # Half of the window length centered on SV positions, where to consider clipped read positions
    confint = 25

    # Keep track of the windows where no clipped reads can be found
    # How many SV positions reported by the SV caller do not contain any clipped reads?
    no_clipped_read_pos = dict()

    # Deletions and insertions of the GiaB dataset
    # del_start, del_end, ins_start = read_Mills2011_positions()

    # consider HMF COLO829 Manta/BPI positions
    # del_start, del_end = read_bpi_positions()

    # read Lumpy positions
    del_start, del_end = get_Lumpy_positions()

    sum_cr = 0

    # Iterate over chromosomes
    for chr in del_start.keys():
        print('Considering chromosome %s' % chr)

        # Local path to clipped read positions
        clipped_read_pos_file = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/clipped_read_pos/' \
                                + chr + '_clipped_read_pos.pbz2'

        logging.info('Reading clipped read positions')
        with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
            clipped_pos_cnt = pickle.load(f)
        logging.info('End of reading')

        # Select clipped read positions with minimum support of min_cr_support clipped reads
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

        # Consider positions referring to start of deletions
        if chr in del_start.keys():
            start_bpj_lst = []
            for pos in set(del_start[chr]):
                # Take the intersection between the position in the [pos-confint,pos+confint] window positions and the
                # clipped read positions
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                # If some clipped read positions are included in the window, add them to the set of candidate breakpoint
                # junctions
                if len(int_set) > 0:
                    start_bpj_lst.extend(int_set)
                else:
                    # If there are no clipped read positions, record the SV position in the no_clipped_read_pos
                    # dictionary
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        # Consider positions referring to end of deletions
        if chr in del_end.keys():
            end_bpj_lst = []
            for pos in set(del_end[chr]):
                # Take the intersection between the position in the [pos-confint,pos+confint] window positions and the
                # clipped read positions
                int_set = set(list(range(pos - confint, pos + confint))) & set(clipped_pos)
                # If some clipped read positions are included in the window, add them to the set of candidate breakpoint
                # junctions
                if len(int_set) > 0:
                    end_bpj_lst.extend(int_set)
                else:
                    # If there are no clipped read positions, record the SV position in the no_clipped_read_pos
                    # dictionary
                    if chr in no_clipped_read_pos.keys():
                        no_clipped_read_pos[chr].append(pos)
                    else:
                        no_clipped_read_pos[chr] = [pos]

        # print(set(start_bpj_lst) & set(end_bpj_lst))
        # assert len(set(start_bpj_lst) & set(end_bpj_lst)) == 0
        labels = []
        for pos in clipped_pos:
            # Is a clipped read position both in the list of start and end candidate positions?
            if pos in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_start_and_end')
            # The clipped read position belongs only to the set of start breakpoint candidates
            elif pos in start_bpj_lst and pos not in end_bpj_lst:
                labels.append('del_start')
            # The clipped read position belongs only to the set of end breakpoint candidates
            elif pos not in start_bpj_lst and pos in end_bpj_lst:
                labels.append('del_end')
            # The clipped read position does not belong to any deletion endpoint positions
            else:
                labels.append('no_del')

        # A label should be assigned to each clipped read position
        assert len(clipped_pos) == len(labels)

        # (Locally) save the labels per chromosome
        outFile = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/labels/' \
                  + chr + '_labels.npy.gz'
        with gzip.GzipFile(outFile, "w") as f:
            np.save(file=f, arr=labels)
        f.close()

    # File with SV positions without clipped reads
    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/200418/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "w") as f:
        pickle.dump(no_clipped_read_pos, f)
    f.close()

    print("Number of clipped positions: %d" % sum_cr)
    print("No clipped positions: %d" % len(list(no_clipped_read_pos.values())))


def get_Lumpy_positions():
    '''
    This function returns two dictionaries: del_start and del_end. They contain respectively start and end positions
    of the deletions (SVTYPE:DEL) present in the Lumpy output file.
    :return: del_start and del_end dictionaries
    '''

    del_start = dict()
    del_end = dict()

    # Local Lumpy VCF output file for the Z424 dataset
    # vcf_file = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/VCF/Z424/Lumpy.vcf'
    # On the HPC
    vcf_file = '/hpc/cog_bioinf/ridder/users/cshneider/Breast_Cancer_Pilot_outside_IAP/' + \
               'Z424/Lumpy/somatic_SVs/Z424/lumpy_Z424.vcf'
    vcf_in = VariantFile(vcf_file)

    for rec in vcf_in.fetch():
        if rec.info['SVTYPE'] == 'DEL' and rec.chrom not in ['Y', 'MT']:
            # print (str(rec.chrom)+' '+str(rec.pos)+' '+str(rec.stop))

            if rec.chrom in del_start.keys():
                del_start[rec.chrom].append(rec.pos)
            else:
                del_start[rec.chrom] = [rec.pos]

            if rec.chrom in del_end.keys():
                del_end[rec.chrom].append(rec.stop)
            else:
                del_end[rec.chrom] = [rec.stop]
    # print(del_start.keys())
    return del_start, del_end


# NanoSV specific functions:


def initialize_nanosv_vcf(sampleName):
    vcf_files = dict()

    if HPC_MODE:

        if sampleName == 'NA12878':

            vcf_dir = '/hpc/cog_bioinf/kloosterman/shared/nanosv_comparison/NA12878'

            for mapper in ['bwa', 'minimap2', 'ngmlr', 'last']:

                vcf_files[mapper] = dict()

                vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv.sorted.vcf'
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

                vcf_files[mapper]['nanosv_sniffles_settings'] = vcf_dir + '/' + mapper + '/' + \
                                                                mapper + '_nanosv_with_sniffles_settings.sorted.vcf'
                assert os.path.isfile(vcf_files[mapper]['nanosv_sniffles_settings'])

                if mapper in ['bwa', 'ngmlr']:
                    vcf_files[mapper]['sniffles'] = vcf_dir + '/' + mapper + '/' + mapper + '_sniffles.sorted.vcf'
                    assert os.path.isfile(vcf_files[mapper]['sniffles'])

                    vcf_files[mapper]['sniffles_nanosv_settings'] = vcf_dir + '/' + mapper + '/' + \
                                                                    mapper + '_sniffles_with_nanosv_settings.sorted.vcf'
                    assert os.path.isfile(vcf_files[mapper]['sniffles_nanosv_settings'])

    else:

        if sampleName == 'NA12878':

            vcf_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth/NA12878/SV'
            vcf_files = dict()

            for mapper in ['bwa', 'last']:

                vcf_files[mapper] = dict()

                vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv_pysam.sorted.vcf'
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

                if mapper in ['bwa']:
                    vcf_files[mapper]['sniffles'] = vcf_dir + '/' + mapper + '/' + mapper + '_sniffles.sorted.vcf'
                    assert os.path.isfile(vcf_files[mapper]['sniffles'])

        elif sampleName == 'Patient1' or sampleName == 'Patient2':

            vcf_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth/' + sampleName + '/SV'
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
    vcf_files = initialize_nanosv_vcf(sampleName)

    if sampleName == 'NA12878' or sampleName == 'Patient1' or sampleName == 'Patient2':

        filename = vcf_files['last']['nanosv']
        vcf_in = VariantFile(filename, 'r')

        sv = []

        # create sv list with SVRecord objects
        for rec in vcf_in.fetch():
            resultBP = re.match(__bpRE__, rec.alts[0])
            if resultBP:
                svrec = SVRecord_nanosv(rec)
                sv.append(svrec)

        # Select good quality (no LowQual) deletions (DEL)
        sv = [svrec for svrec in sv if svrec.svtype == 'DEL'
              # if 'LowQual' not in list(svrec.filter)]
              if 'PASS' in list(svrec.filter)]

        # How many distinct FILTERs?
        filter_set = set([f for svrec in sv for f in svrec.filter])
        # print(filter_set)

        # How many VCF record with a specific FILTER?
        filter_list = sorted(filter_set)
        s = pd.Series([sum(list(map(lambda x: int(f in x),
                                    [list(svrec.filter) for svrec in sv])))
                       for f in filter_list],
                      index=filter_list)
        s = s.append(pd.Series([len(sv)], index=['Total']))
        s = s.sort_values()
        # Print pd.Series with stats on FILTERs
        # print(s)

        return sv


def create_labels_nanosv_vcf(sampleName, ibam):
    '''
    This function writes the label files based on the nanosv VCF file information

    :param sampleName: str, name of the sample considered
    :param chrName: str, chromosome name
    :param ibam: str, path of the BAM file in input
    :return: None
    '''

    # Lines to write in the BED file
    lines = []

    def closest_loc(pos, pos_list):
        pos_array = np.asarray(pos_list)
        deltas = np.abs(pos_array - pos)
        idx = np.argmin(deltas)
        return (pos_list[idx], deltas[idx])

    # Load SV list
    sv_list = read_nanosv_vcf(sampleName)
    # Select deletions
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL']
    # list of chromosomes
    chr_list = set([var.chrom for var in sv_list])

    # print('Plotting CI distribution')
    # plot_ci_dist(sv_list, sampleName)

    print(chr_list)
    print('# DELs:%d' % len(sv_list))

    cnt = Counter([sv.chrom for sv in sv_list])
    s = pd.Series([v for v in cnt.values()], index=cnt.keys())
    print(s)

    assert sum(s) == len(sv_list)

    confint = 100

    for chrName in chr_list:

        sv_list_chr = [var for var in sv_list if var.chrom == chrName]

        chrLen = get_chr_len(ibam, chrName)

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # Remove positions with windows falling off chromosome boundaries
        cr_pos = [pos for pos in cr_pos if win_hlen <= pos <= (chrLen - win_hlen)]

        # print(sorted(cr_pos))

        # Using IntervalTree for interval search
        t = IntervalTree()

        print('# SVs in Chr: %d' % len(sv_list_chr))

        for var in sv_list_chr:
            # cipos[0] and ciend[0] are negative in the VCF file
            id_start = var.svtype + '_start'
            id_end = var.svtype + '_end'

            # id_start = '_'.join((var.chrom, str(var.start+var.cipos[0]),  str(var.start+var.cipos[1])))
            # id_end = '_'.join((var.chrom, str(var.end + var.ciend[0]), str(var.end+var.ciend[1])))
            assert var.start < var.end

            # print('var start -> %s:%d CIPOS: (%d, %d)' % (chrName, var.start, var.cipos[0], var.cipos[1]))
            # print('var end -> %s:%d CIEND: (%d, %d)' % (chrName, var.end, var.ciend[0], var.ciend[1]))

            t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
            t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

            # t[var.start - confint:var.start + confint + 1] = var.svtype + '_start'
            # t[var.end - confint:var.end + confint + 1] = var.svtype + '_end'

        label = [sorted(t[p - win_hlen: p + win_hlen + 1]) for p in cr_pos]

        crpos_full_ci, crpos_partial_ci = get_crpos_win_with_ci_overlap(sv_list_chr, cr_pos)
        print('Clipped read positions with complete CI overlap: %s' % crpos_full_ci)
        print('Clipped read positions with partial CI overlap: %s' % crpos_partial_ci)
        crpos_ci_isec = set(crpos_full_ci) & set(crpos_partial_ci)
        print('Intersection: %s' % crpos_ci_isec)

        print('# CRPOS in CI: %d' % len([l for l in label if len(l) != 0]))

        count_zero_hits = 0
        count_multiple_hits = 0

        label_ci = []
        label_ci_full_overlap = []

        for elem, pos in zip(label, cr_pos):
            if len(elem) == 1:
                # print(elem)
                label_ci.append(elem[0].data)
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
                label_ci.append('noSV')
                label_ci_full_overlap.append('noSV')
            elif len(elem) > 1:
                count_multiple_hits += 1
                label_ci.append('UK')
                label_ci_full_overlap.append('UK')
                # if pos in crpos_full_ci:
                #     label_ci_full_overlap.append('Multiple_Full')
                #     #print('Multiple full: %s -> %s' % ( [d for s,e,d in elem], set([d for s,e,d in elem]) ) )
                #     #for s, e, d in elem:
                #     #    print('%d -> %d %d %s' % (pos, s, e, d))
                # #else:
                #     #label_ci_full_overlap.append('Multiple_Partial')

        print('CR positions: %d' % len(cr_pos))
        print('Label length: %d' % len(label))
        assert len(label_ci) == len(cr_pos)
        assert len(label_ci_full_overlap) == len(cr_pos)

        print('Label_CI: %s' % Counter(label_ci))
        print('Label_CI_full_overlap: %s' % Counter(label_ci_full_overlap))
        print('Zero hits:%d' % count_zero_hits)
        print('Multiple hits:%d' % count_multiple_hits)

        hit_set = set([elem for l in label for elem in l])
        # print(hit_set)
        hits = [(var.start, var.end) for var in sv_list_chr
                if var.start + var.cipos[0] in [start for start, end, lab in hit_set]
                if var.start + var.cipos[1] + 1 in [end for start, end, lab in hit_set]
                if var.end + var.ciend[0] in [start for start, end, lab in hit_set]
                if var.end + var.ciend[1] + 1 in [end for start, end, lab in hit_set]]
        # print(hits)
        # no_hits = sorted(no_hits)

        # print('No hits:%d/%d' % (len(sv_list_chr) - len(hits), len(sv_list_chr)))
        # print(no_hits)

        # Pick the clipped read position closest to either var.start or var.end
        label_close = []
        label_BPJ = []
        label_distance = []

        start_positions = [var.start for var in sv_list_chr]
        end_positions = [var.end for var in sv_list_chr]

        # Very slow:
        if len(start_positions) != 0 and len(end_positions) != 0:
            for pos in cr_pos:
                # print('Pos: %d' % pos)
                close_start, close_start_dist = closest_loc(pos, start_positions)
                close_end, close_end_dist = closest_loc(pos, end_positions)
                if close_start_dist <= confint or close_end_dist <= confint:
                    if close_start_dist <= close_end_dist:
                        # print(str(close_start_dist))
                        bpj_id = str(chrName) + '_' + str(close_start)
                        # print(bpj_id)
                        label_close.append('DEL_start')
                        label_BPJ.append(bpj_id)
                        label_distance.append(abs(pos - close_start))
                    else:
                        # print(str(close_end_dist))
                        bpj_id = str(chrName) + '_' + str(close_end)
                        # print(bpj_id)
                        label_close.append('DEL_end')
                        label_BPJ.append(bpj_id)
                        label_distance.append(abs(pos - close_end))
                else:
                    label_close.append('noSV')
                    label_BPJ.append('noSV')
                    label_distance.append('noSV')

            assert len(cr_pos) == len(label_close)

            # print(len(labels_list))
            print('Chr:%s' % chrName)
            print(Counter(map(len, label)))
            print(Counter(label_close))

            if not HPC_MODE:
                channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
            else:
                channel_dir = ''

            output_dir = '/'.join((channel_dir, sampleName, 'label'))
            create_dir(output_dir)

            # print(output_dir)

            with gzip.GzipFile('/'.join((output_dir, chrName + '_label_ci.npy.gz')), "w") as f:
                np.save(file=f, arr=label_ci)
            f.close()
            with gzip.GzipFile('/'.join((output_dir, chrName + '_label_ci_full_overlap.npy.gz')), "w") as f:
                np.save(file=f, arr=label_ci_full_overlap)
            f.close()
            with gzip.GzipFile('/'.join((output_dir, chrName + '_label.npy.gz')), "w") as f:
                np.save(file=f, arr=label_close)
            f.close()
            with gzip.GzipFile('/'.join((output_dir, chrName + '_label_BPJ.npy.gz')), "w") as f:
                np.save(file=f, arr=label_BPJ)
            f.close()
            with gzip.GzipFile('/'.join((output_dir, chrName + '_label_distance.npy.gz')), "w") as f:
                np.save(file=f, arr=label_distance)
            f.close()

        else:
            print('Chromosome %s is not in the SV file' % chrName)

    # Write BED file with labelled CI positions
    outfile = sampleName + '_nanosv_vcf_ci_labelled.bed.gz'
    f = gzip.open(outfile, 'wb')
    try:
        for l in set(lines):
            f.write(l)
    finally:
        f.close()


def write_sv_without_cr(sampleName, ibam):
    '''

    :param sampleName:
    :param ibam:
    :return:
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

        chrLen = get_chr_len(ibam, chrName)

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)
        # Remove positions with windows falling off chromosome boundaries
        cr_pos = [pos for pos in cr_pos if win_hlen <= pos <= (chrLen - win_hlen)]

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


def get_crpos_win_with_ci_overlap(sv_list, cr_pos):
    '''

    :param sv_list: list, list of SVs
    :param cr_pos: list, list of clipped read positions
    :return: list, list of clipped read positions whose window completely overlap either the CIPOS interval
    or the CIEND interval
    '''
    # Tree with windows for CR positions
    t_cr = IntervalTree()
    for pos in cr_pos:
        t_cr[pos - win_hlen:pos + win_hlen + 1] = pos
        # t_cr[pos - 100:pos + 100 + 1] = pos

    cr_full_overlap_cipos = []
    cr_partial_overlap_cipos = []

    rg_overlap = [sorted(t_cr[var.start + var.cipos[0]: var.start + var.cipos[1] + 1]) for var in sv_list]
    # print('Range overlap: %s' % rg_overlap)

    for rg, start, end in zip(rg_overlap,
                              [var.start + var.cipos[0] for var in sv_list],
                              [var.start + var.cipos[1] for var in sv_list]):
        for elem in rg:
            elem_start, elem_end, elem_data = elem
            if start >= elem_start and end <= elem_end:
                # print('CIPOS->Full: %s\t%d\t%d' % (elem, start, end))
                cr_full_overlap_cipos.append(elem_data)
            else:
                # print('CIPOS->Partial: %s\t%d\t%d' % (elem, start, end))
                cr_partial_overlap_cipos.append(elem_data)

    cr_full_overlap_ciend = []
    cr_partial_overlap_ciend = []

    rg_overlap = [sorted(t_cr[var.end + var.ciend[0]: var.end + var.ciend[1] + 1]) for var in sv_list]
    # print('Range overlap: %s' % rg_overlap)

    for rg, start, end in zip(rg_overlap,
                              [var.end + var.ciend[0] for var in sv_list],
                              [var.end + var.ciend[1] for var in sv_list]):
        for elem in rg:
            elem_start, elem_end, elem_data = elem
            if start >= elem_start and end <= elem_end:
                # print('CIEND->Full: %s\t%d\t%d' % (elem, start, end))
                cr_full_overlap_ciend.append(elem_data)
            else:
                # print('CIEND->Partial: %s\t%d\t%d' % (elem, start, end))
                cr_partial_overlap_ciend.append(elem_data)

    # print('Intersection CIPOS: %s' % list(set(cr_full_overlap_cipos) & set(cr_partial_overlap_cipos)))
    # print('Intersection CIEND: %s' % list(set(cr_full_overlap_ciend) & set(cr_partial_overlap_ciend)))

    cr_full_overlap = sorted(cr_full_overlap_cipos + cr_full_overlap_ciend)
    cr_partial_overlap = sorted(cr_partial_overlap_cipos + cr_partial_overlap_ciend)

    return sorted(list(set(cr_full_overlap))), sorted(list(set(cr_partial_overlap) - set(cr_full_overlap)))


def plot_ci_dist(sv_list, sampleName):
    '''
    Save the plots of the distributions of the confidence intervals reported by NanoSV
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
        # 0:Delly, 1:GRIDSS, 2:nanosv, 3:Lumpy, 4:Manta
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

def read_bed_sv(sampleName):
    # lumpy-Mills2011_manta_nanosv
    # inbed = '/Users/lsantuari/Documents/IGV/Screenshots/NA12878/overlaps/lumpy-Mills2011_manta_nanosv.bed'

    # manta_nanosv
    inbed = '/Users/lsantuari/Documents/IGV/Screenshots/' + sampleName + '/' + sampleName + '_manta_nanosv_vcf_ci.bed'

    assert os.path.isfile(inbed)
    sv_dict = defaultdict(list)
    with(open(inbed, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom = str(columns[0])
            if columns[3][:3] == "DEL":
                sv_dict[chrom].append((int(columns[1]), int(columns[2]), columns[3]))

    # print(sv_dict)
    return sv_dict


def create_labels_bed(sampleName, ibam):
    print('sample = %s' % sampleName)
    print('window = %d' % win_len)

    sv_list = read_bed_sv(sampleName)
    chr_list = sv_list.keys()

    for chrName in chr_list:

        sv_list_chr = sv_list[chrName]

        chrLen = get_chr_len(ibam, chrName)

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # Remove positions with windows falling off chromosome boundaries
        cr_pos = [pos for pos in cr_pos if win_hlen <= pos <= (chrLen - win_hlen)]

        # print(sorted(cr_pos))

        # Using IntervalTree for interval search
        t = IntervalTree()

        print('# Breakpoints in Chr: %d' % len(sv_list_chr))

        for start, end, lab in sv_list_chr:
            t[start:end + 1] = lab

        label = [sorted(t[p - win_hlen: p + win_hlen + 1]) for p in cr_pos]

        crpos_full_ci, crpos_partial_ci = get_crpos_win_with_bed_overlap(sv_list_chr, cr_pos)
        print('Clipped read positions with complete CI overlap: %s' % crpos_full_ci)
        print('Clipped read positions with partial CI overlap: %s' % crpos_partial_ci)
        crpos_ci_isec = set(crpos_full_ci) & set(crpos_partial_ci)
        print('Intersection: %s' % crpos_ci_isec)

        print('# CRPOS in BED: %d' % len([l for l in label if len(l) != 0]))

        count_zero_hits = 0
        count_multiple_hits = 0

        label_ci_full_overlap = []

        for elem, pos in zip(label, cr_pos):
            if len(elem) == 1:
                # print(elem)
                if pos in crpos_full_ci:
                    label_ci_full_overlap.append(elem[0].data)
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

        print('CR positions: %d' % len(cr_pos))
        print('Label length: %d' % len(label))
        assert len(label_ci_full_overlap) == len(cr_pos)

        print('Label_CI_full_overlap: %s' % Counter(label_ci_full_overlap))
        print('Zero hits:%d' % count_zero_hits)
        print('Multiple hits:%d' % count_multiple_hits)

        if not HPC_MODE:
            channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
        else:
            channel_dir = ''

        output_dir = '/'.join((channel_dir, sampleName, 'label'))
        create_dir(output_dir)

        # print(output_dir)

        with gzip.GzipFile('/'.join((output_dir, chrName + '_label_ci_full_overlap.npy.gz')), "w") as f:
            np.save(file=f, arr=label_ci_full_overlap)
        f.close()

    else:
        print('Chromosome %s is not in the SV file' % chrName)


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
        # t_cr[pos - 100:pos + 100 + 1] = pos

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
    if HPC_MODE:
        # To fill
        survivor_vcf = ''
    else:
        if sampleName == 'NA12878':
            context = 'trio'
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


def create_labels_from_SURVIVOR_merge_VCF(sampleName):
    if HPC_MODE:
        channel_dir = ''
        vcf_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Breast_Cancer_Pilot_SV/Manta/survivor_merge.vcf'
    else:
        channel_dir = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/ChannelMaker'
        vcf_file = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/SV/survivor_merge.vcf'

    vec = 'clipped_read_pos'
    confint = 100

    vcf_in = VariantFile(vcf_file)
    samples_list = list((vcf_in.header.samples))
    samples = [w.split('_')[0].split('/')[1] for w in samples_list]

    # sv_bedtool = BedTool(vcf_file)
    # print(sv_bedtool)

    sv = []

    # create sv list with SVRecord_SUR objects
    for rec in vcf_in.fetch():
        # print(rec)
        # print(rec.info['SUPP_VEC'][16])
        # avoid SVs on chromosomes Y and MT
        if rec.chrom not in ['Y', 'MT'] and rec.info['CHR2'] not in ['Y', 'MT']:
            # print(rec)
            sv.append(SVRecord_SUR(rec))

    svtypes = set(var.svtype for var in sv)

    # for t in svtypes:
    #    for sample_idx in range(len(sv[0].supp_vec)):
    #        pass
    # print([var.supp_vec[0][0] for var in sv if var.svtype == t])
    # print(
    #      'Sample %s, SVTYPE:%s = %d' % (
    #      sample_idx, t,
    #      sum([int(var.supp_vec[sample_idx]) for var in sv if var.svtype == t])
    #      )
    # )

    # chromosome to consider
    # chr = '2'
    # SVTYPE
    # type = 'DEL'

    cpos = dict()
    for tn in ['Tumor', 'Normal']:
        cpos[tn] = dict()

    # for chr in set([var.chrom for var in sv]):
    #
    #     print('Loading CR positions for Chr %s' % chr)
    #     # Load files
    #     fn = '/'.join((channel_dir, sampleName, 'Tumor', vec, chr + '_' + vec + '.pbz2'))
    #     with bz2file.BZ2File(fn, 'rb') as f:
    #         cpos['Tumor'][chr] = pickle.load(f)
    #     fn = '/'.join((channel_dir, sampleName, 'Normal', vec, chr + '_' + vec + '.pbz2'))
    #     with bz2file.BZ2File(fn, 'rb') as f:
    #         cpos['Normal'][chr] = pickle.load(f)

    chr_list = set([var.chrom for var in sv])

    # for chr in chr_list:
    #
    #     print('Considering Chr %s' % chr)
    #
    #     clipped_pos = [k for k, v in cpos['Tumor'][chr].items() if v >= min_CR_support]
    #
    #     #for cr in clipped_pos:
    #     #    for type in set([var.svtype for var in sv]):
    #
    #     # Positions with at least cnt support
    #     cpos_set1 = set([elem for elem, cnt in cpos['Tumor'][chr].items() if cnt >= min_CR_support])
    #     cpos_set2 = set([elem for elem, cnt in cpos['Normal'][chr].items() if cnt >= 1])
    #
    #     set_diff = cpos_set1 - cpos_set2
    #     set_int = cpos_set1 & cpos_set2
    #
    #     print('Tumor: %d, Normal: %d' % (len(cpos_set1), len(cpos_set2)))
    #     print('Difference clipped reads Tumor - Normal: %d' % len(set_diff))
    #     print('Intersection clipped reads Tumor/Normal: %d' % len(set_int))
    #
    #     for type in svtypes:
    #
    #         sv_pos = [var.start for var in sv if var.chrom == chr
    #                   and var.svtype == type and var.supp_vec[16] == '1']
    #
    #         cnt_hist = []
    #         crpos_cnt_per_int = []
    #         for i in range(confint):
    #             crpos_cnt_per_int = [len(set(range(p - i, p + i)) & cpos_set1) for p in sv_pos]
    #             zero_cnt = [cnt for elem, cnt in Counter(crpos_cnt_per_int).items() if elem == 0]
    #             if len(zero_cnt) > 0:
    #                 cnt_hist.append(zero_cnt[0])
    #         print([p for p, c in zip(sv_pos, crpos_cnt_per_int) if c == 0])
    #
    #         #print(str(Counter(crpos_cnt_per_int).most_common(3)))
    #
    #         df = pd.DataFrame({"Count": cnt_hist})
    #
    #         print('# SV %s:%d' % (type, len(sv_pos)))
    #         print(df[-1:])

    for chr in chr_list:

        print('Loading CR positions for Chr %s' % chr)
        # Load files
        if HPC_MODE:
            fn = '/'.join((sampleName, 'Tumor', vec, chr + '_' + vec + '.pbz2'))
        else:
            fn = '/'.join((channel_dir, sampleName, 'Tumor', vec, chr + '_' + vec + '.pbz2'))
        with bz2file.BZ2File(fn, 'rb') as f:
            cpos['Tumor'][chr] = pickle.load(f)

        if HPC_MODE:
            fn = '/'.join((sampleName, 'Normal', vec, chr + '_' + vec + '.pbz2'))
        else:
            fn = '/'.join((channel_dir, sampleName, 'Normal', vec, chr + '_' + vec + '.pbz2'))
        with bz2file.BZ2File(fn, 'rb') as f:
            cpos['Normal'][chr] = pickle.load(f)

        cr_pos = [elem for elem, cnt in cpos['Tumor'][chr].items() if cnt >= min_cr_support]
        labels = []

        # for type in svtypes:
        #    pass

        # Using IntervalTree for interval search
        t = IntervalTree()

        if HPC_MODE:
            sample_for_index = sampleName.split('/')[1]
        else:
            sample_for_index = sampleName

        for var in sv:
            if var.supp_vec[samples.index(sample_for_index)] == '1':
                if var.chrom == chr:
                    t[var.start - confint:var.start + confint] = var.svtype + '_start'
                if var.chrom2 == chr:
                    t[var.end - confint:var.end + confint] = var.svtype + '_end'
        # print(t)
        # print(len(cr_pos))
        labels_list = [sorted(t[p]) for p in cr_pos]

        # print(Counter(map(len, labels_list)))

        def test_noSV(x):
            if len(x) == 0:
                return ['noSV']
            else:
                my_list = []
                for elem in x:
                    begin, end, data = elem
                    my_list.append(data)
                return my_list

        label = list(map(test_noSV, labels_list))
        print(Counter([x for elem in label for x in elem]))
        output_dir = '/'.join((sampleName, 'label'))
        create_dir(output_dir)

        with gzip.GzipFile('/'.join((output_dir, chr + '_label.npy.gz')), "w") as f:
            np.save(file=f, arr=label)
        f.close()

        # print([lambda x: ['noSV'] if len(x)==0 else x for x in labels_list])

        # for p in cr_pos:
        #     my_label = []
        #     for var in sv:
        #         if var.chrom == chr and var.svtype == 'DEL' and var.supp_vec[16] == '1':
        #
        #             if p in range(var.start - confint, var.start + confint):
        #                 print('%s_start at pos %d in range (%d,%d)' % (var.svtype, p,
        #                                                                 var.start - confint,
        #                                                                 var.start + confint))
        #                 my_label.append(var.svtype+'_start')
        #             if p in range(var.end - confint, var.end + confint):
        #                 print('%s_end at pos %d in range (%d,%d)' % (var.svtype, p,
        #                                                                 var.end - confint,
        #                                                                 var.end + confint))
        #                 my_label.append(var.svtype + '_end')
        #     labels.append(my_label)
        # print(labels)


def load_NoCR_positions():
    '''
    This function provides an overview of SV positions without clipped read support that are stored in the
    no_clipped_read_positions file.
    :return: Noneq
    '''

    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "r") as f:
        no_clipped_read_pos = pickle.load(f)
    f.close()
    print(list(no_clipped_read_pos))


def clipped_read_positions_to_bed(sampleName):
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
    # sv_list = read_nanosv_vcf(sampleName)
    # nanoSV & Manta SVs
    sv_list = get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName)

    # Select deletions
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL' if sv.chrom == sv.chrom2 if sv.start < sv.end]

    lines = []
    for sv in sv_list:
        lines.append(bytes(sv.chrom + '\t' + str(sv.start + sv.cipos[0]) + '\t' \
                           + str(sv.start + sv.cipos[1] + 1) + '\t' + 'DEL_start' + '\n', 'utf-8'))
        lines.append(bytes(sv.chrom + '\t' + str(sv.end + sv.ciend[0]) + '\t' \
                           + str(sv.end + sv.ciend[1] + 1) + '\t' + 'DEL_end' + '\n', 'utf-8'))

    # outfile = sampleName + '_nanosv_vcf_ci.bed.gz'
    outfile = sampleName + '_manta_nanosv_vcf_ci.bed.gz'
    f = gzip.open(outfile, 'wb')
    try:
        for l in lines:
            f.write(l)
    finally:
        f.close()


def get_gc_bigwig():
    bw = pyBigWig.open("/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/UCSC/hg19/hg19.gc5Base.bw")
    return bw


def get_mappability_bigwig():
    bw = pyBigWig.open("/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Mappability/GRCh37.151mer.bw")
    return bw


def channel_maker(ibam, chrList, sampleName, SVmode, trainingMode, outFile):
    '''
    This function loads the channels and for each clipped read position with at least min_cr_support clipped reads, it
    creates a vstack with 22 channel vectors (11 for the Tumor and 11 for the Normal sample) of width equal to
    twice the window length (win_hlen*2).
    :param ibam: BAM file used only to get the length of the chromosome from the header
    :param chrName: chromosome to consider
    :param sampleName: name of the sample pair (HMF, GiaB, or other)
    :param trainingMode: set to True only if it is used to generate the NoSV category for the Training data
    :param outFile: main output file for the list of channel vstacks
    :return: None. It saves a list of channel vstacks. If in 'NoSV' mode, it also saves a list of labels.
    '''

    # List where to store the channel vstacks
    ch_list = []
    # check if the BAM file exists
    assert os.path.isfile(ibam)
    # open the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    # Extract chromosome length from the BAM header
    header_dict = bamfile.header

    chrLen = dict()
    for chrName in chrList:
        chrLen[chrName] = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Get GC BigWig
    bw_gc = get_gc_bigwig()
    # Get Mappability BigWig
    bw_map = get_mappability_bigwig()

    # Set the correct prefix for the path
    if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
        prefix_train = 'Training_' + SVmode + '/'
    else:
        # prefix_train = ''
        # only for ovarian cancer
        prefix_train = 'OC/'

    # Check for file existence
    if not HPC_MODE:

        # Local BAM file for testing
        sample_list = ['T0']

        clipped_read_pos_file = 'clipped_read_pos.pbz2'
        clipped_read_distance_file = 'clipped_read_distance.pbz2'
        clipped_reads_file = 'clipped_reads.pbz2'
        coverage_file = 'coverage.npy.gz'
        split_read_distance_file = 'split_read_distance.pbz2'

        # Check file existence
        assert os.path.isfile(clipped_read_pos_file)
        assert os.path.isfile(clipped_read_distance_file)
        assert os.path.isfile(clipped_reads_file)
        assert os.path.isfile(coverage_file)
        assert os.path.isfile(split_read_distance_file)

    else:

        # Consider a single sample
        sample_list = sampleName.split('_')

        clipped_read_pos_file = dict()
        clipped_read_distance_file = dict()
        clipped_reads_file = dict()
        coverage_file = dict()
        split_read_distance_file = dict()
        clipped_pos_cnt = dict()

        for chrName in chrList:

            # File with clipped read positions, output of the clipped_read_pos script
            # clipped_read_pos_file[chrName] = prefix_train + sample_list[0] + \
            #                                  '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
            clipped_read_pos_file[chrName] = 'clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
            # File with the clipped read distances, output of the clipped_read_distance script
            clipped_read_distance_file[chrName] = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
            # File with the clipped reads, output of the clipped_reads script
            clipped_reads_file[chrName] = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
            # File with the coverage array, output of the coverage script
            coverage_file[chrName] = 'coverage/' + chrName + '_coverage.npy.bz2'
            # File with the split reads and split read distance, output of the split_read_distance script
            split_read_distance_file[chrName] = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

            # Check file existence
            # print('Checking file: %s' % clipped_read_pos_file[chrName])
            # assert os.path.isfile(clipped_read_pos_file[chrName])

            for sample in sample_list:
                # Check file existence
                print('Checking file: %s => %s' % (sample, clipped_read_pos_file[chrName]))
                assert os.path.isfile(prefix_train + sample + '/' + clipped_read_pos_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + clipped_read_distance_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + clipped_reads_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + coverage_file[chrName])
                assert os.path.isfile(prefix_train + sample + '/' + split_read_distance_file[chrName])

            logging.info('Chromosome %s' % chrName)

            if len(sample_list) == 2:

                clipped_pos_cnt_per_sample = dict()
                clipped_pos = dict()

                for sample in sample_list:

                    logging.info('Reading clipped read positions for sample %s' % sample)
                    with bz2file.BZ2File(prefix_train + sample + '/' +
                                         clipped_read_pos_file[chrName], 'rb') as f:
                        clipped_pos_cnt_per_sample[sample] = pickle.load(f)
                    logging.info('End of reading')
                    print('Length of clipped_pos_cnt_per_sample for sample %s: %d' % (sample,
                                                                                      len(clipped_pos_cnt_per_sample[
                                                                                              sample])))

                    # Count the number of clipped read positions with a certain minimum number of clipped reads
                    count_clipped_read_positions(clipped_pos_cnt_per_sample[sample])

                    # if sample == sample_list[0]:
                    #     cr_support = min_cr_support
                    # else:
                    #     cr_support = 1

                    clipped_pos[sample] = [k for k, v in clipped_pos_cnt_per_sample[sample].items()
                                           if v >= min_cr_support]
                    print('Length of clipped_pos_cnt_per_sample ' +
                          ' for sample %s after min support = %d: %d' %
                          (sample, min_cr_support, len(clipped_pos[sample])))

                clipped_pos_keep = set(clipped_pos[sample_list[0]]) - set(clipped_pos[sample_list[1]])
                print('Length of cr_pos_keep: %d' % len(clipped_pos_keep))

                sample = sample_list[0]
                # clipped_pos_cnt[chrName] = {k: v for (k, v) in clipped_pos_cnt_per_sample[sample_list[0]]
                #                             if k in clipped_pos_keep}
                print('Length of clipped_pos_cnt keys: %d, intersection size: %d' %
                      (len(clipped_pos_cnt_per_sample[sample].keys()),
                       len(set(clipped_pos_cnt_per_sample[sample].keys()) & clipped_pos_keep)))

                clipped_pos_cnt[chrName] = {k: clipped_pos_cnt_per_sample[sample][k]
                                            for k in clipped_pos_cnt_per_sample[sample].keys()
                                            if k in clipped_pos_keep}

                # Count the number of clipped read positions with a certain minimum number of clipped reads
                logging.info('Clipped read positions with support only in the Tumor:')
                count_clipped_read_positions(clipped_pos_cnt[chrName])

                clipped_pos[sample] = [pos for pos in clipped_pos[sample]
                                       if win_hlen <= pos <= (chrLen[chrName] - win_hlen)]
                print('Length of cr_pos for sample %s after extremes removed: %d' % (sample,
                                                                                     len(clipped_pos[sample])))

            else:

                logging.info('Reading clipped read positions')
                with bz2file.BZ2File(clipped_read_pos_file[chrName], 'rb') as f:
                    clipped_pos_cnt[chrName] = pickle.load(f)
                logging.info('End of reading')

                # Count the number of clipped read positions with a certain minimum number of clipped reads
                count_clipped_read_positions(clipped_pos_cnt[chrName])

    # Load channel data
    # Dictionaries where to load the channel data
    clipped_read_distance = dict()
    clipped_reads = dict()
    clipped_reads_inversion = dict()
    clipped_reads_duplication = dict()
    clipped_reads_translocation = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    outliers = dict()

    for sample in sample_list:

        prefix = prefix_train + sample + '/' if HPC_MODE else ''

        clipped_pos = dict()

        clipped_read_distance[sample] = dict()
        clipped_reads[sample] = dict()
        clipped_reads_inversion[sample] = dict()
        clipped_reads_duplication[sample] = dict()
        clipped_reads_translocation[sample] = dict()
        coverage[sample] = dict()
        split_reads[sample] = dict()
        split_read_distance[sample] = dict()

        outliers[sample] = dict()

        logging.info('Considering %s' % sample)

        for chrName in chrList:
            logging.info('Considering %s' % chrName)

            logging.info('Reading clipped read distances')
            with bz2file.BZ2File(prefix + clipped_read_distance_file[chrName], 'rb') as f:
                clipped_read_distance[sample][chrName] = pickle.load(f)
            logging.info('End of reading')

            logging.info('Reading clipped reads')
            with bz2file.BZ2File(prefix + clipped_reads_file[chrName], 'rb') as f:
                clipped_reads[sample][chrName], clipped_reads_inversion[sample][chrName], \
                clipped_reads_duplication[sample][chrName], clipped_reads_translocation[sample][chrName] = pickle.load(
                    f)
            logging.info('End of reading')

            logging.info('Reading coverage')
            with bz2file.BZ2File(prefix + coverage_file[chrName], 'rb') as f:
                coverage[sample][chrName] = np.load(file=f)
            logging.info(
                'End of reading, coverage length: %d out of %d' % (len(coverage[sample][chrName]), chrLen[chrName]))

            logging.info('Reading split read distances')
            with bz2file.BZ2File(prefix + split_read_distance_file[chrName], 'rb') as f:
                split_read_distance[sample][chrName], split_reads[sample][chrName] = pickle.load(f)
            logging.info('End of reading')

            logging.info('Finding outliers')
            outliers[sample][chrName] = dict()
            for direction in ['forward', 'reverse']:
                outliers[sample][chrName][direction] = dict()
                for clipped_arrangement in ['left', 'right', 'unclipped']:
                    points = np.array(list(chain.from_iterable(
                        clipped_read_distance[sample][chrName][direction][clipped_arrangement].values()
                    )))
                    outlier_vec = is_outlier(points)
                    outliers[sample][chrName][direction][clipped_arrangement] = \
                        set(points[np.where(outlier_vec)].flatten())
                    # print(outliers)
            logging.info('Outliers found')

    for chrName in chrList:

        # If in 'NoSV' mode, consider all the clipped read positions (minimum clipped read support equal to 0)
        if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
            clipped_pos[chrName] = [k for k, v in clipped_pos_cnt[chrName].items()]
        else:
            clipped_pos[chrName] = [k for k, v in clipped_pos_cnt[chrName].items() if v >= min_cr_support]

    # print(clipped_pos)

    # Dictionaries where to store the channel arrays as generated from the dictionaries
    clipped_read_distance_array = dict()
    clipped_read_distance_num = dict()
    clipped_read_distance_median = dict()
    clipped_read_distance_outlier = dict()

    clipped_reads_array = dict()
    clipped_reads_inversion_array = dict()
    clipped_reads_duplication_array = dict()
    clipped_reads_translocation_array = dict()

    coverage_array = dict()

    split_read_distance_array = dict()
    split_read_distance_num = dict()
    split_read_distance_median = dict()

    split_reads_array = dict()

    gc_array = dict()
    mappability_array = dict()

    # Log info every n_r times
    n_r = 10 ** 3
    # print(n_r)
    last_t = time()

    # See how to proceed from here

    total_clipped_pos = []
    total_chr_for_clipped_pos = []

    for chrName in chrList:
        total_clipped_pos.extend(clipped_pos[chrName])
        total_chr_for_clipped_pos.extend([chrName] * len(clipped_pos[chrName]))
    print('Count of total_chr_for_clipped_pos: %s' % Counter(total_chr_for_clipped_pos))

    for i, outzipped in enumerate(zip(total_chr_for_clipped_pos, total_clipped_pos), start=1):

        chrName = outzipped[0]
        center_pos = outzipped[1]

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d clipped positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Check if center_pos is within the chromosome boundaries
        if win_hlen <= center_pos <= (chrLen[chrName] - win_hlen):

            start_win = center_pos - win_hlen
            end_win = center_pos + win_hlen

            # Build arrays for the numpy vstack
            for sample in sample_list:

                # clipped read distance
                clipped_read_distance_array[sample] = dict()
                clipped_read_distance_num[sample] = dict()
                clipped_read_distance_median[sample] = dict()
                clipped_read_distance_outlier[sample] = dict()

                for direction in ['forward', 'reverse']:
                    clipped_read_distance_array[sample][direction] = dict()
                    clipped_read_distance_num[sample][direction] = dict()
                    clipped_read_distance_median[sample][direction] = dict()
                    clipped_read_distance_outlier[sample][direction] = dict()

                for direction in ['forward', 'reverse']:
                    # for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                    for clipped_arrangement in ['left', 'right', 'unclipped']:
                        clipped_read_distance_array[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                       dtype=int)
                        clipped_read_distance_num[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                     dtype=int)
                        clipped_read_distance_median[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                        dtype=int)
                        clipped_read_distance_outlier[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                         dtype=int)

                        for pos in range(start_win, end_win):
                            if pos in clipped_read_distance[sample][chrName][direction][clipped_arrangement].keys():
                                clipped_read_distance_array[sample][direction][clipped_arrangement][pos - start_win] = \
                                    sum(clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                clipped_read_distance_num[sample][direction][clipped_arrangement][pos - start_win] = \
                                    len(clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                clipped_read_distance_median[sample][direction][clipped_arrangement][pos - start_win] = \
                                    statistics.median(
                                        clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])

                                if direction == 'forward':

                                    clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                    (pos - start_win):] = \
                                        clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                        (pos - start_win):] + \
                                        len(set(
                                            clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                            & outliers[sample][chrName][direction][clipped_arrangement])

                                elif direction == 'reverse':

                                    clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                    :(pos - start_win)] = \
                                        clipped_read_distance_outlier[sample][direction][clipped_arrangement][
                                        :(pos - start_win)] + \
                                        len(set(
                                            clipped_read_distance[sample][chrName][direction][clipped_arrangement][pos])
                                            & outliers[sample][chrName][direction][clipped_arrangement])

                        # print(clipped_read_distance_array[direction][clipped_arrangement])

                # clipped reads
                clipped_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    clipped_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads[sample][chrName][split_direction].keys():
                            clipped_reads_array[sample][split_direction][pos - start_win] = \
                                clipped_reads[sample][chrName][split_direction][pos]

                # clipped reads inversions
                clipped_reads_inversion_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_inversion_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_inversion[sample][chrName][mate_position].keys():
                            clipped_reads_inversion_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_inversion[sample][chrName][mate_position][pos]

                # clipped reads duplication
                clipped_reads_duplication_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_duplication_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_duplication[sample][chrName][mate_position].keys():
                            clipped_reads_duplication_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_duplication[sample][chrName][mate_position][pos]

                # clipped reads traslocation
                clipped_reads_translocation_array[sample] = dict()
                for orientation in ['opposite', 'same']:
                    clipped_reads_translocation_array[sample][orientation] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_translocation[sample][chrName][orientation].keys():
                            clipped_reads_translocation_array[sample][orientation][pos - start_win] = \
                                clipped_reads_translocation[sample][chrName][orientation][pos]

                # coverage
                coverage_array[sample] = coverage[sample][chrName][start_win:end_win]
                assert len(coverage_array[sample]) == win_len

                # split read distance
                split_read_distance_array[sample] = dict()
                split_read_distance_num[sample] = dict()
                split_read_distance_median[sample] = dict()

                for split_direction in ['left', 'right']:
                    split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    split_read_distance_num[sample][split_direction] = np.zeros(win_len, dtype=int)
                    split_read_distance_median[sample][split_direction] = np.zeros(win_len, dtype=int)

                    if pos in split_read_distance[sample][chrName][split_direction].keys():
                        split_read_distance_array[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][chrName][split_direction][pos])
                        split_read_distance_num[sample][split_direction][pos - start_win] = \
                            len(split_read_distance[sample][chrName][split_direction][pos])
                        split_read_distance_median[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][chrName][split_direction][pos])

                # split reads
                split_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in split_reads[sample][chrName][split_direction].keys():
                            split_reads_array[sample][split_direction][pos - start_win] = \
                                split_reads[sample][chrName][split_direction][pos]

            # gc_array = bw_gc.values('chr' + chrName, start_win, end_win)
            # assert len(gc_array) == win_len
            mappability_array = bw_map.values(chrName, start_win, end_win)
            assert len(mappability_array) == win_len

            # Fill the numpy vstack
            vstack_list = []
            for sample in sample_list:
                # logging.info("Considering sample %s" % sample)

                vstack_list.append(coverage_array[sample])

                for clipped_arrangement in ['left', 'right']:
                    vstack_list.append(clipped_reads_array[sample][clipped_arrangement])

                for mate_position in ['before', 'after']:
                    vstack_list.append(clipped_reads_inversion_array[sample][mate_position])
                for mate_position in ['before', 'after']:
                    vstack_list.append(clipped_reads_duplication_array[sample][mate_position])
                for orientation in ['opposite', 'same']:
                    vstack_list.append(clipped_reads_translocation_array[sample][orientation])

                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['left', 'right', 'unclipped']:
                        # vstack_list.append(
                        #     clipped_read_distance_array[sample][direction][clipped_arrangement])
                        # vstack_list.append(
                        #     clipped_read_distance_num[sample][direction][clipped_arrangement])
                        # vstack_list.append(
                        #     clipped_read_distance_median[sample][direction][clipped_arrangement])
                        vstack_list.append(
                            clipped_read_distance_outlier[sample][direction][clipped_arrangement])
                for direction in ['left', 'right']:
                    vstack_list.append(split_reads_array[sample][direction])
                for direction in ['left', 'right']:
                    vstack_list.append(split_read_distance_array[sample][direction])
                    vstack_list.append(split_read_distance_num[sample][direction])
                    vstack_list.append(split_read_distance_median[sample][direction])

            # vstack_list.append(gc_array)
            vstack_list.append(mappability_array)

            # append one hot encoded sequence for the genomic region
            for nuc in ['A', 'T', 'C', 'G', 'N']:
                one_hot_n = get_one_hot_sequence(chrName, start_win, end_win, nuc, HPC_MODE)
                assert len(one_hot_n) == win_len
                vstack_list.append(one_hot_n)

            # logging.info("Shape of channel matrix: %s" % str(ch_vstack.shape))
            ch_vstack = np.vstack(vstack_list)
            ch_list.append(ch_vstack)

    outDir = os.path.dirname(outFile)
    labelDir = outDir + '/label/'
    create_dir(labelDir)

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    # Write labels for noSV category
    # if trainingMode and sampleName == 'noSV':
    if trainingMode and (sampleName == 'N1' or sampleName == 'N2'):
        label = ['noSV'] * len(ch_list)
        with gzip.GzipFile(labelDir + sampleName + '_' + chrName + '_label.npy.gz', "w") as f:
            np.save(file=f, arr=label)
        f.close()

    logging.info('Number of windows: %d' % len(ch_list))


def main():
    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    # Path on the HPC for the test BAM file
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + 'T0_dedup.bam'

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='2',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output")
    parser.add_argument('-s', '--sample', type=str, default='NA12878',
                        help="Specify sample")
    parser.add_argument('-m', '--svmode', type=str, default='INDEL',
                        help="Specify SV type")
    parser.add_argument('-t', '--train', type=bool, default=True,
                        help="Specify if training mode is active")
    parser.add_argument('-l', '--logfile', default='channel_maker.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()

    # Get chromosome length
    # chrlen = get_chr_len(args.bam, args.chr)

    if args.svmode == 'TRA':
        chrList = args.chr.split('_')
    else:
        chrList = [args.chr]

    channel_maker(ibam=args.bam, chrList=chrList, sampleName=args.sample, SVmode=args.svmode,
                  trainingMode=args.train, outFile=args.out)

    # for sampleName in ['NA12878', 'Patient1', 'Patient2']:
    #     create_labels_nanosv_vcf(sampleName=sampleName, ibam=args.bam)

    # create_labels_bed(sampleName=args.sample, ibam=args.bam)

    # Generate labels for the datasets of real data (HMF or GiaB)
    # generate_labels()
    # read_SURVIVOR_merge_VCF(sampleName=args.sample)

    # print(read_SURVIVOR_merge_VCF(sampleName=args.sample))

    # write_sv_without_cr(sampleName=args.sample, ibam=args.bam)

    # clipped_read_positions_to_bed(sampleName=args.sample)
    # nanosv_vcf_to_bed(sampleName=args.sample)

    # get_nanosv_manta_sv_from_SURVIVOR_merge_VCF(sampleName=args.sample)

    # print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))
    print('Elapsed time channel_maker_real = %f' % (time() - t0))


if __name__ == '__main__':
    main()
