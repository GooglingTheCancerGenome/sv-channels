# Imports

import re
import pysam
import statistics
from pysam import VariantFile
from collections import Counter
from intervaltree import Interval, IntervalTree
# from collections import defaultdict

import numpy as np
from scipy.spatial import distance

import argparse
import bz2file
import gzip
import os, errno
import pickle
from time import time
import logging
import csv

import pandas as pd

# from pybedtools import BedTool
# from ggplot import *

# Flag used to set either paths on the local machine or on the HPC
HPC_MODE = False


class SVRecord:

    def __init__(self, record):
        if type(record) != pysam.VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')

        self.chrom = record.chrom
        self.chrom2 = record.info['CHR2']
        self.start = record.pos
        self.end = record.stop
        self.supp_vec = record.info['SUPP_VEC']
        self.svtype = record.info['SVTYPE']


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
        # extract strand/orientation, position, and (possibly) inserted string
        # length from record eg [9:1678896[N.
        # result is result from re.match with the explicit BP regexp

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

        setupREs()

        altstr = str(record.alts[0])
        resultBP = re.match(__bpRE__, altstr)

        if resultBP:
            assert resultBP.group(2) == resultBP.group(4)
            ct, chr2, pos2, indellen = self.locFromBkpt(str(record.ref), resultBP.group(1),
                                                        resultBP.group(2), resultBP.group(3), resultBP.group(4),
                                                        resultBP.group(5))
        return (ct, chr2, pos2, indellen)


__bpRE__ = None
__symbolicRE__ = None


def setupREs():
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

    min_CR_support = 3

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

    cr_pos = [elem for elem, cnt in cpos.items() if cnt >= min_CR_support]

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
    to the
    :return:
    '''

    # Only clipped read positions supported by at least min_cr_support clipped reads are considered
    min_cr_support = 3
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

            vcf_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth/'+sampleName+'/SV'
            vcf_files = dict()

            for mapper in ['last']:
                vcf_files[mapper] = dict()

                vcf_files[mapper]['nanosv'] = vcf_dir + '/' + mapper + '/' + mapper + '_nanosv.sorted.vcf'
                assert os.path.isfile(vcf_files[mapper]['nanosv'])

    return vcf_files


def read_nanosv_vcf(sampleName):

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

        sv = [svrec for svrec in sv if svrec.svtype =='DEL'
              if 'LowQual' not in list(svrec.filter)]

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
        #print(s)

        return sv


def create_labels_nanosv_vcf(sampleName, chrName, ibam):

    def closest_loc(pos, pos_list):
        pos_array = np.asarray(pos_list)
        deltas = np.abs(pos_array - pos)
        idx = np.argmin(deltas)
        return (pos_list[idx], deltas[idx])

    # Load SV list
    sv_list = read_nanosv_vcf(sampleName)
    sv_list = [sv for sv in sv_list if sv.svtype == 'DEL']
    chr_list = set([var.chrom for var in sv_list])
    print(chr_list)
    print('# DELs:%d' % len(sv_list))

    cnt = Counter([sv.chrom for sv in sv_list])
    s = pd.Series([v for v in cnt.values()], index=cnt.keys())
    print(s)

    assert sum(s) == len(sv_list)

    win_hlen = 25
    win_len = win_hlen * 2

    confint = 100

    for chrName in chr_list:

        chrLen = get_chr_len(ibam, chrName)

        # if chrName == '17':

        # Load CR positions
        cr_pos = load_clipped_read_positions(sampleName, chrName)

        # Remove positions with windows falling off chromosome boundaries
        cr_pos = [pos for pos in cr_pos if win_hlen <= pos <= (chrLen - win_hlen)]

        print(sorted(cr_pos))

        # Using IntervalTree for interval search
        t = IntervalTree()

        for var in sv_list:
            if var.chrom == chrName:
                # cipos[0] and ciend[0] are negative in the VCF file
                id_start =  var.svtype + '_start'
                id_end = var.svtype + '_end'

                #id_start = '_'.join((var.chrom, str(var.start+var.cipos[0]),  str(var.start+var.cipos[1])))
                #id_end = '_'.join((var.chrom, str(var.end + var.ciend[0]), str(var.end+var.ciend[1])))

                t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
                t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

                #t[var.start - confint:var.start + confint + 1] = var.svtype + '_start'
                #t[var.end - confint:var.end + confint + 1] = var.svtype + '_end'

        label = [sorted(t[p]) for p in cr_pos]

        label_ci = []
        for elem in label:
            if len(elem) == 1:
                #print(elem)
                label_ci.append(elem[0].data)
            else:
                label_ci.append('noSV')

        print('CR positions: %d' % len(cr_pos))
        print('Label length: %d' % len(label))
        assert len(label_ci) == len(cr_pos)
        print(Counter(label_ci))

        hit_set = set([elem for l in label for elem in l])

        no_hits = [ (var.start + var.cipos[0], var.start + var.cipos[1]) for var in sv_list if var.chrom == chrName
                   if var.start + var.cipos[0] not in [start for start,end,lab in hit_set]]
        no_hits = sorted(no_hits)

        #print(len(no_hits))
        #print(no_hits)

        label_close = []
        label_BPJ = []
        label_distance = []

        start_positions = [var.start for var in sv_list if var.chrom == chrName]
        end_positions = [var.end for var in sv_list if var.chrom == chrName]

        # Very slow:
        if len(start_positions) != 0 and len(end_positions) != 0:
            for pos in cr_pos:
                #print('Pos: %d' % pos)
                close_start, close_start_dist = closest_loc(pos, start_positions)
                close_end, close_end_dist = closest_loc(pos, end_positions)
                if close_start_dist <= confint or close_end_dist <= confint:
                    if close_start_dist <= close_end_dist:
                        #print(str(close_start_dist))
                        bpj_id = str(chrName) + '_' + str(close_start)
                        #print(bpj_id)
                        label_close.append('DEL_start')
                        label_BPJ.append(bpj_id)
                        label_distance.append(abs(pos - close_start))
                    else:
                        #print(str(close_end_dist))
                        bpj_id = str(chrName) + '_' + str(close_end)
                        #print(bpj_id)
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


# END: NanoSV specific functions

def read_SURVIVOR_merge_VCF(sampleName):
    if HPC_MODE:
        channel_dir = ''
        vcf_file = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Breast_Cancer_Pilot_SV/Manta/survivor_merge.vcf'
    else:
        channel_dir = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/ChannelMaker'
        vcf_file = '/Users/lsantuari/Documents/Data/Breast_Cancer_Pilot/SV/survivor_merge.vcf'

    vec = 'clipped_read_pos'
    min_CR_support = 3
    confint = 100

    vcf_in = VariantFile(vcf_file)
    samples_list = list((vcf_in.header.samples))
    samples = [w.split('_')[0].split('/')[1] for w in samples_list]

    # sv_bedtool = BedTool(vcf_file)
    # print(sv_bedtool)

    sv = []

    # create sv list with SVRecord objects
    for rec in vcf_in.fetch():
        # print(rec)
        # print(rec.info['SUPP_VEC'][16])
        # avoid SVs on chromosomes Y and MT
        if rec.chrom not in ['Y', 'MT'] and rec.info['CHR2'] not in ['Y', 'MT']:
            # print(rec)
            sv.append(SVRecord(rec))

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

        cr_pos = [elem for elem, cnt in cpos['Tumor'][chr].items() if cnt >= min_CR_support]
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
    :return: None
    '''

    no_cr_File = '/Users/lsantuari/Documents/Data/HMF/ChannelMaker_results/HMF_Tumor/labels/' \
                 + 'no_clipped_read_positions.pk.gz'
    with gzip.GzipFile(no_cr_File, "r") as f:
        no_clipped_read_pos = pickle.load(f)
    f.close()
    print(list(no_clipped_read_pos))


def channel_maker(ibam, chrName, sampleName, trainingMode, outFile):
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
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Use a 50 bp window length
    # TODO: the window length should be an input parameter
    win_hlen = 25
    win_len = win_hlen * 2
    # Minimum clipped read support to consider
    min_cr_support = 3

    # Set the correct prefix for the path
    if trainingMode and ( sampleName == 'N1' or sampleName == 'N2' ):
        prefix_train = 'Training/'
    else:
        prefix_train = ''

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
        sample_list = [sampleName]

        # File with clipped read positions, output of the clipped_read_pos script
        clipped_read_pos_file = prefix_train + sampleName + \
                                '/clipped_read_pos/' + chrName + '_clipped_read_pos.pbz2'
        # File with the clipped read distances, output of the clipped_read_distance script
        clipped_read_distance_file = 'clipped_read_distance/' + chrName + '_clipped_read_distance.pbz2'
        # File with the clipped reads, output of the clipped_reads script
        clipped_reads_file = 'clipped_reads/' + chrName + '_clipped_reads.pbz2'
        # File with the coverage array, output of the coverage script
        coverage_file = 'coverage/' + chrName + '_coverage.npy.bz2'
        # File with the split reads and split read distance, output of the split_read_distance script
        split_read_distance_file = 'split_read_distance/' + chrName + '_split_read_distance.pbz2'

        # Check file existence
        print('Checking file: %s' % clipped_read_pos_file)
        assert os.path.isfile(clipped_read_pos_file)

        for sample in sample_list:
            assert os.path.isfile(prefix_train + sample + '/' + clipped_read_distance_file)
            assert os.path.isfile(prefix_train + sample + '/' + clipped_reads_file)
            assert os.path.isfile(prefix_train + sample + '/' + coverage_file)
            assert os.path.isfile(prefix_train + sample + '/' + split_read_distance_file)

    logging.info('Chromosome %s' % chrName)

    logging.info('Reading clipped read positions')
    with bz2file.BZ2File(clipped_read_pos_file, 'rb') as f:
        clipped_pos_cnt = pickle.load(f)
    logging.info('End of reading')

    # Count the number of clipped read positions with a certain minimum number of clipped reads
    count_clipped_read_positions(clipped_pos_cnt)

    # Load channel data
    # Dictionaries where to load the channel data
    clipped_read_distance = dict()
    clipped_reads = dict()
    clipped_reads_inversion = dict()
    clipped_reads_duplication = dict()
    coverage = dict()
    split_reads = dict()
    split_read_distance = dict()

    for sample in sample_list:
        prefix = prefix_train + sample + '/' if HPC_MODE else ''

        logging.info('Considering %s' % sample)
        logging.info('Reading clipped read distances')
        with bz2file.BZ2File(prefix + clipped_read_distance_file, 'rb') as f:
            clipped_read_distance[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading clipped reads')
        with bz2file.BZ2File(prefix + clipped_reads_file, 'rb') as f:
            clipped_reads[sample], clipped_reads_inversion[sample], clipped_reads_duplication[sample] = pickle.load(f)
        logging.info('End of reading')

        logging.info('Reading coverage')
        with bz2file.BZ2File(prefix + coverage_file, 'rb') as f:
            coverage[sample] = np.load(file=f)
        logging.info('End of reading, coverage length: %d out of %d' % (len(coverage[sample]), chrLen))

        logging.info('Reading split read distances')
        with bz2file.BZ2File(prefix + split_read_distance_file, 'rb') as f:
            split_read_distance[sample], split_reads[sample] = pickle.load(f)
        logging.info('End of reading')

    # If in 'NoSV' mode, consider all the clipped read positions (minimum clipped read support equal to 0)
    if trainingMode and ( sampleName == 'N1' or sampleName == 'N2' ):
        clipped_pos = [k for k, v in clipped_pos_cnt.items()]
    else:
        clipped_pos = [k for k, v in clipped_pos_cnt.items() if v >= min_cr_support]

    # print(clipped_pos)

    # Dictionaries where to store the channel arrays as generated from the dictionaries
    clipped_read_distance_array = dict()
    clipped_read_distance_num = dict()
    clipped_read_distance_median = dict()

    clipped_reads_array = dict()
    clipped_reads_inversion_array = dict()
    clipped_reads_duplication_array = dict()

    coverage_array = dict()

    split_read_distance_array = dict()
    split_read_distance_num = dict()
    split_read_distance_median = dict()

    split_reads_array = dict()

    # Log info every n_r times
    n_r = 10 ** 3
    # print(n_r)
    last_t = time()

    for i, center_pos in enumerate(clipped_pos, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d clipped positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Check if center_pos is within the chromosome boundaries
        if win_hlen <= center_pos <= (chrLen - win_hlen):

            start_win = center_pos - win_hlen
            end_win = center_pos + win_hlen

            # Build arrays for the numpy vstack
            for sample in sample_list:

                # clipped read distance
                clipped_read_distance_array[sample] = dict()
                clipped_read_distance_num[sample] = dict()
                clipped_read_distance_median[sample] = dict()

                for direction in ['forward', 'reverse']:
                    clipped_read_distance_array[sample][direction] = dict()
                    clipped_read_distance_num[sample][direction] = dict()
                    clipped_read_distance_median[sample][direction] = dict()

                for direction in ['forward', 'reverse']:
                    # for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
                    for clipped_arrangement in ['left', 'right']:
                        clipped_read_distance_array[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                       dtype=int)
                        clipped_read_distance_num[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                     dtype=int)
                        clipped_read_distance_median[sample][direction][clipped_arrangement] = np.zeros(win_len,
                                                                                                        dtype=int)
                        for pos in range(start_win, end_win):
                            if pos in clipped_read_distance[sample][direction][clipped_arrangement].keys():
                                clipped_read_distance_array[sample][direction][clipped_arrangement][pos - start_win] = \
                                    sum(clipped_read_distance[sample][direction][clipped_arrangement][pos])
                                clipped_read_distance_num[sample][direction][clipped_arrangement][pos - start_win] = \
                                    len(clipped_read_distance[sample][direction][clipped_arrangement][pos])
                                clipped_read_distance_median[sample][direction][clipped_arrangement][pos - start_win] = \
                                    statistics.median(
                                        clipped_read_distance[sample][direction][clipped_arrangement][pos])

                        # print(clipped_read_distance_array[direction][clipped_arrangement])

                # clipped reads
                clipped_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    clipped_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads[sample][split_direction].keys():
                            clipped_reads_array[sample][split_direction][pos - start_win] = \
                                clipped_reads[sample][split_direction][pos]

                # clipped reads inversions
                clipped_reads_inversion_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_inversion_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_inversion[sample][mate_position].keys():
                            clipped_reads_inversion_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_inversion[sample][mate_position][pos]

                # clipped reads duplication
                clipped_reads_duplication_array[sample] = dict()
                for mate_position in ['before', 'after']:
                    clipped_reads_duplication_array[sample][mate_position] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in clipped_reads_duplication[sample][mate_position].keys():
                            clipped_reads_duplication_array[sample][mate_position][pos - start_win] = \
                                clipped_reads_duplication[sample][mate_position][pos]

                # coverage
                coverage_array[sample] = coverage[sample][start_win:end_win]
                assert len(coverage_array[sample]) == win_len

                # split read distance
                split_read_distance_array[sample] = dict()
                split_read_distance_num[sample] = dict()
                split_read_distance_median[sample] = dict()

                for split_direction in ['left', 'right']:
                    split_read_distance_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    split_read_distance_num[sample][split_direction] = np.zeros(win_len, dtype=int)
                    split_read_distance_median[sample][split_direction] = np.zeros(win_len, dtype=int)

                    if pos in split_read_distance[sample][split_direction].keys():
                        split_read_distance_array[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][split_direction][pos])
                        split_read_distance_num[sample][split_direction][pos - start_win] = \
                            len(split_read_distance[sample][split_direction][pos])
                        split_read_distance_median[sample][split_direction][pos - start_win] = \
                            sum(split_read_distance[sample][split_direction][pos])

                # split reads
                split_reads_array[sample] = dict()
                for split_direction in ['left', 'right']:
                    split_reads_array[sample][split_direction] = np.zeros(win_len, dtype=int)
                    for pos in range(start_win, end_win):
                        if pos in split_reads[sample][split_direction].keys():
                            split_reads_array[sample][split_direction][pos - start_win] = \
                                split_reads[sample][split_direction][pos]

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

                for direction in ['forward', 'reverse']:
                    for clipped_arrangement in ['left', 'right']:
                        vstack_list.append(
                            clipped_read_distance_array[sample][direction][clipped_arrangement])
                        vstack_list.append(
                            clipped_read_distance_num[sample][direction][clipped_arrangement])
                        vstack_list.append(
                            clipped_read_distance_median[sample][direction][clipped_arrangement])
                for direction in ['left', 'right']:
                    vstack_list.append(split_reads_array[sample][direction])
                for direction in ['left', 'right']:
                    vstack_list.append(split_read_distance_array[sample][direction])
                    vstack_list.append(split_read_distance_num[sample][direction])
                    vstack_list.append(split_read_distance_median[sample][direction])

            # logging.info("Shape of channel matrix: %s" % str(ch_vstack.shape))
            ch_vstack = np.vstack(vstack_list)
            ch_list.append(ch_vstack)

    outDir = os.path.dirname(outFile)
    labelDir = outDir +'/label/'
    create_dir(labelDir)

    # Save the list of channel vstacks
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=ch_list)
    f.close()

    # Write labels for noSV category
    #if trainingMode and sampleName == 'noSV':
    if trainingMode and ( sampleName == 'N1' or sampleName == 'N2' ):
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
    parser.add_argument('-s', '--sample', type=str, default='Patient2',
                        help="Specify sample")
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
    #chrlen = get_chr_len(args.bam, args.chr)

    #channel_maker(ibam=args.bam, chrName=args.chr, sampleName=args.sample,
    #               trainingMode=args.train, outFile=args.out)

    create_labels_nanosv_vcf(sampleName=args.sample, chrName=args.chr, ibam=args.bam)

    # Generate labels for the datasets of real data (HMF or GiaB)
    # generate_labels()
    # read_SURVIVOR_merge_VCF(sampleName=args.sample)

    print('Elapsed time channel_maker_real on BAM %s and Chr %s = %f' % (args.bam, args.chr, time() - t0))


if __name__ == '__main__':
    main()