import os
import argparse
import glob

from pysam import VariantFile

import sys

# add scripts folder in path
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.realpath(__file__)
            )
        )
    )
)
from genome_wide.label_classes import SVRecord

# Threadpool configuration
# avoid bcolz error: nthreads cannot be larger than environment variable
os.environ["NUMEXPR_MAX_THREADS"] = "128"

import bcolz
import numpy as np
import tensorflow as tf
from random import sample
from intervaltree import IntervalTree
from collections import defaultdict, Counter

import gzip
import json

# set GPU options
# allow_growth allows fair share of GPU memory across processes
gpu_options = tf.GPUOptions(allow_growth=True)
sess = tf.Session(config=tf.ConfigProto(gpu_options=gpu_options))
tf.keras.backend.set_session(sess)


def get_chr_array(channels_dir, chrom):
    carray_fname = glob.glob(os.path.join(channels_dir, 'chr_array',
                                          '*_' + chrom + '_carray'))
    assert len(carray_fname) == 1, 'Not a single carray folder found:' + carray_fname
    carray_fname = carray_fname[0]
    assert os.path.exists(carray_fname), carray_fname + ' not found'

    return bcolz.open(rootdir=carray_fname)


def get_nchannels(genome_array):
    dict_item = next(iter(genome_array.values()))
    return dict_item.shape[1]


def get_truth_set_trees(truth_set_file):

    def read_vcf(invcf):

        # Check file existence
        assert os.path.isfile(invcf), invcf + ' not found!'
        # Dictionary with chromosome keys to store SVs
        sv_list = []

        vcf_in = VariantFile(invcf, 'r')

        for rec in vcf_in.fetch():

            var = SVRecord(rec, '')

            chrom1 = var.chrom
            pos1_start = var.start + var.cipos[0]
            pos1_bp = var.start
            pos1_end = var.start + var.cipos[1] + 1

            chrom2 = var.chrom2
            pos2_start = var.end + var.ciend[0]
            pos2_bp = var.end
            pos2_end = var.end + var.ciend[1] + 1
            svtype = var.svtype

            # choose only deletions?
            if svtype == "DEL":
                sv_list.append((
                    chrom1, pos1_start, pos1_bp, pos1_end,
                    chrom2, pos2_start, pos2_bp, pos2_end
                ))

        print('{} SVs'.format(len(sv_list)))

        return sv_list

    file_basename = os.path.basename(truth_set_file)

    filename, file_ext = os.path.splitext(file_basename)

    bp1_tree = None
    bp2_tree = None
    sv_list = None

    if file_basename == 'Personalis_1000_Genomes_deduplicated_deletions.bed':

        # read the svclassify truth set and return two dictionaries bp1_tree and bp2_tree,
        # keys are chromosome names, values are IntervalTree objects
        # downloaded from:'ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/technical/svclassify_Manuscript/'+
        # 'Supplementary_Information/Personalis_1000_Genomes_deduplicated_deletions.bed'

        print('Reading truth set: {}'.format(file_basename))

        with open(truth_set_file, "r") as f:
            file_content = f.read()
            # split lines
            file_content = file_content.splitlines()
            # remove header. It should not be present in a BED file
            file_content = file_content[1:]
        f.close()

        # split by tab
        sv_list = [l.split('\t') for l in file_content]
        print('{} deletions in the truth set'.format(len(sv_list)))

        bp1_tree = defaultdict(IntervalTree)
        bp2_tree = defaultdict(IntervalTree)

        for c, s, e in sv_list:
            interval_id = '_'.join([c, s, e])
            s = int(s)
            e = int(e)
            bp1_tree[c][s:s + 1] = interval_id
            bp2_tree[c][e:e + 1] = interval_id

    elif file_basename == 'lumpy-Mills2011-call-set.nanosv.sorted.bed':

        # read the lumpy-nanosv truth set and return two dictionaries bp1_tree and bp2_tree,
        # keys are chromosome names, values are IntervalTree objects

        with open(truth_set_file, "r") as f:
            file_content = f.read()
            # split lines
            file_content = file_content.splitlines()
        f.close()

        # split by tab
        sv_list = [l.split('\t') for l in file_content]

        print('{} breakpoint intervals in truth set'.format(len(sv_list)))

        bp1_tree = defaultdict(IntervalTree)
        bp2_tree = defaultdict(IntervalTree)

        chr_list = []
        for l in sv_list:
            c = l[0]
            chr_list.append(c)
            s = l[1]
            e = l[2]
            interval_id = '_'.join([c, s, e])
            s = int(s)
            e = int(e)
            if l[3] == 'DEL_start':
                bp1_tree[c][s:e + 1] = interval_id
            elif l[3] == 'DEL_end':
                bp2_tree[c][s:e + 1] = interval_id

    elif file_ext in ['.vcf', '.gz']:

        sv_list = read_vcf(truth_set_file)

        print('{} SVs in truth set'.format(len(sv_list)))

        bp1_tree = defaultdict(IntervalTree)
        bp2_tree = defaultdict(IntervalTree)

        for c1, s1, p1, e1, c2, s2, p2, e2 in sv_list:
            interval_id = '_'.join([c1, str(s1), str(e1), c2, str(s2), str(e2)])
            bp1_tree[c1][s1:e1 + 1] = interval_id
            bp2_tree[c2][s2:e2 + 1] = interval_id

    # sanity check
    if bp1_tree is not None and bp2_tree is not None:
        for i, j in zip(bp1_tree.items(), bp2_tree.items()):
            k1, v1 = i
            k2, v2 = j
            # print('{} : {}'.format(len(v1), len(v2)))
            assert len(v1) == len(v2)

    return sv_list, bp1_tree, bp2_tree


def generate_positive_set(win, sv_list, genome_array):

    win_hlen = int(win / 2)

    a_size = len(sv_list)

    n_channels = get_nchannels(genome_array)

    a = np.zeros((a_size, win_hlen * 2, n_channels))
    b_s = bcolz.carray(a)
    b_e = bcolz.carray(a)

    for i, sv in enumerate(sv_list):
        if i % 100 == 0 and i != 0:
            print(i)

        if len(sv_list[0]) == 3:
            # from (chrom, bp1_pos, bp2_pos)
            c = sv[0]
            s = int(sv[1])
            e = int(sv[2])
            b_s[i, :, :] = genome_array[c][s - win_hlen:s + win_hlen, :]
            b_e[i, :, :] = genome_array[c][e - win_hlen:e + win_hlen, :]

        elif len(sv_list[0]) == 3 * 2:
            # generate positive set from GRIDSS
            # from (chrom1, bp1_start, bp1_end, chrom2, bp2_start, bp2_end)
            c1 = sv[0]
            s1 = int(sv[2])
            c2 = sv[4]
            s2 = int(sv[6])
            b_s[i, :, :] = genome_array[c1][s1 - win_hlen:s1 + win_hlen, :]
            b_e[i, :, :] = genome_array[c2][s2 - win_hlen:s2 + win_hlen, :]

    print(b_s.shape)
    print(b_e.shape)

    return b_s, b_e


def generate_negative_set(win, sv_list, bp1_tree, bp2_tree, genome_array):
    # 1) generate random (chr, pos)
    # 2) create negative set with shape (n_svs, win_hlen*2, n_channels) for each random position that does not overlap
    # the truth set

    def get_random_positions(bp1_tree, bp2_tree):

        # random chromosomes?
        # chrom_choice = np.random.choice(list(genome_array.keys()), n_svs)
        # keep same number of sv per chromosome?
        chrom_choice = np.array([l[0] for l in sv_list])
        pos_choice = np.zeros(len(sv_list))

        bp_tree = {c: bp1_tree[c] | bp2_tree[c] for c in bp1_tree.keys()}

        for c in np.unique(chrom_choice):
            if c in bp_tree.keys():
                chr_len = genome_array[c].shape[0]
                idx_matches = np.where(chrom_choice == c)[0]
                while True:
                    # print('Chr{}, idx_matches: {}'.format(c, idx_matches))
                    pos_choice[idx_matches] = np.random.randint(0, chr_len, len(idx_matches))

                    matches = [bp_tree[c][p - win_hlen:p + win_hlen] for p in list(pos_choice[idx_matches])]

                    matches_j = [i for i, l in enumerate(matches) if len(l) > 0]

                    if len(matches_j) == 0:
                        break
                    idx_matches = idx_matches[matches_j]
        # print('pos_choice: {}'.format(pos_choice[idx_matches]))
        for c in np.unique(chrom_choice):
            idx_matches = np.where(chrom_choice == c)[0]
            pos_choice[idx_matches] = np.sort(pos_choice[idx_matches])
        res = [[c, int(p)] for c, p in zip(chrom_choice, pos_choice)]

        return res

    win_hlen = int(win / 2)
    n_channels = get_nchannels(genome_array)

    rand_pos = get_random_positions(bp1_tree, bp2_tree)

    # generate negative set

    a_size = len(rand_pos)

    a = np.zeros((a_size, win_hlen * 2, n_channels))
    b = bcolz.carray(a)

    for i, sv in enumerate(rand_pos):
        if i % 100 == 0 and i != 0:
            print(i)
        c = sv[0]
        p = sv[1]

        # print('{} {} {}'.format(i, c, p))
        b[i, :, :] = genome_array[c][p - win_hlen:p + win_hlen, :]

    print(b.shape)

    return b


def positive(args):

    def save_positive_windows(outfile, start_array, end_array):
        np.savez_compressed(outfile,
                            start=np.array(start_array),
                            end=np.array(end_array))
        size_in_bytes = os.path.getsize(outfile)
        print('{} MB'.format(size_in_bytes / 10 ** 6))

    genome_array = {chrom: get_chr_array(args.inputdir, chrom) for chrom in args.chrlist}
    # print(genome_array)

    sv_list, bp1_tree, bp2_tree = get_truth_set_trees(args.truthset)
    # only keep SVs on chromosomes from chrlist
    sv_list = [e for e in sv_list if e[0] in args.chrlist]

    b_s, b_e = generate_positive_set(args.win, sv_list, genome_array)

    save_positive_windows(args.output, b_s, b_e)


def negative(args):

    def save_negative_windows(outfile, array):
        np.savez_compressed(outfile,
                            neg=np.array(array)
                            )
        size_in_bytes = os.path.getsize(outfile)
        print('{} MB'.format(size_in_bytes / 10 ** 6))

    genome_array = {chrom: get_chr_array(args.inputdir, chrom) for chrom in args.chrlist}
    # print(genome_array)

    sv_list, bp1_tree, bp2_tree = get_truth_set_trees(args.truthset)
    # only keep SVs on chromosomes from chrlist
    sv_list = [e for e in sv_list if e[0] in args.chrlist]

    b = generate_negative_set(args.win, sv_list, bp1_tree, bp2_tree, genome_array)

    save_negative_windows(args.output, b)


def main():

    truth_sets = {
        'svclassify': os.path.join('/Users/lsantuari/Documents/Local_GitHub',
                                   'sv-callers/scripts/data/benchmark/in',
                                   'Personalis_1000_Genomes_deduplicated_deletions.bed'),
        'mills_nanosv': os.path.join('/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878',
                                     'lumpy-Mills2011-call-set.nanosv.sorted.bed'),
        'gridss': '/Users/lsantuari/Documents/Data/germline/NA12878/SV/Filtered/gridss.vcf',
        'manta': '/Users/lsantuari/Documents/Data/germline/NA12878/SV/Filtered/manta.vcf'
    }

    parser = argparse.ArgumentParser(
        description='Generate training data',
        usage='''T0_S1_generate_training_data.py <command> [<args>]

    Commands are:
       positive   Generate positive set
       negative   Generate negative set
    ''')

    subparsers = parser.add_subparsers(help='sub-command help')

    # create the parser for the "positive" command
    parser_pos = subparsers.add_parser('positive', help='positive help')

    parser_pos.add_argument('-chrlist', nargs='+', default=['17'],
                            help="List of chromosomes to consider")
    parser_pos.add_argument('-win', type=int, default=200,
                            help="Window length")
    parser_pos.add_argument('-truthset', nargs='+',
                            default=truth_sets['gridss'],
                            help="Truth set file")
    parser_pos.add_argument('-inputdir', type=str,
                            default='/Users/lsantuari/Documents/Processed/channel_maker_output/T1',
                            help="Specify input data directory")
    parser_pos.add_argument('-output', type=str,
                            default='positive.npz',
                            help="Specify output")
    parser_pos.set_defaults(func=positive)

    # create the parser for the "negative" command
    parser_neg = subparsers.add_parser('negative', help='negative help')

    parser_neg.add_argument('-chrlist', nargs='+', default=['17'],
                            help="List of chromosomes to consider")
    parser_neg.add_argument('-win', type=int, default=200,
                            help="Window length")
    parser_neg.add_argument('-truthset', nargs='+',
                            default=truth_sets['gridss'],
                            help="Truth set file")
    parser_neg.add_argument('-inputdir', type=str,
                            default='/Users/lsantuari/Documents/Processed/channel_maker_output/T1',
                            help="Specify input data directory")
    parser_neg.add_argument('-output', type=str,
                            default='negative.npz',
                            help="Specify output")
    parser_neg.set_defaults(func=negative)

    # args = parser.parse_args()
    # args.func(args)

    # test
    args = parser.parse_args('positive'.split())
    args.func(args)
    args = parser.parse_args('negative'.split())
    args.func(args)


if __name__ == '__main__':
    main()
