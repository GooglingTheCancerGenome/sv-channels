#!/usr/bin/env python3

"""
Split VCF into SVs that overlaps SR positions and SVs that do not
"""

import argparse
import gzip
import os
import sys
from collections import defaultdict
from time import time

sys.path.append('../genome_wide/')

from functions import *
from intervaltree import IntervalTree
from label_windows import *


def parse_cl_args(in_args, caller):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    if caller == 'test':
        mypath = os.path.join(pathout, 'test.bedpe')
    else:
        mypath = os.path.join(pathout, caller + '.bedpe')
    description = "Split VCF by overlap on split read positions"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input", type=str,
                        default=mypath,
                        help="VCF input file")
    parser.add_argument("-o", "--output", type=str,
                        default=os.path.join(pathout, caller+'_sr.bedpe'),
                        help="VCF output file")
    parser.add_argument("-n", "--output_nosr", type=str,
                        default=os.path.join(pathout, caller+'_nosr.bedpe'),
                        help="VCF output file")
    parser.add_argument("-w", "--win", type=int,
                        default=200,
                        help="VCF output file")
    parser.add_argument("-s", "--split_reads", type=str,
                        default=os.path.join(pathout, 'split_reads.bedpe.gz'),
                        help="BEDPE file with split read positions")
    parser.add_argument('-sv',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    args = parser.parse_args(in_args)
    return args


def get_srpos_from_bedpe(inbedpe, svtype):
    # Dictionary with chromosome keys to store SVs
    srpos = set()
    with gzip.GzipFile(inbedpe, 'rb') as fin:
        for line in fin.readlines():
            columns = line.decode('utf8').rstrip().split("\t")
            # print(columns)
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])
            type = columns[6]

            if svtype == type:
                srpos.add((chrom1, pos1_start, chrom2, pos2_start))
    print('{} split read positions for {}'.format(len(srpos), svtype))
    return srpos


def create_gtrees(sv_list):
    print('Building SV GenomicTrees...')
    # Tree with windows for candidate positions
    trees_start = defaultdict(IntervalTree)
    trees_end = defaultdict(IntervalTree)
    # Populate tree
    for i, sv in enumerate(sv_list, start=0):
        chrom1, pos1_start, pos1_end, chrom2, pos2_start, pos2_end, svtype = sv
        sv_id = '_'.join(
            (chrom1, str(pos1_start), chrom2, str(pos2_start)))
        trees_start[chrom1][pos1_start:pos1_end] = (i, sv_id)
        trees_end[chrom2][pos2_start:pos2_end] = (i, sv_id)
    return trees_start, trees_end


def search_tree_with_bedpe(cpos, trees_start, trees_end, win_hlen):
    print('Searching SV GenomicTrees with candidate positions...')
    lookup_start = []
    lookup_end = []
    # Log info every n_r times
    n_r = 10**4
    last_t = time()
    for i, p in enumerate(cpos, start=1):
        chrom1, pos1, chrom2, pos2 = p
        lookup_start.append(trees_start[chrom1].envelop(
                            pos1 - win_hlen, pos1 + win_hlen + 1))
        lookup_end.append(
            trees_end[chrom2].envelop(
                            pos2 - win_hlen, pos2 + win_hlen + 1))
    return lookup_start, lookup_end


def main():

    for caller in ['manta', 'gridss', 'lumpy', 'delly', 'test']:

        print('Considering {}:'.format(caller))
        args = parse_cl_args(sys.argv[1:], caller)

        srpos = get_srpos_from_bedpe(os.path.join(pathout, 'split_reads.bedpe.gz'), args.svtype)

        win_hlen = int(
            args.win / 2) if args.win % 2 == 0 else int((args.win + 1) / 2)
        filename, file_extension = os.path.splitext(args.input)
        if file_extension == '.vcf':
            sv_list = read_vcf(args.input)
            sv_list = [i for i in sv_list if i[-1] == args.svtype]
        elif file_extension == '.bedpe':
            sv_list = read_bedpe(args.input, args.svtype)

        trees_start, trees_end = create_gtrees(sv_list)
        lookup_start, lookup_end = search_tree_with_bedpe(
            srpos, trees_start, trees_end, win_hlen)
        idx = []

        sr_sv = set()
        sr_nosv = set()

        for i, x in enumerate(zip(srpos, lookup_start, lookup_end), start=0):
            p, lu_start, lu_end = x
            l1 = len(lu_start)
            l2 = len(lu_end)
            if l1 > 0 and l2 > 0:

                sr_sv.add(p)
                lu_start_set = []
                lu_end_set = []
                for s in lu_start:
                    lu_start_elem_start, lu_start_elem_end, lu_start_elem_data = s
                    lu_start_i, lu_start_elem_svid = lu_start_elem_data
                    lu_start_set.append(lu_start_i)

                for s in lu_end:
                    lu_end_elem_start, lu_end_elem_end, lu_end_elem_data = s
                    lu_end_i, lu_end_elem_svid = lu_end_elem_data
                    lu_end_set.append(lu_end_i)

                olap = set(lu_start_set) & set(lu_end_set)
                if len(olap) > 0:
                    idx.extend(olap)
            elif l1 == 0 and l2 == 0:
                sr_nosv.add(p)

        print('list of indices: {}'.format(len(idx)))
        out = open(args.output, 'w')
        out_nosr = open(args.output_nosr, 'w')
        j = 0
        sr_lines = nosr_lines = 0
        with open(args.input, 'r') as fin:
            for line in fin.readlines():
                if line[0] == '#':
                    out.write(line)
                    out_nosr.write(line)
                else:
                    if j in set(idx):
                        out.write(line)
                        sr_lines += 1
                    else:
                        out_nosr.write(line)
                        nosr_lines += 1
                    j += 1

        print('{} SVs => with {} SR:{}({}%); without {} SR:{}({}%)'.format(
            sr_lines+nosr_lines,
            len(sr_sv),
            sr_lines,
            sr_lines/(sr_lines + nosr_lines) * 100,
            len(sr_nosv),
            nosr_lines,
            nosr_lines / (sr_lines + nosr_lines) * 100,
        ))


if __name__ == "__main__":
    global pathout
    pathout = '/Users/lsantuari/Documents/Processed/SPIDER/split_reads/CHM1_CHM13/bedpe'
    main()
