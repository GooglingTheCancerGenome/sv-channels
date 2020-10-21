#!/usr/bin/env python3

"""
Split VCF into SVs that overlaps SR positions and SVs that do not
"""

import argparse
import os
import sys
import gzip

sys.path.append('../genome_wide/')
from functions import *
from label_windows import *


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """

    path = '../../data/vcf/gridss_out'

    description = "Split VCF by overlap on split read positions"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input", type=str,
                        default=os.path.join(path, 'gridss.vcf'),
                        help="VCF input file")
    parser.add_argument("-o", "--output", type=str,
                        default=os.path.join(path, 'gridss_sr.vcf'),
                        help="VCF output file")
    parser.add_argument("-n", "--output_nosr", type=str,
                        default=os.path.join(path, 'gridss_nosr.vcf'),
                        help="VCF output file")
    parser.add_argument("-w", "--win", type=int,
                        default=200,
                        help="VCF output file")
    parser.add_argument("-s", "--split_reads", type=str,
                        default='../genome_wide/split_reads/split_reads.bedpe.gz',
                        help="BEDPE file with split read positions")
    args = parser.parse_args(in_args)
    return args


def get_srpos_from_bedpe(inbedpe):

    # Check file existence
    assert os.path.isfile(inbedpe), inbedpe + ' not found!'
    # Dictionary with chromosome keys to store SVs
    srpos = []

    with gzip.GzipFile(inbedpe, 'rb') as fin:

        for line in fin.readlines():
            columns = line.decode('utf8').rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])
            svtype = columns[6]

            srpos.append((chrom1, pos1_start, pos1_end, chrom2,
                            pos2_start, pos2_end, svtype))

    # print('counting {} split read positions'.format(len(srpos)))
    # print(sv_list)

    return srpos


def main():

    args = parse_cl_args(sys.argv[1:])

    srpos = get_srpos_from_bedpe(args.split_reads)

    sv_list = read_vcf(args.input)
    trees_start, trees_end = make_gtrees_from_svlist(sv_list)

    # print(srpos)
    srpos_reduced = [(c[0], c[1], c[3], c[4], c[6]) for c in srpos]
    lookup_start, lookup_end = search_tree_with_cpos(srpos_reduced, trees_start, trees_end, int(args.win/2))

    idx = []
    for i, x in enumerate(zip(lookup_start, lookup_end)):

        lu_start, lu_end = x

        l1 = len(lu_start)
        l2 = len(lu_end)

        if l1 > 1 and l2 > 1:

            lu_start_set = set()
            lu_end_set = set()

            for s in lu_start:
                lu_start_elem_start, lu_start_elem_end, lu_start_elem_data = s
                lu_start_elem_svtype, lu_start_elem_svid = lu_start_elem_data
                lu_start_set.add(lu_start_elem_svid)

            for s in lu_end:
                lu_end_elem_start, lu_end_elem_end, lu_end_elem_data = s
                lu_end_elem_svtype, lu_end_elem_svid = lu_end_elem_data
                lu_end_set.add(lu_end_elem_svid)

        if len(lu_start_set & lu_end_set) > 0:
            idx.append(i)

    # print(idx)

    out = open(args.output, 'w')
    out_nosr = open(args.output_nosr, 'w')
    i = 0

    with open(args.input, 'r') as fin:

        for line in fin.readlines():
            if line[0] == '#':
                out.write(line)
                out_nosr.write(line)
            else:
                if i in idx:
                    out.write(line)
                else:
                    out_nosr.write(line)
                i+=1

    out.close()
    out_nosr.close()


if __name__ == "__main__":
    main()
