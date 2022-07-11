#!/usr/bin/env python3

"""
Adapted from https://git.wur.nl/wijfj001/hecaton/-/blob/master/scripts/convert/bedpe_to_vcf.py
Convert BEDPE file containing SV calls to VCF
"""

import argparse
import datetime
import os
import sys

import pandas as pd
import twobitreader as twobit

working_dir = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/manuscript/Figures/F2/'

rev_bases = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert BEDPE file containing CNV calls to VCF"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_bedpe", type=str,
                        default=os.path.join(
                            working_dir, 'simulated-data-results/coverage/cov30/hmz-sv/sv-channels.split_reads.bedpe'),
                        # default=os.path.join(working_dir,
                        #                     'simulated-data-results/coverage/cov30/hmz-sv/hmz-sv.bedpe'),
                        help="Path to BEDPE file")
    parser.add_argument("-b", "--twobit", type=str,
                        default=os.path.join(
                            working_dir, 'helpers/aux_files/hs37d5.2bit'),
                        help="Path to 2bit file")
    parser.add_argument("-o", "--output_vcf", type=str,
                        default=os.path.join(
                            working_dir, 'simulated-data-results/coverage/cov30/hmz-sv/sv-channels.split_reads.vcf'),
                        # default=os.path.join(working_dir,
                        #                     'simulated-data-results/coverage/cov30/hmz-sv/hmz-sv.new.vcf'),
                        help="Name of VCF output file")
    parser.add_argument("-s", "--sample_name", type=str,
                        default='HTZ-SV',
                        help="Name of the sample")
    args = parser.parse_args(in_args)
    return args


def convert_bedpe(input_bedpe_fn, output_vcf_fn, sample_name, genome):
    """Convert BEDPE file containing SV calls to VCF

    :param input_bedpe_fn: Path to input BEDPE file
    :param output_vcf_fn: Path to output VCF file
    :param sample_name: Name of the sample
    :return: 0 (integer)
    """

    # write VCF header to output file
    vcf_header_elems = ["##fileformat=VCFv4.2",
                        "##fileDate={}".format(
                            datetime.date.today().strftime("%Y%m%d")),
                        "##FILTER=<ID=PASS,Description=\"All filters passed\">",
                        "##ALT=<ID=DEL,Description=\"Deletion\">",
                        "##ALT=<ID=INS,Description=\"Novel sequence insertion\">",
                        "##ALT=<ID=INV,Description=\"Inversion\">",
                        "##ALT=<ID=DUP,Description=\"Duplication\">",
                        "##ALT=<ID=BND,Description=\"Inter-chromosomal translocation\">",
                        "##FORMAT=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                        "##FORMAT=<ID=CHR2,Number=1,Type=Float,Description=\"Stop chromosome of the interval\">",
                        "##FORMAT=<ID=END,Number=1,Type=Float,Description=\"Stop position of the interval\">",
                        "##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of partner breakend\">",
                        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the interval\">",
                        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=' +
                        '"Difference in length between REF and ALT alleles">',
                        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">"]

    for k in sorted(genome.keys()):
        vcf_header_elems.append(
            "##contig=<ID={},length={}>".format(k, len(genome[k])))

    vcf_header_elems.append(
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample_name))
    vcf_header = "\n".join(vcf_header_elems)
    output_vcf = open(output_vcf_fn, "w")
    output_vcf.write(vcf_header)
    output_vcf.write("\n")
    calls_index = 1

    if os.path.exists(input_bedpe_fn):
        input_bedpe = pd.read_csv(input_bedpe_fn, sep="\t", header=None)
    else:
        print('{} does not exist!'.format(input_bedpe_fn))
        sys.exit(1)

    # loop through bedpe file writing each variant to the vcf
    for sv_calls in input_bedpe.itertuples(index=False, name='Pandas'):
        chrom = str(sv_calls[0])
        start = str(sv_calls[1])
        chrom2 = str(sv_calls[3])
        end = str(sv_calls[4])
        svlen = sv_calls[4] - sv_calls[1]
        svtype = str(sv_calls[6])
        svtype = svtype.split('_')[1] if svtype not in [
            'DEL', 'INS', 'INV', 'DUP', 'CTX', 'TRA'] else svtype
        qual = str(sv_calls[7]) if len(sv_calls) > 7 else '1'
        filtered = "PASS"
        # extract format field elements
        format_field = "GT"
        sample_field_elems = ["./."]
        sample_field = ":".join(sample_field_elems)
        ref = genome[str(chrom)][sv_calls[1]].upper()

        if svtype not in ['TRA', 'CTX']:
            identifier = sample_name+'_' + str(calls_index)
            alt = '<' + svtype + '>'
            if svtype == 'INS':
                info_field = "SVTYPE={};END={};SVLEN={}".format(
                    svtype, str(sv_calls[1]+1), svlen)
            else:
                info_field = "SVTYPE={};END={};SVLEN={}".format(
                    svtype, end, svlen)

            # create new line for variant
            variant_line_elems = [chrom, start, identifier, ref, alt, qual,
                                  filtered, info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
            calls_index += 1
        else:
            def get_identifier(calls_index, j):
                identifier = '_'.join((sample_name, str(calls_index), str(j)))
                return identifier
            # Writing a breakpoint of an inter-chromosomal translocation as in BND format
            j = 1
            bp1_pos = sv_calls[1]
            bp2_pos = sv_calls[4]
            svlen = abs(bp2_pos-bp1_pos) if chrom == chrom2 else -1
            brkt_fw = ']'
            brkt_bw = '['
            bp1_ref_fw = genome[str(chrom)][bp1_pos].upper()
            bp1_ref_rv = rev_bases[bp1_ref_fw]
            bp1_alt_1 = bp1_ref_fw + brkt_bw + \
                chrom2 + ':' + str(bp2_pos) + brkt_bw
            bp1_alt_2 = brkt_fw + chrom2 + ':' + \
                str(bp2_pos) + brkt_fw + bp1_ref_rv
            # bp1_1
            bp1_1_info_field = "SVTYPE={};MATEID={}".format(
                'BND', get_identifier(calls_index, 4))
            variant_line_elems = [chrom, start, get_identifier(calls_index, j), bp1_ref_fw, bp1_alt_1, qual,
                                  filtered, bp1_1_info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
            j += 1
            # bp1_2
            bp1_2_info_field = "SVTYPE={};MATEID={}".format(
                'BND', get_identifier(calls_index, 3))
            variant_line_elems = [chrom, start, get_identifier(calls_index, j), bp1_ref_rv, bp1_alt_2, qual,
                                  filtered, bp1_2_info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
            j += 1
            bp2_ref_fw = genome[str(chrom2)][bp2_pos].upper()
            bp2_ref_rv = rev_bases[bp2_ref_fw]
            bp2_alt_1 = bp2_ref_rv + brkt_bw + \
                chrom + ':' + str(bp1_pos) + brkt_bw
            bp2_alt_2 = brkt_fw + chrom + ':' + \
                str(bp1_pos) + brkt_fw + bp1_ref_fw
            # bp2_1
            bp2_1_info_field = "SVTYPE={};MATEID={}".format(
                'BND', get_identifier(calls_index, 2))
            variant_line_elems = [chrom2, end, get_identifier(calls_index, j), bp2_ref_rv, bp2_alt_1, qual,
                                  filtered, bp2_1_info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
            j += 1
            # bp2_2
            bp2_2_info_field = "SVTYPE={};MATEID={}".format(
                'BND', get_identifier(calls_index, 1))
            variant_line_elems = [chrom2, end, get_identifier(calls_index, j), bp2_ref_fw, bp2_alt_2, qual,
                                  filtered, bp2_2_info_field, format_field, sample_field]
            variant_line = "\t".join(variant_line_elems)
            output_vcf.write(variant_line)
            output_vcf.write("\n")
            j += 1
            calls_index += 1
    output_vcf.close()
    return 0


def main():
    args = parse_cl_args(sys.argv[1:])
    genome = twobit.TwoBitFile(args.twobit)
    convert_bedpe(args.input_bedpe, args.output_vcf, args.sample_name, genome)


if __name__ == "__main__":
    main()
