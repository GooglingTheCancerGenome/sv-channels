#!/usr/bin/env python3

"""
Adapted from https://git.wur.nl/wijfj001/hecaton/-/blob/master/scripts/convert/bedpe_to_vcf.py
Convert BEDPE file containing SV calls to VCF
"""

import argparse
import datetime
import pandas as pd
import sys
import os


def parse_cl_args(in_args):
    """
    Parse command line arguments

    :param in_args: All command line arguments
    :return: None
    """
    description = "Convert BEDPE file containing CNV calls to VCF"
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i", "--input_bedpe", type=str,
                        default='',
                        help="Path to BEDPE file")
    parser.add_argument("-o", "--output_vcf", type=str,
                        default='',
                        help="Name of VCF output file")
    parser.add_argument("-s", "--sample_name", type=str,
                        default='HTZ-SV',
                        help="Name of the sample")
    args = parser.parse_args(in_args)
    return args


def convert_bedpe(input_bedpe_fn, output_vcf_fn, sample_name):
    """Convert BEDPE file containing SV calls to VCF

    :param input_bedpe_fn: Path to input BEDPE file
    :param output_vcf_fn: Path to output VCF file
    :param sample_name: Name of the sample
    :return: 0 (integer)
    """

    # write VCF header to output file
    vcf_header_elems = ["##fileformat=VCFv4.3",
                        "##fileDate={}".format(datetime.date.today().strftime("%Y%m%d")),
                        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
                        "##INFO=<ID=END,Number=.,Type=Integer,Description=\"End position of the variant described in this region\">",
                        "##ALT=<ID=DEL,Description=\"Deletion\">",
                        "##ALT=<ID=INS,Description=\"Insertion of novel sequence\">",
                        "##ALT=<ID=INV,Description=\"Inversion\">",
                        "##ALT=<ID=DUP,Description=\"Duplication\">",
                        "##ALT=<ID=BND,Description=\"Breakend\">",
                        "##FORMAT=<ID=CHROM2,Number=1,Type=Float,Description=\"Breakpoint2 chromosome\">",
                        "##FORMAT=<ID=END,Number=1,Type=Float,Description=\"Breakpoint2 position\">",
                        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}".format(sample_name)]

    vcf_header = "\n".join(vcf_header_elems)

    with open(output_vcf_fn, "w") as output_vcf:
        output_vcf.write(vcf_header)
        output_vcf.write("\n")
        # load bedpe file
        calls_index = 1
        for t in ['DEL', 'INS', 'INV', 'DUP', 'BND']:
            fn = input_bedpe_fn + '_' + t + '.bedpe'

            if os.path.exists(fn):

                input_bedpe = pd.read_csv(input_bedpe_fn+'_'+t+'.bedpe', sep="\t", header=None)

                # loop through bedpe file writing each variant to the vcf
                for sv_calls in input_bedpe.itertuples(index=False, name='Pandas'):

                    # print(sv_calls)

                    chrom = str(sv_calls[0])
                    start = str(sv_calls[1])
                    chrom2 = str(sv_calls[3])
                    end = str(sv_calls[4])
                    qual = str(sv_calls[7])
                    identifier = 'sv-channels_'+str(calls_index)
                    calls_index += 1
                    ref ='N'
                    alt='<'+t+'>'
                    filtered = "PASS"
                    info_field = "SVTYPE={};CHR2={};END={}".format(t, chrom2, end)
                    # extract format field elements
                    format_field = "GT"
                    sample_field_elems = ["1/1"]
                    sample_field = ":".join(sample_field_elems)
                    # create new line for variant
                    variant_line_elems = [chrom, start, identifier, ref, alt, qual,
                                           filtered, info_field, format_field, sample_field]
                    variant_line = "\t".join(variant_line_elems)
                    output_vcf.write(variant_line)
                    output_vcf.write("\n")
    return 0


def main():
    args = parse_cl_args(sys.argv[1:])
    # convert vcf file
    convert_bedpe(args.input_bedpe, args.output_vcf, args.sample_name)


if __name__ == "__main__":
    main()


