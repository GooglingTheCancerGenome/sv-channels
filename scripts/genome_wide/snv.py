import argparse
import logging
import os
from collections import Counter
from time import time

import numpy as np
import pysam
import twobitreader as twobit

from functions import *


def get_snvs(ibam, itwobit, chrName, max_coverage, outFile):

    def get_snv_number(query_seq_list, reference_base):

        reference_base = reference_base.upper()
        if len(query_seq_list) > 0 and reference_base != 'N':
            cnt = Counter(list(map(lambda x: x.upper(), query_seq_list)))
            return cnt['A'] + cnt['T'] + cnt['C'] + cnt['G'] - cnt[reference_base]
        return 0

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]
    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    reference_sequence = twobit.TwoBitFile(itwobit)
    snv_list = ['BQ', 'SNV', 'MAPQ']
    snv_array = np.zeros(shape=(stop_pos, len(snv_list)), dtype=np.float32)
    snv_dict = {v: n for n, v in enumerate(snv_list)}

    for pileupcolumn in bamfile.pileup(chrName,
                                       start_pos,
                                       stop_pos,
                                       stepper='all'):
        if 0 < pileupcolumn.nsegments < max_coverage and start_pos <= pileupcolumn.pos <= stop_pos:
            quals = pileupcolumn.get_query_qualities()
            if len(quals) > 0:
                snv_array[pileupcolumn.pos, snv_dict['BQ']] = np.median(
                    quals)
            quals = pileupcolumn.get_mapping_qualities()
            if len(quals) > 0:
                snv_array[pileupcolumn.pos, snv_dict['MAPQ']] = np.median(
                    quals)
            try:
                query_seq_list = pileupcolumn.get_query_sequences()
                snv_number = get_snv_number(
                    query_seq_list,
                    reference_sequence[chrName][pileupcolumn.pos])
                snv_array[pileupcolumn.pos, snv_dict['SNV']] = snv_number / pileupcolumn.nsegments \
                    if pileupcolumn.nsegments != 0 else 0

            except AssertionError as error:
                # Output expected AssertionErrors.
                logging.info(error)
                logging.info("Position %s:%d has %d nsegments" % (
                    str(chrName), pileupcolumn.pos, pileupcolumn.nsegments))
                continue

    for i in np.arange(snv_array.shape[1]):
        logging.info("snv array %s: non-zero elements at index %d:%d" %
                     (str(chrName), i, np.argwhere(snv_array[:, i] != 0).shape[0]))

    # Write the output
    np.save(file=outFile, arr=snv_array)
    os.system('gzip -f ' + outFile)


def main():
    parser = argparse.ArgumentParser(description='Get SNV info')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-t',
                        '--twobit',
                        type=str,
                        default='../../data/test.2bit',
                        help="Specify input file (2bit)")
    parser.add_argument('-c',
                        '--chr',
                        type=str,
                        default='12',
                        help="Specify chromosome")
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        default='snv.npy',
                        help="Specify output")
    parser.add_argument('-p',
                        '--outputpath',
                        type=str,
                        default='.',
                        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='snv.log',
                        help='File in which to write logs.')
    parser.add_argument('-pb',
                        '--max_coverage',
                        type=int,
                        default=1000,
                        help='Consider only regions with coverage less than max_coverage to speed up the processing')
    args = parser.parse_args()
    cmd_name = 'snv'
    output_dir = os.path.join(args.outputpath, cmd_name)
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, '_'.join((args.chr, args.logfile)))
    output_file = os.path.join(output_dir, '_'.join((args.chr, args.out)))
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)

    t0 = time()
    get_snvs(ibam=args.bam,
             itwobit=args.twobit,
             chrName=args.chr,
             max_coverage=args.max_coverage,
             outFile=output_file)
    logging.info('Time: SNVs on BAM %s and Chr %s: %f' %
                 (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
