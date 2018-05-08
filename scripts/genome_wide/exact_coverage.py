'''

This file contains the functions to compute the channel for the exact coverage. In the current implementation of
ChannelMaker, the exact coverage was discarded because considered unnecessary, given that the coverage channel is
already present. Computing the exact coverage requires too much time.

'''

#Imports
import argparse
import pysam
import bz2file
import pickle
from time import time
import twobitreader as twobit
from collections import Counter
import logging


def get_exact_coverage(ibam, chrName, outFile):
    '''

        :param ibam: input BAM alignment file
        :param chrName: chromosome name to consider
        :param outFile: output file with the list containing the exact coverage per chromosome position. Exact coverage
        is the number of bases mapped in a certain position that are equal to the reference base.
        :return: None. It only writes the exact coverage in the output file
    '''

    # Path on the local machine of the 2bit version of the human reference genome (hg19)
    # genome = twobit.TwoBitFile('/Users/lsantuari/Documents/Data/GiaB/reference/hg19.2bit')
    # Path on the HPC of the 2bit version of the human reference genome (hg19)
    genome = twobit.TwoBitFile('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes/hg19.2bit')

    # List used to store the exact coverage
    cov = []

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # print(chrLen)

    # Print every n_r alignments processed
    n_r = 10 ** 6
    # print(n_r)
    # Record the current time
    last_t = time()

    # Iterate over the chromosome positions
    for i, pile in enumerate(bamfile.pileup(chrName, start_pos, stop_pos, truncate=True), start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d pileup positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        while pile.pos != start_pos:
            cov.append(0)
            start_pos = start_pos + 1

        # Calculate exact coverage
        excov = Counter([pileupread.alignment.query_sequence[pileupread.query_position]
                         for pileupread in pile.pileups if not pileupread.is_del and not pileupread.is_refskip
                         ])[genome['chr' + chrName][pile.pos].upper()]
        cov.append(excov)

        start_pos = start_pos + 1

    # Save the exact coverage list
    with bz2file.BZ2File(outFile, 'w') as f:
        pickle.dump(cov, f)


def main():

    # Default BAM file for testing
    wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    inputBAM = wd + "T0_dedup.bam"

    parser = argparse.ArgumentParser(description='Create exact coverage channel')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='exact_coverage.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='exact_coverage.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    # Compute the exact coverage
    get_exact_coverage(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
