# Imports

import argparse
import pysam
import bz2file
from time import time
import logging
import numpy as np

minMAPQ = 30


def check_read(read):
    '''

    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if read.reference_name == read.next_reference_name and read.mapping_quality >= minMAPQ \
            and read.is_reverse != read.mate_is_reverse:
        return True

    return False
    #return True

def get_coverage(ibam, chrName, outFile):
    '''
    This function fills the coverage array for the chromosome
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the coverage array
    :return: None
    '''

    # Open BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Get chromosome length from BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen

    # Numpy array to store the coverage
    cov = np.zeros(chrLen, dtype=np.uint32)

    # Log information every n_r base pair positions
    n_r = 10 ** 6
    # print(n_r)
    last_t = time()
    # print(type(last_t))

    # # Iterate over the chromosome positions
    # for i, pile in enumerate(bamfile.pileup(chrName, start_pos, stop_pos, truncate=True), start=1):
    #
    #     if not i % n_r:
    #         now_t = time()
    #         # print(type(now_t))
    #         logging.info("%d pileup positions processed (%f positions / s)" % (
    #             i,
    #             n_r / (now_t - last_t)))
    #         last_t = time()
    #     #print('Pos: %d, Cov: %d' % (pile.pos, pile.n))
    #     try:
    #         cov[pile.pos] = pile.n
    #     except MemoryError:
    #         logging.info("Out of memory for chr %s and BAM file %s !" % (chrName, ibam))

    # Replacing pileup with count_coverage
    cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos, read_callback=check_read)
    cov = np.asarray(cov_A, dtype=np.uint32) + \
          np.asarray(cov_C, dtype=np.uint32) + \
          np.asarray(cov_G, dtype=np.uint32) + \
          np.asarray(cov_T, dtype=np.uint32)
    # print(cov)

    # cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos)
    # cov_nofilter = np.asarray(cov_A, dtype=int) + \
    #       np.asarray(cov_C, dtype=int) + \
    #       np.asarray(cov_G, dtype=int) + \
    #       np.asarray(cov_T, dtype=int)
    # print(cov_nofilter)
    #
    # assert np.all(cov == cov_nofilter)

    # Save coverage numpy array
    with bz2file.BZ2File(outFile, 'w') as f:
        try:
            np.save(file=f, arr=cov)
        except MemoryError:
            logging.info("Out of memory for chr %s and BAM file %s !" % (chrName, ibam))


def main():
    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    # Default chromosome is 17 for the artificial data

    parser = argparse.ArgumentParser(description='Create coverage channel')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='coverage.npy.bz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='coverage.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    get_coverage(ibam=args.bam, chrName=args.chr, outFile=args.out)
    logging.info('Time: coverage on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
