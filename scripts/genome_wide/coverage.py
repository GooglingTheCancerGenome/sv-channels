import argparse
import logging
import os
from time import time

import numpy as np
import pysam

from functions import get_insert_size


def is_discordant(read, bam_mean, bam_stddev):
    dist = read.next_reference_start - read.reference_start
    if dist > 0 and abs(dist) <= bam_mean + bam_stddev:
        return False
    return True


def is_properly_mapped(read):
    '''
    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if read.reference_name == read.next_reference_name \
            and read.is_reverse != read.mate_is_reverse:
        return True

    return False


def get_coverage(ibam, chrName, minMAPQ, outFile):
    '''
    This function fills the coverage array for the chromosome
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the coverage array
    :return: None
    '''
    bamfile = pysam.AlignmentFile(ibam, "rb")
    bam_mean, bam_stddev = get_insert_size(ibam, bamfile, minMAPQ)
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]
    start_pos = 0
    stop_pos = chrLen
    cov = np.zeros((chrLen, 5), dtype=np.int)
    n_r = 10 ** 6
    last_t = time()

    for i, read in enumerate(bamfile.fetch(chrName, start_pos, stop_pos), start=1):
        # Every n_r alignments, write log informations
        if not i % n_r:
            logging.info("%d alignments processed (%f alignments / s)" %
                         (i, n_r / (time() - last_t)))
            last_t = time()

        if not read.is_unmapped and read.mapping_quality >= minMAPQ:
            if is_properly_mapped(read):
                cov[read.reference_start:read.reference_end - 1, 0] += 1
            read_discordant = is_discordant(read, bam_mean, bam_stddev)

            if not read.mate_is_unmapped:
                if read_discordant:
                    if read.is_reverse:
                        cov[read.reference_start:read.reference_end - 1, 2] += 1
                    else:
                        cov[read.reference_start:read.reference_end - 1, 1] += 1
                if not read.is_proper_pair:
                    if read.is_reverse:
                        cov[read.reference_start:read.reference_end - 1, 4] += 1
                    else:
                        cov[read.reference_start:read.reference_end - 1, 3] += 1
    logging.info(cov.shape)

    for i in np.arange(cov.shape[1]):
        logging.info("chromosome %s coverage: non-zero elements at index %d:%d" %
                     (chrName, i, np.argwhere(cov[i, :] != 0).shape[0]))
        logging.info("mean:%f, sd:%f" %
                     (np.mean(cov[:, i]), np.std(cov[:, i])))

    # Save coverage numpy array
    try:
        np.save(file=outFile, arr=cov)
    except MemoryError:
        logging.info("Out of memory for chr %s and BAM file %s !" %
                     (chrName, ibam))
    os.system('gzip -f ' + outFile)
    # To load it
    # cov = np.load(outFile)['coverage']


def main():
    parser = argparse.ArgumentParser(description='Create coverage channel')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chr',
                        type=str,
                        default='22',
                        help="Specify chromosome")
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        default='coverage.npy',
                        help="Specify output")
    parser.add_argument('-m',
                        '--min_mapq',
                        type=int,
                        default=10,
                        help='Minimum read mapping quality')
    parser.add_argument(
        '-p',
        '--outputpath',
        type=str,
        default='.',
        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='coverage.log',
                        help='File in which to write logs.')
    args = parser.parse_args()
    cmd_name = 'coverage'
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
    get_coverage(ibam=args.bam, chrName=args.chr,
                 minMAPQ=args.min_mapq, outFile=output_file)
    logging.info('Time: coverage on BAM %s and Chr %s: %f' %
                 (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
