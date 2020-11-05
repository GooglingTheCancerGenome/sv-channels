import argparse
import logging
import os
import numpy as np
import pysam
from time import time

from functions import get_config_file, get_insert_size


def is_read_proper_paired(read, bam_mean, bam_stddev):

    dist = read.next_reference_start - read.reference_start

    if dist > 0 and abs(dist) <= bam_mean + bam_stddev:
        return True
    else:
        return False


def check_read(read):
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


def check_read_is_proper_paired_forward(read, bam_mean, bam_stddev):
    '''

    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if not read.mate_is_unmapped \
            and is_read_proper_paired(read, bam_mean, bam_stddev) \
            and not read.is_reverse:
        return False

    return True


def check_read_is_proper_paired_reverse(read, bam_mean, bam_stddev):
    '''

    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if not read.mate_is_unmapped \
            and is_read_proper_paired(read, bam_mean, bam_stddev) \
            and read.is_reverse:
        return False

    return True


def get_coverage(ibam, chrName, minMAPQ, outFile):
    '''
    This function fills the coverage array for the chromosome
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the coverage array
    :return: None
    '''

    # Open BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    bam_mean, bam_stddev = get_insert_size(ibam, bamfile, minMAPQ)

    # Get chromosome length from BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen

    # Numpy array to store the coverage
    cov = np.zeros((3, chrLen), dtype=np.uint32)

    read_quality_sum = np.zeros(chrLen, dtype=np.uint32)
    read_quality_count = np.zeros(chrLen, dtype=np.uint32)

    # Log information every n_r base pair positions
    n_r = 10 ** 6
    # print(n_r)
    last_t = time()
    # print(type(last_t))

    # Pysam iterator to fetch the reads
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    for i, read in enumerate(iter, start=1):

        # Every n_r alignments, write log informations
        if not i % n_r:
            # Record the current time
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" %
                         (i, n_r / (now_t - last_t)))

            last_t = time()

        if not read.is_unmapped and read.mapping_quality >= minMAPQ:
            if check_read(read):
                cov[0, read.reference_start:read.reference_end - 1] += 1
            if check_read_is_proper_paired_forward(read, bam_mean, bam_stddev):
                cov[1, read.reference_start:read.reference_end - 1] += 1
            if check_read_is_proper_paired_reverse(read, bam_mean, bam_stddev):
                cov[2, read.reference_start:read.reference_end - 1] += 1

            # add read mapping quality
            read_quality_sum[read.reference_start:read.reference_end -
                                                  1] += read.mapping_quality
            read_quality_count[read.reference_start:read.reference_end -
                                                    1] += 1

    read_quality = np.divide(read_quality_sum,
                             read_quality_count,
                             where=read_quality_count != 0)
    # where there are no reads, use median mapping quality
    read_quality[np.where(read_quality_count == 0)] = np.median(read_quality)

    cov = np.vstack((cov, read_quality))
    logging.info(cov.shape)

    for i in np.arange(cov.shape[0]):
        logging.info('chromosome {} coverage:'+ \
                     'non-zero elements at index {}:{}'.format(chrName,
                                                               i,
                                                               np.argwhere(cov[i, :] != 0).shape[0]))

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

    # Default chromosome is 12 for the artificial data

    parser = argparse.ArgumentParser(description='Create coverage channel')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chr',
                        type=str,
                        default='12',
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

    get_coverage(ibam=args.bam, chrName=args.chr, minMAPQ=args.min_mapq, outFile=output_file)

    logging.info('Time: coverage on BAM %s and Chr %s: %f' %
                 (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
