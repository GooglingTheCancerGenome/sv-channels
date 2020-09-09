import argparse
import logging
import os
import numpy as np
import pysam
from time import time
from zlib import crc32
from functions import get_config_file

config = get_config_file()
minMAPQ = config["DEFAULT"]["MIN_MAPQ"]


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


def check_read_is_proper_paired_forward(read):
    '''

    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ \
            and not read.is_proper_pair and not read.is_reverse:
        return True

    return False


def check_read_is_proper_paired_reverse(read):
    '''

    :param read: AlignedSegment
    :return: True if all these conditions are valid:
        - read and mate are mapped on the same chromosome,
        - read mapping quality is greater than minMAPQ,
        - read and mate are mapped on opposite strands
    '''

    if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ \
            and not read.is_proper_pair and read.is_reverse:
        return True

    return False


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
    cov = np.zeros((5, chrLen), dtype=np.uint32)

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

        if check_read(read):
            cov[0, read.reference_start:read.reference_end - 1] += 1
        if check_read_is_proper_paired_forward(read):
            cov[1, read.reference_start:read.reference_end - 1] += 1
        if check_read_is_proper_paired_reverse(read):
            cov[2, read.reference_start:read.reference_end - 1] += 1

        if not read.is_unmapped and read.mapping_quality >= minMAPQ:

            # add read mapping quality
            read_quality_sum[read.reference_start:read.reference_end -
                                                  1] += read.mapping_quality
            read_quality_count[read.reference_start:read.reference_end -
                                                    1] += 1

            # using hash of query name
            j = int((read.reference_end - read.reference_start) / 2)
            # print(j)
            h = crc32(read.query_name.encode('utf8')) & 0xffffffff
            # print(h)
            if not read.is_reverse:
                cov[3, j] = h if not cov[3, j] else 0
            else:
                cov[4, j] = h if not cov[4, j] else 0

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
    # cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos, read_callback=check_read)
    # cov = np.asarray(cov_A, dtype=int) + \
    #       np.asarray(cov_C, dtype=int) + \
    #       np.asarray(cov_G, dtype=int) + \
    #       np.asarray(cov_T, dtype=int)
    #
    # cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos,
    #                                                     read_callback=check_read_is_proper_paired_forward)
    # cov_disc_f = np.asarray(cov_A, dtype=int) + \
    #       np.asarray(cov_C, dtype=int) + \
    #       np.asarray(cov_G, dtype=int) + \
    #       np.asarray(cov_T, dtype=int)
    #
    # cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos,
    #                                                     read_callback=check_read_is_proper_paired_reverse)
    # cov_disc_r = np.asarray(cov_A, dtype=int) + \
    #       np.asarray(cov_C, dtype=int) + \
    #       np.asarray(cov_G, dtype=int) + \
    #       np.asarray(cov_T, dtype=int)
    #
    # cov = np.vstack((cov, cov_disc_f, cov_disc_r))

    # print(cov)

    # cov_A, cov_C, cov_G, cov_T = bamfile.count_coverage(chrName, start_pos, stop_pos)
    # cov_nofilter = np.asarray(cov_A, dtype=int) + \
    #       np.asarray(cov_C, dtype=int) + \
    #       np.asarray(cov_G, dtype=int) + \
    #       np.asarray(cov_T, dtype=int)
    # print(cov_nofilter)
    #
    # assert np.all(cov == cov_nofilter)

    read_quality = np.divide(read_quality_sum,
                             read_quality_count,
                             where=read_quality_count != 0)
    # where there are no reads, use median mapping quality
    read_quality[np.where(read_quality_count == 0)] = np.median(read_quality)

    cov = np.vstack((cov, read_quality))
    logging.info(cov.shape)

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
    get_coverage(ibam=args.bam, chrName=args.chr, outFile=output_file)
    logging.info('Time: coverage on BAM %s and Chr %s: %f' %
                 (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
