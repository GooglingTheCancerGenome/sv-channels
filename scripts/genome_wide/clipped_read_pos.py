import argparse
import logging
import os
import pickle
import bz2file
from collections import Counter
from time import time

import pysam

import functions as fun


def get_clipped_read_positions(ibam, chrName, outFile):

    assert os.path.isfile(ibam)

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    clipped_pos = []
    clipped_read_1 = set()
    clipped_read_2 = set()

    def add_clipped_read(read):
        if read.is_read1:
            clipped_read_1.add(read.query_name)
        if read.is_read2:
            clipped_read_2.add(read.query_name)

    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    i = 0
    n_r = 10 ** 6
    # print(n_r)
    last_t = time()
    # print(type(last_t))
    for i, read in enumerate(iter, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        if (not read.is_unmapped) and (not read.mate_is_unmapped) and len(read.get_reference_positions()) > 0:
            # assert read.reference_name in clipped_pos.keys()
            # assert read.cigartuples[0][0] in [4, 5] or read.cigartuples[-1][0] in [4, 5]
            if fun.is_left_clipped(read):
                # add_clipped_read(read)
                cpos = read.get_reference_positions()[0] + 1
                clipped_pos.append(cpos)
            if fun.is_right_clipped(read):
                # add_clipped_read(read)
                cpos = read.get_reference_positions()[-1] + 1
                clipped_pos.append(cpos)

    #clipped_pos = list(set(clipped_pos))
    bamfile.close()

    clipped_pos_cnt = Counter(clipped_pos)

    # cPickle data persistence
    outFile_base = os.path.splitext(os.path.basename(outFile))[0]
    with bz2file.BZ2File(outFile, 'wb') as f:
        # obj = (clipped_pos_cnt, clipped_read_1, clipped_read_2)
        pickle.dump(clipped_pos_cnt, f)
    #with bz2file.BZ2File(outFile_base + '_cr.pbz2', 'wb') as f:
        # obj = (clipped_pos_cnt, clipped_read_1, clipped_read_2)
    #    pickle.dump((clipped_read_1, clipped_read_2), f)


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.cr.bam"

    parser = argparse.ArgumentParser(description='Get clipped reads positions')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='clipped_read_pos.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='clipped_read_pos.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    get_clipped_read_positions(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()