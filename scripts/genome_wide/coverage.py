import argparse
import pysam
import bz2file
import pickle
from time import time
import logging
import numpy as np

def get_coverage(ibam, chrName, outFile):

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    cov = dict()

    i = 0
    n_r = 10 ** 6
    # print(n_r)
    last_t = time()
    # print(type(last_t))

    for i, pile in enumerate(bamfile.pileup(chrName, start_pos, stop_pos, truncate=True), start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d pileup positions processed (%f positions / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()
        #print('Pos: %d, Cov: %d' % (pile.pos, pile.n))
        cov[pile.pos] = pile.n

        #while pile.pos != start_pos:
        #    cov.append(0)
        #    start_pos = start_pos + 1
        #cov.append(pile.n)
        #start_pos = start_pos + 1

    # cPickle data persistence

    with bz2file.BZ2File(outFile, 'w') as f:
        pickle.dump(cov, f)


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    parser = argparse.ArgumentParser(description='Create coverage channel')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='coverage.pbz2',
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
    print(time() - t0)


if __name__ == '__main__':
    main()
