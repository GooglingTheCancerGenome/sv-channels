import argparse
import pysam
import bz2
import cPickle as pickle
from time import time


def get_coverage(ibam, chrName, outFile):

    cov = []

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen

    for pile in bamfile.pileup(chrName, start_pos , stop_pos, truncate=True):
        while pile.pos != start_pos:
            cov.append(0)
            start_pos = start_pos + 1
        cov.append(pile.n)
        start_pos = start_pos + 1

    # cPickle data persistence
    with bz2.BZ2File(outFile, 'w') as f:
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

    args = parser.parse_args()

    t0 = time()
    get_coverage(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time()-t0)

if __name__ == '__main__':
    main()