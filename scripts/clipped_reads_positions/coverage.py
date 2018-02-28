import argparse
import pysam
import numpy as np
from time import time
import functions as fun

def get_coverage(bamfile, chr_name, win_center, window_size=fun.window_size):

    #fun.logging.debug('position: %d', win_center)

    cov = np.zeros(window_size)
    start_pos = win_center - int(window_size/2)
    stop_pos = win_center + int(window_size/2)

    for pile in bamfile.pileup(chr_name, start_pos, stop_pos, truncate=True):
        cov[pile.pos - start_pos] = pile.n

    return cov

def get_coverage_list(ibam, chr_name, outFile):

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    chr_len = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chr_name][0]

    clipped_pos = fun.get_clipped_reads_positions(chr_name)

    cov_list = []

    i = 1
    for win_center in clipped_pos:
        if int(fun.window_size / 2) <= win_center <= (chr_len - int(fun.window_size / 2)):
            if i % 10000 == 0:
                fun.logging.debug('count: %d', i)
            cov_list.append(get_coverage(bamfile, chr_name, win_center))
            i += 1

    #print(cov_list)
    # cPickle data persistence
    with fun.bz2.BZ2File(outFile, 'w') as f:
        fun.pickle.dump(cov_list, f)


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
    get_coverage_list(ibam=args.bam, chr_name=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
