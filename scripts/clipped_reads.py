import argparse
import pysam
import bz2
import cPickle as pickle
from time import time
import twobitreader as twobit
from collections import Counter

# Return if a read is clipped on the left
def is_left_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False

# Return if a read is clipped on the right
def is_right_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False

def get_clipped_reads(ibam, chrName, outFile):

    clipped_left = dict()
    clipped_right = dict()

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    #print([i['LN'] for i in header_dict['SQ'] if i['SN']])
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    #print(chrLen)

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    for read in iter:
        if not read.is_unmapped and len(read.get_reference_positions()) > 0 and read.cigartuples is not None:
            if is_left_clipped(read):
                #print(str(read))
                #print('Clipped at the start: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                #print('Pos:%d, clipped_pos:%d' % (read.reference_start, read.get_reference_positions()[0]))
                ref_pos = read.get_reference_positions()[0] + 1
                if ref_pos not in clipped_left.keys():
                    clipped_left[ref_pos] = 1
                else:
                    clipped_left[ref_pos] += 1
            if is_right_clipped(read):
                #print('Clipped at the end: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                #print('Pos:%d, clipped_pos:%d' %(read.reference_end, read.get_reference_positions()[-1]))
                ref_pos = read.get_reference_positions()[-1] + 1
                if ref_pos not in clipped_left.keys():
                    clipped_right[ref_pos] = 1
                else:
                    clipped_right[ref_pos] += 1

    # cPickle data persistence
    with bz2.BZ2File(outFile, 'w') as f:
        pickle.dump((clipped_left, clipped_right), f)


def main():

    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    parser = argparse.ArgumentParser(description='Create channels with number of left/right clipped reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='clipped_reads.pbz2',
                        help="Specify output")

    args = parser.parse_args()

    t0 = time()
    get_clipped_reads(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time()-t0)

if __name__ == '__main__':
    main()