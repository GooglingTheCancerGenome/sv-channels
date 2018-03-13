import argparse
import pysam
import bz2file
import pickle
from time import time
import logging


# Return if a read is clipped on the left
def is_left_clipped(read):
    if read.cigartuples[0][0] in [4, 5]:
        return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    if read.cigartuples[-1][0] in [4, 5]:
        return True
    return False


def get_clipped_reads(ibam, chrName, outFile):

    clipped_reads = dict()
    for split_direction in ['left', 'right']:
        clipped_reads[split_direction] = dict()

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    # print([i['LN'] for i in header_dict['SQ'] if i['SN']])
    # print([i['SN'] for i in header_dict['SQ'] if i['SN']])
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    # print(chrLen)

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    i = 0
    n_r = 10**6
    #print(n_r)
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

        if not read.is_unmapped and not read.mate_is_unmapped:
            if is_left_clipped(read):
                # print(str(read))
                # print('Clipped at the start: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                # print('Pos:%d, clipped_pos:%d' % (read.reference_start, read.get_reference_positions()[0]))
                # print('start:'+str(read.get_reference_positions()[0])+'=='+str(read.reference_start))
                ref_pos = read.reference_start
                if ref_pos not in clipped_reads['left'].keys():
                    clipped_reads['left'][ref_pos] = 1
                else:
                    clipped_reads['left'][ref_pos] += 1
            if is_right_clipped(read):
                # print('Clipped at the end: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                # print('Pos:%d, clipped_pos:%d' %(read.reference_end, read.get_reference_positions()[-1]))
                # print('end: '+str(read.get_reference_positions()[-1]) + '==' + str(read.reference_end))
                ref_pos = read.reference_end + 1
                if ref_pos not in clipped_reads['right'].keys():
                    clipped_reads['right'][ref_pos] = 1
                else:
                    clipped_reads['right'][ref_pos] += 1

    # cPickle data persistence
    with bz2file.BZ2File(outFile, 'wb') as f:
        pickle.dump(clipped_reads, f)


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
    parser.add_argument('-l', '--logfile', default='clipped_reads.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    get_clipped_reads(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
