import argparse
import pysam
import bz2
import cPickle as pickle
from time import time


# Return if a read is clipped on the left
def is_left_clipped(read):
    if not read.is_unmapped and not read.mate_is_unmapped:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    if not read.is_unmapped and not read.mate_is_unmapped:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False

# Return if a read is clipped on the right
def deletion(read):
    if not read.is_unmapped and not read.mate_is_unmapped and read.is_reverse != read.mate_is_reverse:
        if read.reference_name == read.next_reference_name:
            if abs(read.reference_start - read.next_reference_start) >= 300:
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
    stop_pos = 1000000
    # print(chrLen)

    A, C, G, T = bamfile.count_coverage(chrName, start_pos, stop_pos, quality_threshold=20, read_callback=is_left_clipped)
    v = map(sum, zip(A, C, G, T))
    #d = dict()
    #for i in range(len(v)-1):
    #    if v[i] <= v[i+1]:
    #        v[i] = 0L
    #print(len(v))

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    for read in iter:
        if not read.is_unmapped and not read.mate_is_unmapped:
            if is_left_clipped(read):
                # print(str(read))
                # print('Clipped at the start: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                # print('Pos:%d, clipped_pos:%d' % (read.reference_start, read.get_reference_positions()[0]))
                # print('start:'+str(read.get_reference_positions()[0])+'=='+str(read.reference_start))
                ref_pos = read.get_reference_positions()[0]
                if ref_pos not in clipped_reads['left'].keys():
                    clipped_reads['left'][ref_pos] = 1
                else:
                    clipped_reads['left'][ref_pos] += 1

    for k in clipped_reads['left'].keys():
        if clipped_reads['left'][k] != v[k]:
            print(str(clipped_reads['left'][k]) + '!=' + str(v[k]) + ' at ' + str(k) )

    # cPickle data persistence
    with bz2.BZ2File(outFile, 'w') as f:
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

    args = parser.parse_args()

    t0 = time()
    get_clipped_reads(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
