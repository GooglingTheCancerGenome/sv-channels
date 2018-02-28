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

# Return chromosome and starting position of a supplementary alignment (split reads)
def get_suppl_aln(read):

    if len(read.get_tag('SA')) > 0:
        #print(read.get_tag('SA'))
        supp_aln = read.get_tag('SA').split(';')[0]
        sa_info = supp_aln.split(',')
        #print(supp_aln)
        #print(sa_info)
        chr_sa = sa_info[0]
        start_sa = int(sa_info[1])
        return chr_sa, start_sa
    else:
        return None

def get_split_read_distance(ibam, chrName, outFile):

    split_read_distance = dict()
    for split_direction in ['left', 'right']:
        split_read_distance[split_direction] = dict()

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    # print(chrLen)

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    for read in iter:
        if not read.is_unmapped and not read.mate_is_unmapped:  # and read.reference_start >= start:
            if is_left_clipped(read) and read.has_tag('SA'):
                chr, pos = get_suppl_aln(read)
                if chr == read.reference_name:
                    # print('Left split')
                    # print(str(read))
                    refpos = read.get_reference_positions()[0]
                    if pos not in split_read_distance['left'].keys():
                        split_read_distance['left'][pos] = [abs(refpos - pos)]
                    else:
                        split_read_distance['left'][pos].append(abs(refpos - pos))
            elif is_right_clipped(read) and read.has_tag('SA'):
                chr, pos = get_suppl_aln(read)
                if chr == read.reference_name:
                    # print('Right split')
                    # print(str(read))
                    refpos = read.get_reference_positions()[-1]
                    if pos not in split_read_distance['right'].keys():
                        split_read_distance['right'][pos] = [abs(pos - refpos)]
                    else:
                        split_read_distance['right'][pos].append(abs(pos - refpos))

    #print(split_read_distance)
    # cPickle data persistence
    with bz2.BZ2File(outFile, 'w') as f:
        pickle.dump(split_read_distance, f)


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    parser = argparse.ArgumentParser(description='Create channels with split read distance for left/right split reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='split_read_distance.pbz2',
                        help="Specify output")

    args = parser.parse_args()

    t0 = time()
    get_split_read_distance(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()