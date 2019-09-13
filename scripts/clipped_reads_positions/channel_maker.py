import pysam
import numpy as np
import argparse
import gzip
import os
import logging
from time import time


def initialize_variables():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    bam1 = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"
    bam2 = wd + "reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/" + "G1_dedup.bam"

    return bam1, bam2


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


# Return chromosome and starting position of a supplementary alignment (split reads)
def get_suppl_aln(read):
    if len(read.get_tag('SA')) > 0:
        # print(read.get_tag('SA'))
        supp_aln = read.get_tag('SA').split(';')[0]
        sa_info = supp_aln.split(',')
        # print(supp_aln)
        # print(sa_info)
        chr_sa = sa_info[0]
        start_sa = int(sa_info[1])
        return chr_sa, start_sa
    else:
        return None


def load_bam(bamFile):
    assert os.path.isfile(bamFile)
    return pysam.AlignmentFile(bamFile, "rb")


def get_chr_len(chrName, bamFile):
    return [i['LN'] for i in bamFile.header['SQ'] if i['SN'] == chrName][0]


def get_coverage(start_win, end_win, chrName, bamFile):
    cov = np.zeros(abs(end_win - start_win), dtype=int)
    for pile in bamFile.pileup(chrName, start_win, end_win, truncate=True):
        cov[pile.pos-start_win] = pile.n
    return cov


def create_channel_matrix(centerPos, chrName, bamFile1, bamFile2):

    def fill_vectors(iter, sample, start_win, end_win, clipped_reads,
                     clipped_reads_distance, split_reads, split_reads_distance):

        for read in iter:
            if not read.is_unmapped and not read.mate_is_unmapped:
                if is_left_clipped(read):
                    pos = read.reference_start - start_win
                    clipped_reads[sample]['left'][pos] += 1
                    if chrName == read.next_reference_name and read.reference_start <= read.next_reference_start:
                        clipped_reads_distance[sample]['left'][pos] += abs(pos - read.next_reference_start)
                    if read.has_tag('SA'):
                        split_reads[sample]['left'][pos] += 1
                        chr_sa, start_sa = get_suppl_aln(read)
                        if chr_sa == chrName:
                            split_reads_distance[sample]['left'][pos] += abs(pos - start_sa)
                if is_right_clipped(read):
                    pos = read.reference_end - start_win
                    if pos <= abs(end_win - end_win)-1:
                        clipped_reads[sample]['right'][pos] += 1
                        if chrName == read.next_reference_name and read.reference_start <= read.next_reference_start:
                            clipped_reads_distance[sample]['right'][pos] += abs(pos - read.next_reference_start)
                        if read.has_tag('SA'):
                            split_reads[sample]['right'][pos] += 1
                            chr_sa, start_sa = get_suppl_aln(read)
                            if chr_sa == chrName:
                                split_reads_distance[sample]['right'][pos] += abs(pos - start_sa)

    win_len = 100
    start_win = centerPos - win_len
    end_win = centerPos + win_len

    iter1 = bamFile1.fetch(chrName, start_win, end_win, multiple_iterators=True)
    iter2 = bamFile2.fetch(chrName, start_win, end_win)

    clipped_reads = dict()
    clipped_reads_distance = dict()
    split_reads = dict()
    split_reads_distance = dict()

    sampleList = ['bam1', 'bam2']

    for sample in sampleList:
        clipped_reads[sample] = dict()
        clipped_reads_distance[sample] = dict()
        split_reads[sample] = dict()
        split_reads_distance[sample] = dict()

        for direction in ['left', 'right']:
            clipped_reads[sample][direction] = np.zeros(abs(end_win - start_win), dtype=int)
            split_reads[sample][direction] = np.zeros(abs(end_win - start_win), dtype=int)
            split_reads_distance[sample][direction] = np.zeros(abs(end_win - start_win), dtype=int)
            clipped_reads_distance[sample][direction] = np.zeros(abs(end_win - start_win), dtype=int)
            #print(clipped_reads_distance[sample][direction])

    for sample, iter in zip(sampleList, [iter1, iter2]):
        #print(sample)
        fill_vectors(iter, sample, start_win, end_win, clipped_reads,
                     clipped_reads_distance, split_reads, split_reads_distance)

    ch_vstack = np.vstack((
        get_coverage(start_win, end_win, chrName, bamFile1),
        clipped_reads['bam1']['left'],
        clipped_reads['bam1']['right'],
        clipped_reads_distance['bam1']['left'],
        clipped_reads_distance['bam1']['right'],
        split_reads['bam1']['left'],
        split_reads['bam1']['right'],
        split_reads_distance['bam1']['left'],
        split_reads_distance['bam1']['right'],
        get_coverage(start_win, end_win, chrName, bamFile2),
        clipped_reads['bam2']['left'],
        clipped_reads['bam2']['right'],
        clipped_reads_distance['bam2']['left'],
        clipped_reads_distance['bam2']['right'],
        split_reads['bam2']['left'],
        split_reads['bam2']['right'],
        split_reads_distance['bam2']['left'],
        split_reads_distance['bam2']['right']
    ))

    #print(ch_vstack)
    return ch_vstack


def channel_maker(ibam1, ibam2, chrName, outFile, labelFile):
    cr_support = 3

    bam1_file = load_bam(ibam1)
    bam2_file = load_bam(ibam2)

    start_pos = 0
    stop_pos = get_chr_len(chrName, bam1_file)

    iter1 = bam1_file.fetch(chrName, start_pos, stop_pos, multiple_iterators=True)

    n_r = 10 ** 6
    last_t = time()

    clipped_positions = dict()
    window_positions = set()
    channel_list = []

    for i, read in enumerate(iter1, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        if not read.is_unmapped and not read.mate_is_unmapped:
            if is_left_clipped(read):
                pos = read.reference_start
                if pos not in clipped_positions.keys():
                    clipped_positions[pos] = 1
                else:
                    clipped_positions[pos] += 1
                if pos not in window_positions and clipped_positions[pos] >= cr_support:
                    #print(read)
                    channel_list.append(create_channel_matrix(pos, chrName, bam1_file, bam2_file))
                    window_positions.add(pos)
            if is_right_clipped(read):
                pos = read.reference_end + 1
                if pos not in clipped_positions.keys():
                    clipped_positions[pos] = 1
                else:
                    clipped_positions[pos] += 1
                if pos not in window_positions and clipped_positions[pos] >= cr_support:
                    #print(read)
                    channel_list.append(create_channel_matrix(pos, chrName, bam1_file, bam2_file))
                    window_positions.add(pos)

    logging.info("Number of windows: %d" % len(channel_list))
    with gzip.GzipFile(outFile, "w") as f:
        np.save(file=f, arr=channel_list)
    f.close()


def main():

    bam1, bam2 = initialize_variables()

    parser = argparse.ArgumentParser(description='Create channels from BAM files')
    parser.add_argument('-b1', '--bam1', type=str,
                        default=bam1,
                        help="Specify first input file (BAM)")
    parser.add_argument('-b2', '--bam2', type=str,
                        default=bam2,
                        help="Specify second input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='channel_maker.npy.gz',
                        help="Specify output for channels")
    parser.add_argument('-ol', '--labels', type=str, default='channel_maker.npy.gz',
                        help="Specify output for labels")
    parser.add_argument('-l', '--logfile', default='channel_maker.log',
                        help='Specify log file')

    args = parser.parse_args()

    log_filename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=log_filename,
        level=logging.INFO)

    t0 = time()
    channel_maker(ibam1=args.bam1, ibam2=args.bam2, chrName=args.chr, outFile=args.out, labelFile=args.labels)
    print(time() - t0)


if __name__ == '__main__':
    main()
