import argparse
import pysam
import bz2file
import pickle
from time import time
import functions as fun
import logging


def get_clipped_read_distance(ibam, chrName, clipped_reads, outFile):

    '''
    # Load clipped reads sets
    logging.info("Loading %s" % clipped_reads)
    with bz2file.BZ2File(clipped_reads, 'rb') as f:
        clipped_pos_cnt, clipped_read_1, clipped_read_2 = pickle.load(f)
    logging.info("Loaded")
    '''

    clipped_read_distance = dict()
    for direction in ['forward', 'reverse']:
        clipped_read_distance[direction] = dict()
        for clipped_arrangement in ['left', 'right']:
            clipped_read_distance[direction][clipped_arrangement] = dict()

    #for direction in ['forward', 'reverse']:
    #    clipped_read_distance[direction] = dict()
    #for direction in ['forward', 'reverse']:
    #    for clipped_arrangement in ['c2c', 'nc2c', 'c2nc', 'nc2nc']:
    #        clipped_read_distance[direction][clipped_arrangement] = dict()

    def get_distance(direction, read, dist):

        if fun.is_left_clipped(read):
            pos = read.reference_start
            if pos not in clipped_read_distance[direction]['left'].keys():
                clipped_read_distance[direction]['left'][pos] = [dist]
            else:
                clipped_read_distance[direction]['left'][pos].append(dist)
        elif fun.is_right_clipped(read):
            pos = read.reference_end + 1
            if pos not in clipped_read_distance[direction]['right'].keys():
                clipped_read_distance[direction]['right'][pos] = [dist]
            else:
                clipped_read_distance[direction]['right'][pos].append(dist)

        '''
        mate_is_clipped = read.is_read1 and read.query_name in clipped_read_2 or \
                          read.is_read2 and read.query_name in clipped_read_1

        if fun.is_clipped(read) and mate_is_clipped:
            if fun.is_left_clipped(read):
                pos = read.get_reference_positions()[0]
                if pos not in clipped_read_distance[direction]['c2c'].keys():
                    clipped_read_distance[direction]['c2c'][pos] = [d]
                else:
                    clipped_read_distance[direction]['c2c'][pos].append(d)
            elif fun.is_right_clipped(read):
                pos = read.get_reference_positions()[-1]
                if pos not in clipped_read_distance[direction]['c2c'].keys():
                    clipped_read_distance[direction]['c2c'][pos] = [d]
                else:
                    clipped_read_distance[direction]['c2c'][pos].append(d)

        elif fun.is_clipped(read) and not mate_is_clipped:
            if fun.is_left_clipped(read):
                pos = read.get_reference_positions()[0]
                if pos not in clipped_read_distance[direction]['c2nc'].keys():
                    clipped_read_distance[direction]['c2nc'][pos] = [d]
                else:
                    clipped_read_distance[direction]['c2nc'][pos].append(d)
            elif fun.is_right_clipped(read):
                pos = read.get_reference_positions()[-1]
                if pos not in clipped_read_distance[direction]['c2nc'].keys():
                    clipped_read_distance[direction]['c2nc'][pos] = [d]
                else:
                    clipped_read_distance[direction]['c2nc'][pos].append(d)
        elif not fun.is_clipped(read) and mate_is_clipped:
            if read.reference_start not in clipped_read_distance[direction]['nc2c'].keys():
                clipped_read_distance[direction]['nc2c'][read.reference_start] = [d]
            else:
                clipped_read_distance[direction]['nc2c'][read.reference_start].append(d)
        elif not fun.is_clipped(read) and not mate_is_clipped:
            if read.reference_start not in clipped_read_distance[direction]['nc2nc'].keys():
                if read.reference_start not in clipped_read_distance[direction]['nc2nc'].keys():
                    clipped_read_distance[direction]['nc2nc'][read.reference_start] = [d]
                else:
                    clipped_read_distance[direction]['nc2nc'][read.reference_start].append(d)
        '''

    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header

    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    start_pos = 0
    stop_pos = chrLen
    # print(chrLen)

    iter = bamfile.fetch(chrName, start_pos, stop_pos)

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

        if not read.is_unmapped and not read.mate_is_unmapped and len(read.get_reference_positions()) > 0:
            # mate = fun.get_read_mate(read, bamfile)
            if read.reference_name == read.next_reference_name:
                dist = abs(read.reference_start - read.next_reference_start)

                if not read.is_reverse and read.mate_is_reverse and read.reference_start <= read.next_reference_start:
                    # pass
                    get_distance('forward', read, dist)
                elif read.is_reverse and not read.mate_is_reverse and read.reference_start >= read.next_reference_start:
                    # pass
                    get_distance('reverse', read, dist)

    # print(clipped_read_distance)
    # cPickle data persistence
    with bz2file.BZ2File(outFile, 'w') as f:
        pickle.dump(clipped_read_distance, f)


def main():
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    parser = argparse.ArgumentParser(description='Create channels with distance between clipped/non-clipped reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-r', '--reads', type=str, default='clipped_read_pos_cr.pbz2',
                        help="Specify clipped read position file")
    parser.add_argument('-o', '--out', type=str, default='clipped_read_distance.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='clipped_read_distance.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    get_clipped_read_distance(ibam=args.bam, chrName=args.chr, clipped_reads=args.reads, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()
