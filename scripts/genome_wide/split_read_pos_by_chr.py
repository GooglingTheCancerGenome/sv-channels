# Imports
import argparse
import logging
import os
import json
import gzip
from collections import Counter, defaultdict
from time import time
import pysam
from functions import *

strand_str = {True: '-', False: '+'}


def get_split_read_positions(ibam, chrName, outFile):
    '''
    
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the dictionary of split read positions
    :return: None. Outputs a dictionary with the positions of split read positions as keys and
    the number of split reads per position as values
    '''''

    def overlap(read1, read2):
        if read1.reference_start <= read2.reference_start <= read1.reference_end or \
                read1.reference_start <= read2.reference_end <= read1.reference_end:
            return True
        else:
            return False

    def get_distance(chrA, refpos, chrB, pos, chrlen_dict, chrlist):

        dist = 0
        if chrA == chrB:
            dist = abs(refpos - pos)
        else:
            if chrName < chr:
                chr1 = chrA
                chr2 = chrB
            else:
                chr1 = chrB
                chr2 = chrA
            dist += chrlen_dict[chr1] - refpos
            i = chrlist.index(chr1) + 1
            while chrlist[i] != chr2:
                dist += chrlen_dict[chrlist[i]]
                i += 1
            dist += pos

        return dist

    def append_coord(split_pos_coord, chrName, refpos, chr, pos):

        if chrName == chr:
            if refpos < pos:
                split_pos_coord.append((chrName, refpos, chr, pos))
            else:
                split_pos_coord.append((chr, pos, chrName, refpos))
        elif chrName < chr:
            split_pos_coord.append((chrName, refpos, chr, pos))
        elif chrName > chr:
            split_pos_coord.append((chr, pos, chrName, refpos))

        return split_pos_coord

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    config = get_config_file()
    minMAPQ = config["DEFAULT"]["MIN_MAPQ"]
    min_support = config["DEFAULT"]["MIN_SR_SUPPORT"]

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    chrlen_dict = {i['SN']: i['LN'] for i in header_dict['SQ']}
    chrlist = sorted(chrlen_dict.keys())

    # List to store the split read positions
    split_pos = []
    split_pos_coord = []

    discordant_reads_pos = []
    discordant_reads_coord = []

    reads_in_cluster = defaultdict()

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # Pysam iterator to fetch the reads
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    # Print every n_r alignments processed
    n_r = 10 ** 6
    # Record the current time
    last_t = time()

    n_indels = 0
    n_split = 0
    n_discordant = 0
    max_cigar_del = 0

    for i, read in enumerate(iter, start=1):

        # Every n_r alignments, write log informations
        if not i % n_r:
            # Record the current time
            now_t = time()
            # print(type(now_t))
            logging.info("%d alignments processed (%f alignments / s)" % (
                i,
                n_r / (now_t - last_t)))
            last_t = time()

        # Both read and mate should be mapped, read should have a minimum mapping quality
        # if (not read.is_unmapped) and (not read.mate_is_unmapped) and read.mapping_quality >= minMAPQ:
        if (not read.is_unmapped) and read.mapping_quality >= minMAPQ:

            if has_indels(read):
                # print(read)

                dels_start, dels_end, ins = get_indels(read)
                dels = dels_start + dels_end + ins
                split_pos.extend(dels)

                for start, end in zip(dels_start, dels_end):

                    n_indels += 1

                    # Calculate DEL size and find largest DEL size encoded by CIGAR 'D' character
                    del_size = end - start
                    if del_size > max_cigar_del:
                        max_cigar_del = del_size

                    split_pos_coord = append_coord(split_pos_coord, chrName, start, chrName, end)
                for pos in ins:
                    n_indels += 1
                    split_pos_coord = append_coord(split_pos_coord, chrName, pos, chrName, pos+1)

            if read.has_tag('SA'):
                n_split += 1
                chr, pos, strand = get_suppl_aln(read)
                if is_right_clipped(read):
                    refpos = read.reference_end
                    if chr == read.reference_name:
                        split_pos.append(refpos)
                        split_pos.append(pos)
                        # print('{}:{}-{}'.format(chrName, refpos, pos))
                        split_pos_coord = append_coord(split_pos_coord, chrName, refpos, chr, pos)
                    else:
                        split_pos.append(refpos)
                elif is_left_clipped(read):
                    refpos = read.reference_start
                    if chr == read.reference_name:
                        pass
                    else:
                        split_pos.append(refpos)
                        # print('{}:{}-{}'.format(chrName, refpos, pos))
                        split_pos_coord = append_coord(split_pos_coord, chrName, refpos, chr, pos)

            # if not read.mate_is_unmapped and ( not read.is_proper_pair or is_clipped(read) ):
            #     n_discordant += 1
            #     refpos = read.reference_end + 1 if not read.is_reverse else read.reference_start
            #
            #     strand_id = '_'.join([read.reference_name, strand_str[read.is_reverse],
            #                           read.next_reference_name, strand_str[read.mate_is_reverse]])
            #
            #     if strand_id in reads_in_cluster.keys():
            #
            #         chr1, pos1, chr2, pos2, cnt, read_ref = reads_in_cluster[strand_id]
            #
            #         if overlap(read, read_ref):
            #
            #             if (not read.is_reverse and refpos > pos1) or \
            #                     (read.is_reverse and refpos < pos1):
            #                 pos1 = refpos
            #
            #             if (read.next_reference_start > pos2 and not read.mate_is_reverse) or \
            #                     (read.next_reference_start < pos2 and read.mate_is_reverse):
            #                 pos2 = read.next_reference_start
            #
            #             cnt += 1
            #             # print('Adding %s_%d_%s_%d_%d' % (chr1, pos1, chr2, pos2, cnt))
            #             reads_in_cluster[strand_id] = (chr1, pos1, chr2, pos2, cnt, read)
            #
            #     else:
            #         reads_in_cluster[strand_id] = (read.reference_name, refpos,
            #                                        read.next_reference_name, read.next_reference_start, 1, read)
            # else:
            #
            #     for strand_id in reads_in_cluster.keys():
            #         chr1, pos1, chr2, pos2, cnt, read = reads_in_cluster[strand_id]
            #         if cnt >= min_support:
            #             discordant_reads_pos.extend([pos1]*cnt)
            #             discordant_reads_coord = append_coord(discordant_reads_coord, chr1, pos1, chr2, pos2)
            #
            #     reads_in_cluster = defaultdict()

    # Close the BAM file
    bamfile.close()

    # Count the number of split reads per position
    split_pos_cnt = Counter(split_pos)
    split_pos_coord = set(split_pos_coord)

    logging.info('Largest CIGAR "D" DEL={}'.format(max_cigar_del))

    logging.info('INDELs={}, split_reads={}, discordant_reads={}'.format(n_indels,
                                                                         n_split,
                                                                         n_discordant))

    logging.info('Number of unique split read positions: %d' % len(
        [p for p, c in split_pos_cnt.items() if c >= min_support])
                 )
    logging.info('Number of unique pair of split read positions: %d' % len(split_pos_coord))

    discordant_reads_cnt = Counter(discordant_reads_pos)
    discordant_reads_coord = set(discordant_reads_coord)

    logging.info('Number of unique discordant positions: %d' % len(
        [p for p, c in discordant_reads_cnt.items() if c >= min_support])
                 )
    logging.info('Number of unique pair of discordant positions: %d' % len(discordant_reads_coord))

    total_reads_cnt = Counter(split_pos + discordant_reads_pos)

    total_reads_coord = list(set(split_pos_coord | discordant_reads_coord))

    positions_with_min_support = [p for p, c in total_reads_cnt.items() if c >= min_support]
    logging.info('Number of positions with min %d support: %d' % (min_support, len(positions_with_min_support)))
    logging.info('Number of unique pair of total positions: %d' % len(total_reads_coord))

    positions_with_min_support_set = set(positions_with_min_support)
    total_reads_coord_min_support = [(chr1, pos1, chr2, pos2) for chr1, pos1, chr2, pos2
                                     in total_reads_coord
                                     if pos1 in positions_with_min_support_set or
                                     pos2 in positions_with_min_support_set]

    logging.info('Number of total pairs of positions with min support: %d' % len(total_reads_coord_min_support))

    data = (positions_with_min_support, total_reads_coord_min_support)
    # Write
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))

    # to load it:
    # with gzip.GzipFile(outFile, 'r') as fin:
    #     total_reads_cnt, total_reads_coord = json.loads(fin.read().decode('utf-8'))


def main():

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"
    # wd = '/Users/lsantuari/Documents/mount_points/hpc_mnt/Datasets/CretuStancu2017/Patient1/'
    # inputBAM = wd + 'Patient1.bam'

    # Default chromosome is 17 for the artificial data

    # Parse the arguments of the script
    parser = argparse.ArgumentParser(description='Get split reads positions')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='split_read_pos.json.gz',
                        help="Specify output")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-l', '--logfile', default='split_read_pos.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    cmd_name = 'split_read_pos'
    output_dir = os.path.join(args.outputpath, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, '_'.join((args.chr, args.logfile)))
    output_file = os.path.join(output_dir, '_'.join((args.chr, args.out)))

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    get_split_read_positions(ibam=args.bam, chrName=args.chr, outFile=output_file)
    logging.info('Time: split read positions on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
