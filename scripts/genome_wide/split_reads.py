# Imports
import argparse
import logging
import os
import json
import gzip
from collections import Counter, defaultdict
from time import time
import pysam
from intervaltree import IntervalTree
from functions import *

strand_str = {True: '-', False: '+'}


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


def get_split_read_positions(ibam, chr_list, outFile):
    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    config = get_config_file()
    minMAPQ = config["DEFAULT"]["MIN_MAPQ"]
    min_support = config["DEFAULT"]["MIN_SR_SUPPORT"]

    # List to store the split read positions
    split_pos_coord = []

    # Tree with windows for candidate positions
    gtrees = defaultdict(IntervalTree)

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")

    chrlen_dict = get_chr_len_dict(ibam)

    n_indels = 0
    n_split = 0
    n_discordant = 0
    max_cigar_del = 0

    split_reads = dict()
    split_read_distance = dict()
    for chrom in chr_list:
        split_reads[chrom] = dict()
        split_read_distance[chrom] = dict()
        for split_direction in ['left_F', 'left_R', 'right_F', 'right_R']:
            split_reads[chrom][split_direction] = defaultdict(int)
            split_read_distance[chrom][split_direction] = defaultdict(list)

    right_split_pos = defaultdict(list, {k: [] for k in chr_list})
    left_split_pos = defaultdict(list, {k: [] for k in chr_list})

    # Pysam iterator to fetch the reads
    iter = bamfile.fetch()

    # Print every n_r alignments processed
    n_r = 10 ** 6
    # Record the current time
    last_t = time()

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

        if not read.is_unmapped and read.mapping_quality >= minMAPQ and \
                read.reference_name in chr_list:

            if has_indels(read):
                # print(read)

                dels_start, dels_end, ins = get_indels(read)
                # dels = dels_start + dels_end + ins
                left_split_pos[read.reference_name].extend(dels_end)
                right_split_pos[read.reference_name].extend(dels_start)

                for start, end in zip(dels_start, dels_end):

                    n_indels += 1

                    # Calculate DEL size and find largest DEL size encoded by CIGAR 'D' character
                    del_size = end - start
                    if del_size > max_cigar_del:
                        max_cigar_del = del_size

                    split_pos_coord = append_coord(split_pos_coord,
                                                   read.reference_name, start,
                                                   read.reference_name, end)

            if read.has_tag('SA'):

                chr_SA, pos_SA, strand_SA = get_suppl_aln(read)
                if is_left_clipped(read) and is_right_clipped(read):
                    # print('Read clipped on both sides:\n{}'.format(read))
                    continue
                elif is_right_clipped(read) and read.reference_name == chr_SA and \
                        read.reference_end < pos_SA and strand_str[read.is_reverse] == strand_SA:

                    n_split += 1

                    # print('split {} at position {}:{} to {}:{}'.format(
                    #     read.query_name,
                    #     read.reference_name,
                    #     read.reference_end,
                    #     chr_SA,
                    #     pos_SA
                    #     ))

                    split_pos_coord = append_coord(split_pos_coord,
                                                   read.reference_name,
                                                   read.reference_end,
                                                   chr_SA,
                                                   pos_SA)

                    if not read.is_reverse:
                        split_reads[read.reference_name]['right_F'][read.reference_end] += 1
                        split_reads[read.reference_name]['left_F'][pos_SA] += 1
                    else:
                        split_reads[read.reference_name]['right_R'][read.reference_end] += 1
                        split_reads[read.reference_name]['left_R'][pos_SA] += 1

                    right_split_pos[read.reference_name].append(read.reference_end)
                    left_split_pos[read.reference_name].append(pos_SA)

                    dist = abs(read.reference_end - pos_SA)
                    if not read.is_reverse:
                        split_read_distance[read.reference_name]['right_F'][read.reference_end].append(dist)
                        split_read_distance[read.reference_name]['left_F'][pos_SA].append(dist)
                    else:
                        split_read_distance[read.reference_name]['right_R'][read.reference_end].append(dist)
                        split_read_distance[read.reference_name]['left_R'][pos_SA].append(dist)

    # Close the BAM file
    bamfile.close()

    right_split_pos_cnt = dict.fromkeys(chr_list)
    left_split_pos_cnt = dict.fromkeys(chr_list)

    # Count the number of split reads per position
    for chrom in chr_list:
        right_split_pos_cnt[chrom] = Counter(right_split_pos[chrom])
        left_split_pos_cnt[chrom] = Counter(left_split_pos[chrom])

    split_pos_coord = set(split_pos_coord)

    logging.info('Largest CIGAR "D" DEL={}'.format(max_cigar_del))

    logging.info('INDELs={}, split_reads={}, discordant_reads={}'.format(n_indels,
                                                                         n_split,
                                                                         n_discordant))

    for chrom in chr_list:
        logging.info('Number of unique left split read positions on Chr{}: {}'.format(chrom,
                                                                                      len(
                                                                                          [p for p, c in
                                                                                           left_split_pos_cnt[
                                                                                               chrom].items() if
                                                                                           c >= min_support]
                                                                                      )
                                                                                      ))
        logging.info('Number of unique right split read positions on Chr{}: {}'.format(chrom,
                                                                                       len(
                                                                                           [p for p, c in
                                                                                            right_split_pos_cnt[
                                                                                                chrom].items() if
                                                                                            c >= min_support]
                                                                                       )
                                                                                       ))

    logging.info('Number of unique pair of split read positions: %d' % len(split_pos_coord))

    # discordant_reads_cnt = Counter(discordant_reads_pos)
    # discordant_reads_coord = set(discordant_reads_coord)
    #
    # logging.info('Number of unique discordant positions: %d' % len(
    #     [p for p, c in discordant_reads_cnt.items() if c >= min_support])
    #             )
    # logging.info('Number of unique pair of discordant positions: %d' % len(discordant_reads_coord))

    total_reads_cnt_ls = dict.fromkeys(chr_list)
    total_reads_cnt_rs = dict.fromkeys(chr_list)

    for chrom in chr_list:
        total_reads_cnt_ls[chrom] = Counter(left_split_pos[chrom])
        total_reads_cnt_rs[chrom] = Counter(right_split_pos[chrom])

    total_reads_coord = list(set(split_pos_coord))  # | discordant_reads_coord))

    positions_with_min_support_ls = dict.fromkeys(chr_list)
    positions_with_min_support_rs = dict.fromkeys(chr_list)

    for chrom in chr_list:
        positions_with_min_support_ls[chrom] = [p for p, c in total_reads_cnt_ls[chrom].items() if c >= min_support]
        positions_with_min_support_rs[chrom] = [p for p, c in total_reads_cnt_rs[chrom].items() if c >= min_support]

        logging.info('Number of LS positions on Chr%s with min %d support: %d' % (chrom, min_support,
                                                                                  len(positions_with_min_support_ls[
                                                                                          chrom])))
        logging.info('Number of RS positions on Chr%s with min %d support: %d' % (chrom, min_support,
                                                                                  len(positions_with_min_support_rs[
                                                                                          chrom])))

    logging.info('Number of unique pair of total positions: %d' % len(total_reads_coord))

    positions_with_min_support_set = dict.fromkeys(chr_list)

    for chrom in chr_list:
        positions_with_min_support_set[chrom] = set(positions_with_min_support_ls[chrom] +
                                                    positions_with_min_support_rs[chrom])

    total_reads_coord_min_support = [(chr1, pos1, chr2, pos2) for chr1, pos1, chr2, pos2
                                     in total_reads_coord
                                     if pos1 in positions_with_min_support_set[chr1] or
                                     pos2 in positions_with_min_support_set[chr2]]

    logging.info('Number of total pairs of positions with min support: %d' % len(total_reads_coord_min_support))

    data = (positions_with_min_support_ls, positions_with_min_support_rs, total_reads_coord_min_support,
            split_reads, split_read_distance)
    # Write
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))


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
    # inputBAM = "/Users/lsantuari/Documents/mount_points/hpc_giab/RMNISTHS_30xdownsample.bam"

    # Default chromosome is 17 for the artificial data

    # Parse the arguments of the script
    parser = argparse.ArgumentParser(description='Get split reads positions')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chrlist', nargs='+', default=['17'],
                        help="List of chromosomes to consider")
    parser.add_argument('-o', '--out', type=str, default='split_reads.json.gz',
                        help="Specify output")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-l', '--logfile', default='split_reads.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    cmd_name = 'split_reads'
    output_dir = os.path.join(args.outputpath, cmd_name)
    create_dir(output_dir)
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, args.out)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    get_split_read_positions(ibam=args.bam, chr_list=args.chrlist, outFile=output_file)
    logging.info('Time: split read positions on BAM %s: %f' % (args.bam, (time() - t0)))


if __name__ == '__main__':
    main()
