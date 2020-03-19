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
from cigar import Cigar

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


def get_split_read_positions(ibam, chr_list, outFile, outBedpe):

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    config = get_config_file()
    minMAPQ = config["DEFAULT"]["MIN_MAPQ"]
    min_support = config["DEFAULT"]["MIN_SR_SUPPORT"]

    # List to store the split read positions
    split_pos_coord = dict()
    for k in ['INDEL_INS', 'INDEL_DEL', 'DEL', 'INS', 'INV', 'DUP', 'TRA']:
        split_pos_coord[k] = []

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

    clipped_pos_dict = dict()
    split_pos = dict()
    split_pos_cnt = dict()
    total_reads_cnt = dict()
    positions_with_min_support = dict()

    for k in ['right', 'left']:
        clipped_pos_dict[k] = defaultdict(list, {k: [] for k in chr_list})
        split_pos[k] = defaultdict(list, {k: [] for k in chr_list})
        split_pos_cnt[k] = dict.fromkeys(chr_list)
        total_reads_cnt[k] = dict.fromkeys(chr_list)
        positions_with_min_support[k] = dict.fromkeys(chr_list)

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

        if is_left_clipped(read):
            clipped_pos_dict['left'][read.reference_name].append(read.reference_start + 1)
        elif is_right_clipped(read):
            clipped_pos_dict['right'][read.reference_name].append(read.reference_end)

        if not read.is_unmapped and read.mapping_quality >= minMAPQ and \
                read.reference_name in chr_list:

            if has_indels(read):
                # print(read)

                dels_start, dels_end, ins = get_indels(read)
                # dels = dels_start + dels_end + ins
                split_pos['left'][read.reference_name].extend(dels_end)
                split_pos['right'][read.reference_name].extend(dels_start)

                for pos in ins:
                    n_indels += 1

                    split_pos_coord['INDEL_INS'] = append_coord(split_pos_coord['INDEL_INS'],
                                                                read.reference_name, start,
                                                                read.reference_name, end)

                for start, end in zip(dels_start, dels_end):

                    n_indels += 1

                    # Calculate DEL size and find largest DEL size encoded by CIGAR 'D' character
                    del_size = end - start
                    if del_size > max_cigar_del:
                        max_cigar_del = del_size

                    split_pos_coord['INDEL_DEL'] = append_coord(split_pos_coord['INDEL_DEL'],
                                                                read.reference_name, start,
                                                                read.reference_name, end)

            if read.has_tag('SA') and is_left_clipped(read) != is_right_clipped(read):

                sa_entry = get_suppl_aln(read)

                if sa_entry is not None:

                    n_split += 1

                    chr_SA, pos_SA, strand_SA, cigar_sa = sa_entry
                    cigar_sa = Cigar(cigar_sa)
                    cigar_sa_list = list(cigar_sa.items())

                    if cigar_sa_list[-1][1] in ['S', 'H']:
                        # print('SA_RS')
                        # print(read.query_name)
                        # print('{} {} {} {}'.format(chr_SA, pos_SA, strand_SA, cigar_sa))

                        cigar_sa_len = 0
                        # print(cigar_sa_list)
                        for l, c in cigar_sa_list:
                            # query consuming CIGAR characters
                            if c in ['M', 'D', 'N', '=', 'X']:
                                cigar_sa_len += l
                        pos_SA += cigar_sa_len

                        # print('{} {} {} {}'.format(chr_SA, pos_SA, strand_SA, cigar_sa))
                        # print('---')

                    elif cigar_sa_list[0][1] in ['S', 'H']:
                        pos_SA += 1

                    if chr_SA in chr_list:

                        clipped_string = 'right' if is_right_clipped(read) else 'left'
                        clipped_orient = 'F' if not read.is_reverse else 'R'
                        clipped_ch = '_'.join([clipped_string, clipped_orient])

                        clipped_pos = read.reference_end if is_right_clipped(read) else read.reference_start + 1

                        sv_type = 'ND'

                        if clipped_string == 'right' and read.reference_name == chr_SA:
                            if read.reference_end < pos_SA:
                                sv_type = 'DEL'
                            else:
                                sv_type = 'DUP'
                        elif clipped_string == 'left' and read.reference_name == chr_SA:
                            if read.reference_start > pos_SA:
                                sv_type = 'DEL'
                            else:
                                sv_type = 'DUP'
                        elif read.reference_name == chr_SA and strand_str[read.is_reverse] != strand_SA:
                            sv_type = 'INV'
                        elif read.reference_name != chr_SA:
                            sv_type = 'TRA'

                        dist = abs(clipped_pos - pos_SA) if read.reference_name == chr_SA else 0

                        # if sv_type == 'DEL':
                        #     print('{} {} {} => {}:{}-{} {} {} {}'.format(
                        #         sv_type, clipped_string, clipped_orient,
                        #         read.reference_name, clipped_pos,
                        #         pos_SA,
                        #         strand_str[read.is_reverse], strand_SA,
                        #         read.query_name
                        #     ))

                        split_pos_coord[sv_type] = append_coord(split_pos_coord[sv_type],
                                                                read.reference_name,
                                                                clipped_pos,
                                                                chr_SA,
                                                                pos_SA)

                        split_reads[read.reference_name][clipped_ch][clipped_pos] += 1

                        split_pos[clipped_string][read.reference_name].append(clipped_pos)
                        # print('adding {}:{}'.format(read.reference_name, clipped_pos))
                        split_read_distance[read.reference_name][clipped_ch][clipped_pos].append(dist)

    # Close the BAM file
    bamfile.close()

    # Look for INS positions:
    for chrom in chr_list:
        for k in ['right', 'left']:
            clipped_pos_dict[k][chrom] = set(clipped_pos_dict[k][chrom])

    # based on artificial INS
    for chrom in chr_list:
        for p in clipped_pos_dict['right'][chrom]:
            # is the right clipped position close to the left clipped position of a neighboring read?
            if len(set(range(p - 1, p + 1, 1)) & clipped_pos_dict['left'][chrom]) > 0:
                split_pos_coord['INS'] = append_coord(split_pos_coord['INS'],
                                                      chrom,
                                                      p,
                                                      chrom,
                                                      p + 1)
        for p in clipped_pos_dict['left'][chrom]:
            # is the right clipped position close to the left clipped position of a neighboring read?
            if len(set(range(p - 1, p + 1, 1)) & clipped_pos_dict['right'][chrom]) > 0:
                split_pos_coord['INS'] = append_coord(split_pos_coord['INS'],
                                                      chrom,
                                                      p,
                                                      chrom,
                                                      p + 1)

    # Count the number of split reads per position
    for chrom in chr_list:
        for k in ['right', 'left']:
            split_pos_cnt[k][chrom] = Counter(split_pos[k][chrom])

    for k in split_pos_coord.keys():
        split_pos_coord[k] = set(split_pos_coord[k])

    logging.info('Largest CIGAR "D" DEL={}'.format(max_cigar_del))

    logging.info('INDELs={}, split_reads={}, discordant_reads={}'.format(n_indels,
                                                                         n_split,
                                                                         n_discordant))

    for chrom in chr_list:
        for k in ['right', 'left']:
            logging.info('Number of unique ' + k + ' split read positions on Chr{}: {}'.format(chrom,
                                                                                               len(
                                                                                                   [p for p, c in
                                                                                                    split_pos_cnt[k][
                                                                                                        chrom].items()
                                                                                                    if
                                                                                                    c >= min_support]
                                                                                               )
                                                                                               ))
    # for chr1,pos1,chr2,pos2 in split_pos_coord['DEL']:
    #     print('{}:{}-{}:{}'.format(chr1, pos1, chr2, pos2))

    for k in split_pos_coord.keys():
        logging.info('Number of unique pair of split read positions ' + k + ': %d' % len(split_pos_coord[k]))

    for chrom in chr_list:
        for k in ['right', 'left']:
            total_reads_cnt[k][chrom] = Counter(split_pos[k][chrom])

    total_reads_coord = dict.fromkeys(split_pos_coord.keys())
    for k in total_reads_coord.keys():
        total_reads_coord[k] = list(set(split_pos_coord[k]))  # | discordant_reads_coord))

    for chrom in chr_list:
        for k in ['right', 'left']:
            positions_with_min_support[k][chrom] = [p for p, c in total_reads_cnt[k][chrom].items()
                                                    if c >= min_support]

            logging.info('Number of ' + k + '-split positions' +
                         ' on Chr%s with min %d support: %d' % (chrom, min_support,
                                                                len(
                                                                    positions_with_min_support[
                                                                        k][
                                                                        chrom])))
    for k in total_reads_coord.keys():
        logging.info('Number of unique pair of total positions ' + k + ': %d' % len(total_reads_coord[k]))

    positions_with_min_support_set = dict.fromkeys(chr_list)

    for chrom in chr_list:
        positions_with_min_support_set[chrom] = set(positions_with_min_support['left'][chrom] +
                                                    positions_with_min_support['right'][chrom])

    total_reads_coord_min_support = dict.fromkeys(total_reads_coord.keys())

    for k in total_reads_coord_min_support.keys():
        if k == 'INS':
            # INS positions are not based on split positions
            total_reads_coord_min_support[k] = total_reads_coord[k]
        else:
            total_reads_coord_min_support[k] = [(chr1, pos1, chr2, pos2) for chr1, pos1, chr2, pos2
                                                in total_reads_coord[k]
                                                if pos1 in positions_with_min_support_set[chr1] or
                                                pos2 in positions_with_min_support_set[chr2]
                                                ]
    for k in total_reads_coord_min_support.keys():
        logging.info(
            'Number of total pairs of ' + k + ' positions with min support: %d' % len(total_reads_coord_min_support[k]))

    data = (positions_with_min_support['left'], positions_with_min_support['right'], total_reads_coord_min_support,
            split_reads, split_read_distance)

    # Write JSON
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))

    #Write BEDPE
    with gzip.open(outBedpe, 'wt') as fout:
        for k in total_reads_coord_min_support.keys():
            for chr1, pos1, chr2, pos2 in total_reads_coord_min_support[k]:
                fout.write('\t'.join([chr1, str(pos1), str(pos1 + 1),
                                      chr2, str(pos2), str(pos2 + 1), k]) + '\n')


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
    parser.add_argument('-ob', '--outbedpe', type=str, default='split_reads.bedpe.gz',
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
    output_file_bedpe = os.path.join(output_dir, args.outbedpe)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    get_split_read_positions(ibam=args.bam, chr_list=args.chrlist, outFile=output_file, outBedpe=output_file_bedpe)

    logging.info('Time: split read positions on BAM %s: %f' % (args.bam, (time() - t0)))


if __name__ == '__main__':
    main()
