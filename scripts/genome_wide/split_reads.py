import argparse
import gzip
import json
import logging
import os
from collections import Counter, defaultdict
from time import time

import pysam
from cigar import Cigar

from functions import *

strand_str = {True: '-', False: '+'}


def append_coord(split_pos_coord, chr1, pos1, chr2, pos2, strand_info):
    if chr1 == chr2:
        if pos2 < pos1:
            split_pos_coord.append((chr2, pos2, chr1, pos1, strand_info))
        else:
            split_pos_coord.append((chr1, pos1, chr2, pos2, strand_info))
    elif chr2 < chr1:
        split_pos_coord.append((chr2, pos2, chr1, pos1, strand_info))
    elif chr2 > chr1:
        split_pos_coord.append((chr1, pos1, chr2, pos2, strand_info))
    return split_pos_coord


def get_split_read_positions(ibam, chr_list, min_mapq, min_sr_support, outFile, outBedpe):
    # List to store the split read positions
    split_pos_coord = dict()
    sv_type_list = ['INDEL_INS', 'INDEL_DEL',
                    'DEL', 'INS', 'INV', 'DUP', 'CTX', 'ND']
    for k in sv_type_list:
        split_pos_coord[k] = []

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    bam_mean, bam_stddev = get_insert_size(ibam, bamfile, min_mapq)
    n_indels = 0
    n_split = 0
    n_discordant = 0
    max_cigar_del = 0
    split_reads = dict()
    split_read_distance = dict()

    for chrom in chr_list:
        split_reads[chrom] = dict()
        split_read_distance[chrom] = dict()
        for split_direction in ['left_F', 'left_R', 'right_F', 'right_R', 'both_F', 'both_R']:
            split_reads[chrom][split_direction] = defaultdict(int)
            split_read_distance[chrom][split_direction] = defaultdict(list)

    clipped_pos_dict = dict()
    split_pos = dict()
    split_pos_cnt = dict()
    total_reads_cnt = dict()
    positions_with_min_support = dict()

    for k in ['right', 'left', 'both']:
        clipped_pos_dict[k] = defaultdict(list, {k: [] for k in chr_list})
        split_pos[k] = defaultdict(list, {k: [] for k in chr_list})
        split_pos_cnt[k] = dict.fromkeys(chr_list)
        total_reads_cnt[k] = dict.fromkeys(chr_list)
        positions_with_min_support[k] = dict.fromkeys(chr_list)

    n_r = 10 ** 6
    last_t = time()

    for i, read in enumerate(bamfile.fetch(), start=1):
        # Every n_r alignments, write log informations
        if not i % n_r:
            logging.info("%d alignments processed (%f alignments / s)" %
                         (i, n_r / (time() - last_t)))
            last_t = time()

        if is_left_clipped(read):
            clipped_pos_dict['left'][read.reference_name].append(
                read.reference_start + 1)
        elif is_right_clipped(read):
            clipped_pos_dict['right'][read.reference_name].append(
                read.reference_end)

        if not read.is_unmapped and read.mapping_quality >= min_mapq and \
                read.reference_name in chr_list:
            if has_indels(read):
                dels_start, dels_end, ins = get_indels(read)
                split_pos['left'][read.reference_name].extend(dels_end)
                split_pos['right'][read.reference_name].extend(dels_start)
                for pos in ins:
                    n_indels += 1
                    split_pos_coord['INDEL_INS'] = append_coord(
                        split_pos_coord['INDEL_INS'], read.reference_name, pos,
                        read.reference_name, pos + 1, '+-')
                for start, end in zip(dels_start, dels_end):
                    n_indels += 1
                    # Calculate DEL size and find largest DEL size encoded by CIGAR 'D' character
                    del_size = end - start
                    if del_size > max_cigar_del:
                        max_cigar_del = del_size
                    split_pos_coord['INDEL_DEL'] = append_coord(
                        split_pos_coord['INDEL_DEL'], read.reference_name, start,
                        read.reference_name, end, '+-')

            if read.has_tag('SA'):
                sa_entry = get_suppl_aln(read)
                if sa_entry is not None:
                    n_split += 1
                    chr_SA, pos_SA, strand_SA, cigar_sa = sa_entry
                    cigar_sa = Cigar(cigar_sa)
                    cigar_sa_list = list(cigar_sa.items())

                    if cigar_sa_list[-1][1] in ['S', 'H'] and not cigar_sa_list[0][1] in ['S', 'H']:
                        cigar_sa_len = 0
                        for l, c in cigar_sa_list:
                            # query consuming CIGAR characters
                            if c in ['M', 'D', 'N', '=', 'X']:
                                cigar_sa_len += l
                        pos_SA += cigar_sa_len
                    elif cigar_sa_list[0][1] in ['S', 'H'] and not cigar_sa_list[-1][1] in ['S', 'H']:
                        pos_SA += 1
                    elif cigar_sa_list[0][1] in ['S', 'H'] and cigar_sa_list[-1][1] in ['S', 'H']:
                        pos_SA += 1

                    if chr_SA in chr_list:
                        if is_right_clipped(read) and not is_left_clipped(read):
                            clipped_string = 'right'
                        if not is_right_clipped(read) and is_left_clipped(read):
                            clipped_string = 'left'
                        elif is_right_clipped(read) and is_left_clipped(read):
                            clipped_string = 'both'
                        clipped_orient = 'F' if not read.is_reverse else 'R'
                        clipped_ch = '_'.join([clipped_string, clipped_orient])
                        clipped_pos = read.reference_end if is_right_clipped(
                            read) else read.reference_start + 1
                        sv_type = 'ND'
                        if clipped_string == 'right' and read.reference_name == chr_SA:
                            if read.reference_start < pos_SA:
                                sv_type = 'DEL'
                            else:
                                sv_type = 'DUP'
                        elif clipped_string == 'left' and read.reference_name == chr_SA:
                            if read.reference_start > pos_SA:
                                sv_type = 'DUP'
                            else:
                                sv_type = 'DEL'
                        elif clipped_string == 'both' and read.reference_name == chr_SA:
                            sv_type = 'DEL'

                        if read.reference_name == chr_SA and \
                                strand_str[read.is_reverse] == strand_str[read.mate_is_reverse]:
                            sv_type = 'INV'

                        if read.reference_name != chr_SA:
                            sv_type = 'CTX'

                        dist = abs(
                            clipped_pos - pos_SA) if read.reference_name == chr_SA else 0
                        dist = (dist - bam_mean) / \
                            bam_stddev if read.reference_name == chr_SA else 0

                        if strand_str[read.is_reverse] == strand_SA:
                            strand_info = '+-'
                        else:
                            strand_info = strand_str[read.is_reverse]+strand_SA

                        split_pos_coord[sv_type] = append_coord(split_pos_coord[sv_type],
                                                                read.reference_name,
                                                                clipped_pos,
                                                                chr_SA,
                                                                pos_SA,
                                                                strand_info)

                        split_reads[read.reference_name][clipped_ch][clipped_pos] += 1
                        split_pos[clipped_string][read.reference_name].append(
                            clipped_pos)
                        split_read_distance[read.reference_name][clipped_ch][clipped_pos].append(
                            dist)
    bamfile.close()

    # Look for INS positions:
    for chrom in chr_list:
        for k in ['right', 'left']:
            clipped_pos_dict_cnt = Counter(clipped_pos_dict[k][chrom])
            clipped_pos_dict[k][chrom] = {
                int(key) for key, val in clipped_pos_dict_cnt.items() if val >= 3}

    # based on artificial INS
    for chrom in chr_list:
        for p in clipped_pos_dict['right'][chrom]:
            # is the right clipped position close to the left clipped position of a neighboring read?
            if len(set(range(p - 1, p + 1, 1)) and clipped_pos_dict['left'][chrom]) > 0:
                split_pos_coord['INS'] = append_coord(
                    split_pos_coord['INS'], chrom, p, chrom, p + 1, '+-')
        for p in clipped_pos_dict['left'][chrom]:
            # is the right clipped position close to the left clipped position of a neighboring read?
            if len(set(range(p - 1, p + 1, 1)) and clipped_pos_dict['right'][chrom]) > 0:
                split_pos_coord['INS'] = append_coord(
                    split_pos_coord['INS'], chrom, p, chrom, p + 1, '+-')
    # Count the number of split reads per position
    for chrom in chr_list:
        for k in ['right', 'left']:
            split_pos_cnt[k][chrom] = Counter(split_pos[k][chrom])

    for k in split_pos_coord.keys():
        split_pos_coord[k] = set(split_pos_coord[k])

    logging.info('Largest CIGAR "D" DEL={}'.format(max_cigar_del))
    logging.info('INDELs={}, split_reads={}, discordant_reads={}'.format(
        n_indels, n_split, n_discordant))

    for chrom in chr_list:
        for k in ['right', 'left']:
            logging.info("Number of unique %s-split read positions on Chr%s: %d" % (k, str(chrom),
                                                                                    len([p for p, c in split_pos_cnt[k][chrom].items() if c >= min_sr_support])))

    for k in split_pos_coord.keys():
        logging.info('Number of unique pair of split read positions ' +
                     k + ': %d' % len(split_pos_coord[k]))

    for chrom in chr_list:
        for k in ['right', 'left']:
            total_reads_cnt[k][chrom] = Counter(split_pos[k][chrom])

    total_reads_coord = dict.fromkeys(split_pos_coord.keys())
    for k in total_reads_coord.keys():
        total_reads_coord[k] = list(set(
            split_pos_coord[k]))  # | discordant_reads_coord))

    for chrom in chr_list:
        for k in ['right', 'left']:
            positions_with_min_support[k][chrom] = [
                p for p, c in total_reads_cnt[k][chrom].items()
                if c >= min_sr_support
            ]
            logging.info("Number of %s-split positions on Chr%s with min %d support: %d" %
                         (k, str(chrom), min_sr_support, len(positions_with_min_support[k][chrom])))

    for k in total_reads_coord.keys():
        logging.info("Number of unique pair of total positions %s: %d" %
                     (k, len(total_reads_coord[k])))
    positions_with_min_support_set = dict.fromkeys(chr_list)

    for chrom in chr_list:
        positions_with_min_support_set[chrom] = set(
            positions_with_min_support['left'][chrom] +
            positions_with_min_support['right'][chrom])
    total_reads_coord_min_support = dict.fromkeys(total_reads_coord.keys())

    for k in total_reads_coord_min_support.keys():
        if k == 'INS':
            # INS positions are not based on split positions
            total_reads_coord_min_support[k] = total_reads_coord[k]
        else:
            total_reads_coord_min_support[k] = [
                (chr1, pos1, chr2, pos2, strand_info)
                for chr1, pos1, chr2, pos2, strand_info in total_reads_coord[k]
                if pos1 in positions_with_min_support_set[chr1] or pos2 in positions_with_min_support_set[chr2]]
    for k in total_reads_coord_min_support.keys():
        logging.info("Number of total pairs of %s positions with min support: %d" % (
            k, len(total_reads_coord_min_support[k])))

    data = (positions_with_min_support['left'],
            positions_with_min_support['right'], total_reads_coord_min_support,
            split_reads, split_read_distance)

    # Write JSON
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))

    # Write BEDPE
    with gzip.open(outBedpe, 'wt') as fout:
        for k in total_reads_coord_min_support.keys():
            for chr1, pos1, chr2, pos2, strand_info in total_reads_coord_min_support[k]:
                assert len(strand_info) == 2
                fout.write('\t'.join([
                    chr1,
                    str(pos1),
                    str(pos1 + 1), chr2,
                    str(pos2),
                    str(pos2 + 1), k, '*', strand_info[0], strand_info[1]
                ]) + '\n')


def main():
    parser = argparse.ArgumentParser(description='Get split reads positions')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chrlist',
                        type=str,
                        default='12,22',
                        help="Comma separated list of chromosomes to consider")
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        default='split_reads.json.gz',
                        help="Specify output")
    parser.add_argument('-ob',
                        '--outbedpe',
                        type=str,
                        default='split_reads.bedpe.gz',
                        help="Specify output")
    parser.add_argument(
        '-p',
        '--outputpath',
        type=str,
        default='.',
        help="Specify output path"
    )
    parser.add_argument('-l',
                        '--logfile',
                        default='split_reads.log',
                        help='File in which to write logs.')
    parser.add_argument('-m',
                        '--min_mapq',
                        type=int,
                        default=10,
                        help='Minimum read mapping quality')
    parser.add_argument('-s',
                        '--min_sr_support',
                        type=int,
                        default=1,
                        help='Minimum number of split reads')

    args = parser.parse_args()

    # Log file
    cmd_name = 'split_reads'
    output_dir = os.path.join(args.outputpath, cmd_name)
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, args.out)
    output_file_bedpe = os.path.join(output_dir, args.outbedpe)
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    get_split_read_positions(ibam=args.bam,
                             chr_list=args.chrlist.split(','),
                             min_mapq=args.min_mapq,
                             min_sr_support=args.min_sr_support,
                             outFile=output_file,
                             outBedpe=output_file_bedpe)
    logging.info('Time: split read positions on BAM %s: %f' %
                 (args.bam, (time() - t0)))


if __name__ == '__main__':
    main()
