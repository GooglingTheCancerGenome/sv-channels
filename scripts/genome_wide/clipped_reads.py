import argparse
import gzip
import json
import logging
import os
from collections import defaultdict
from time import time

import pysam

from functions import *


def get_clipped_reads(ibam, chr_list, minMAPQ, outFile):
    '''
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the dictionary of clipped reads
    :return: None
    '''
    # Dictionary to store number of clipped reads per position
    clipped_reads = dict()
    clipped_reads_inversion = dict()
    clipped_reads_duplication = dict()
    clipped_reads_translocation = dict()

    for chrom in chr_list:
        clipped_reads[chrom] = dict()
        # For left- and right-clipped reads
        for split_direction in [
            'left_F', 'left_R', 'right_F', 'right_R', 'disc_right_F',
            'disc_right_R', 'disc_left_F', 'disc_left_R', 'D_left_F',
            'D_left_R', 'D_right_F', 'D_right_R', 'I_F', 'I_R'
        ]:
            clipped_reads[chrom][split_direction] = defaultdict(int)

        # Dictionary to store number of clipped reads per position for
        # INVersion:
        # reads that are clipped AND mapped on the same chromosome AND with the same orientation (FF or RR)
        # Two channels: mate is mapped before or after the read
        clipped_reads_inversion[chrom] = dict()

        # DUPlication:
        # 1) reads that are right-clipped AND mapped on the same chromosome
        # AND read is forward AND mate is reverse AND mate is mapped before read

        # 2) reads that are left-clipped AND mapped on the same chromosome
        # AND read is reverse AND mate is forward AND mate is mapped after read
        clipped_reads_duplication[chrom] = dict()

        # TRAslocation:
        # Two channels: reads with mate mapped to a different chromosome and with
        # 1: opposite orientation
        # 2: same orientation
        clipped_reads_translocation[chrom] = dict()

        # Mate is mapped before or after?
        for mate_position in ['before', 'after', 'before_split', 'after_split']:
            clipped_reads_inversion[chrom][mate_position] = defaultdict(int)
            clipped_reads_duplication[chrom][mate_position] = defaultdict(int)

        for orientation in ['opposite', 'same', 'opposite_split', 'same_split']:
            clipped_reads_translocation[chrom][orientation] = defaultdict(int)

    # Open BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Log information every n_r reads
    n_r = 10 ** 6
    last_t = time()

    for i, read in enumerate(bamfile.fetch(), start=1):
        if read.reference_name in chr_list:
            if not i % n_r:
                logging.info("%d alignments processed (%f alignments / s)" %
                             (i, n_r / (time() - last_t)))
                last_t = time()

            if not read.is_unmapped and read.mapping_quality >= minMAPQ:
                if has_indels(read):
                    dels_start, dels_end, ins = get_indels(read)
                    for del_pos in dels_start:
                        if not read.is_reverse:
                            clipped_reads[
                                read.reference_name]['D_left_F'][del_pos] += 1
                        else:
                            clipped_reads[
                                read.reference_name]['D_left_R'][del_pos] += 1
                    for del_pos in dels_end:
                        if not read.is_reverse:
                            clipped_reads[
                                read.reference_name]['D_right_F'][del_pos] += 1
                        else:
                            clipped_reads[
                                read.reference_name]['D_right_R'][del_pos] += 1

                    for ins_pos in ins:
                        if not read.is_reverse:
                            clipped_reads[
                                read.reference_name]['I_F'][ins_pos] += 1
                        else:
                            clipped_reads[
                                read.reference_name]['I_R'][ins_pos] += 1

            # Both read and mate should be mapped, with mapping quality greater than minMAPQ
            if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ:
                if is_left_clipped(read):
                    ref_pos = read.reference_start + 1
                elif is_right_clipped(read):
                    ref_pos = read.reference_end

                if read.reference_name == read.next_reference_name:
                    if read.is_reverse != read.mate_is_reverse:
                        # Read is left-clipped
                        if is_left_clipped(read):
                            if not has_suppl_aln(read):
                                if not read.is_reverse:
                                    clipped_reads[read.reference_name]['left_F'][
                                        ref_pos] += 1
                                    if not read.is_proper_pair:
                                        clipped_reads[read.reference_name][
                                            'disc_left_F'][ref_pos] += 1
                                else:
                                    clipped_reads[read.reference_name]['left_R'][
                                        ref_pos] += 1
                                    if not read.is_proper_pair:
                                        clipped_reads[read.reference_name][
                                            'disc_left_R'][ref_pos] += 1

                            # DUPlication, channel 2
                            # Read is mapped on the Reverse strand and mate is mapped on the Forward strand
                            if read.is_reverse and not read.mate_is_reverse \
                                    and read.reference_start < read.next_reference_start:  # Mate is mapped after read
                                if not has_suppl_aln(read):
                                    clipped_reads_duplication[
                                        read.reference_name]['after'][ref_pos] += 1
                                else:
                                    clipped_reads_duplication[
                                        read.reference_name]['after_split'][ref_pos] += 1

                        # Read is right-clipped
                        elif is_right_clipped(read):
                            if not has_suppl_aln(read):
                                if not read.is_reverse:
                                    clipped_reads[read.reference_name]['right_F'][
                                        ref_pos] += 1
                                    if not read.is_proper_pair:
                                        clipped_reads[read.reference_name][
                                            'disc_right_F'][ref_pos] += 1
                                else:
                                    clipped_reads[read.reference_name]['right_R'][
                                        ref_pos] += 1
                                    if not read.is_proper_pair:
                                        clipped_reads[read.reference_name][
                                            'disc_right_R'][ref_pos] += 1

                            # DUPlication, channel 1
                            # Read is mapped on the Forward strand and mate is mapped on the Reverse strand
                            if not read.is_reverse and read.mate_is_reverse:
                                # Mate is mapped before read
                                if read.reference_start > read.next_reference_start:
                                    if not has_suppl_aln(read):
                                        clipped_reads_duplication[read.reference_name]['before'][ref_pos] += 1
                                    else:
                                        clipped_reads_duplication[read.reference_name]['before_split'][ref_pos] += 1

                        # The following if statement takes care of the inversion channels
                        # Read and mate are mapped on the same strand: either Forward-Forward or Reverse-Reverse

                    elif read.is_reverse == read.mate_is_reverse:
                        if is_clipped(read) and not has_suppl_aln(read):
                            # Mate is mapped before read
                            if read.reference_start > read.next_reference_start:
                                if not has_suppl_aln(read):
                                    clipped_reads_inversion[read.reference_name][
                                        'before'][ref_pos] += 1
                                else:
                                    clipped_reads_inversion[read.reference_name][
                                        'before_split'][ref_pos] += 1
                            # Mate is mapped after read
                            else:
                                if not has_suppl_aln(read):
                                    clipped_reads_inversion[
                                        read.reference_name]['after'][ref_pos] += 1
                                else:
                                    clipped_reads_inversion[
                                        read.reference_name]['after_split'][ref_pos] += 1

                else:
                    if is_clipped(read):
                        if read.is_reverse != read.mate_is_reverse:
                            if not has_suppl_aln(read):
                                clipped_reads_translocation[
                                    read.reference_name]['opposite'][ref_pos] += 1
                            else:
                                clipped_reads_translocation[
                                    read.reference_name]['opposite_split'][ref_pos] += 1
                        else:
                            if not has_suppl_aln(read):
                                clipped_reads_translocation[
                                    read.reference_name]['same'][ref_pos] += 1
                            else:
                                clipped_reads_translocation[
                                    read.reference_name]['same_split'][ref_pos] += 1

    # Write clipped reads dictionaries
    data = (clipped_reads, clipped_reads_inversion, clipped_reads_duplication,
            clipped_reads_translocation)
    with gzip.GzipFile(outFile, 'w') as fout:
        fout.write(json.dumps(data).encode('utf-8'))


def main():
    parser = argparse.ArgumentParser(
        description='Create channels with number of left/right clipped reads')
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
                        default='clipped_reads.json.gz',
                        help="Specify output")
    parser.add_argument('-m',
                        '--min_mapq',
                        type=int,
                        default=10,
                        help='Minimum read mapping quality')
    parser.add_argument(
        '-p',
        '--outputpath',
        type=str,
        default='.',
        help="Specify output path")
    parser.add_argument('-l',
                        '--logfile',
                        default='clipped_reads.log',
                        help='File in which to write logs.')
    args = parser.parse_args()
    cmd_name = 'clipped_reads'
    output_dir = os.path.join(args.outputpath, cmd_name)
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, args.out)
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)
    t0 = time()
    get_clipped_reads(ibam=args.bam,
                      chr_list=args.chrlist.split(','),
                      minMAPQ=args.min_mapq,
                      outFile=output_file)
    logging.info('Time: clipped reads on BAM %s: %f' %
                 (args.bam, (time() - t0)))


if __name__ == '__main__':
    main()
