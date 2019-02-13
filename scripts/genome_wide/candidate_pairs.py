# Imports
import argparse
import logging
import os
import pickle
import bz2file
from collections import Counter, defaultdict
from time import time
import pysam
import functions as fun


class Breakpoint:

    def __init__(self, chr, pos, strand):
        self.chr = chr
        self.pos = pos
        self.strand = strand

    def id(self):
        return '_'.join([self.chr, str(self.pos), self.strand])


class StructuralVariant:

    def __init__(self, bp1, bp2):

        if bp1.chr == bp2.chr:
            if bp1.pos <= bp2.pos:
                self.tuple = (bp1, bp2)
            else:
                self.tuple = (bp2, bp1)
        elif bp1.chr < bp2.chr:
            self.tuple = (bp1, bp2)
        else:
            self.tuple = (bp2, bp1)


strand = {False: '+', True: '-'}


def get_candidate_breakpoint_pairs(ibam, chrName, outFile):
    '''
    
    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the dictionary of clipped read positions
    :return: None. Outputs a dictionary with the positions of clipped read positions as keys and
    the number of clipped reads per position as values
    '''''

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    # Minimum read mapping quality to consider
    minMAPQ = 30

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # List to store the clipped read positions
    candidate_pairs = set()

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # Pysam iterator to fetch the reads
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

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

        # Both read and mate should be mapped, read should have a minimum mapping quality
        # if (not read.is_unmapped) and (not read.mate_is_unmapped) and read.mapping_quality >= minMAPQ:
        if not read.is_unmapped and read.mapping_quality >= minMAPQ:

            if fun.has_indels(read):
                # print(read)
                dels_start, dels_end, ins = fun.get_indels(read)
                for deletion in zip(dels_start, dels_end):
                    candidate_pairs.add(StructuralVariant(Breakpoint(chrName, deletion[0], strand[read.is_reverse]),
                                                          Breakpoint(chrName, deletion[1], strand[read.is_reverse])))
                for insertion in ins:
                    candidate_pairs.add(StructuralVariant(Breakpoint(chrName, insertion, strand[read.is_reverse]),
                                                          Breakpoint(chrName, insertion + 1, strand[read.is_reverse])))

            if fun.is_left_clipped(read):
                cpos = read.reference_start
                if fun.has_suppl_aln(read):
                    chr_sa, pos_sa, strand_sa = fun.get_suppl_aln(read)
                    candidate_pairs.add(StructuralVariant(Breakpoint(chrName, cpos, strand[read.is_reverse]),
                                                          Breakpoint(chr_sa, pos_sa, strand_sa)))
                else:
                    if not read.mate_is_unmapped:
                        candidate_pairs.add(StructuralVariant(Breakpoint(chrName, cpos, strand[read.is_reverse]),
                                                              Breakpoint(read.next_reference_name,
                                                                         read.next_reference_start,
                                                                         strand[read.mate_is_reverse])))

            if fun.is_right_clipped(read):
                cpos = read.reference_end + 1
                if fun.has_suppl_aln(read):
                    chr_sa, pos_sa, strand_sa = fun.get_suppl_aln(read)
                    candidate_pairs.add(StructuralVariant(Breakpoint(chrName, cpos, strand[read.is_reverse]),
                                                          Breakpoint(chr_sa, pos_sa, strand_sa)))
                else:
                    if not read.mate_is_unmapped:
                        candidate_pairs.add(StructuralVariant(Breakpoint(chrName, cpos, strand[read.is_reverse]),
                                                              Breakpoint(read.next_reference_name,
                                                                         read.next_reference_start,
                                                                         strand[read.mate_is_reverse])))

    # Close the BAM file
    bamfile.close()

    inspect_pairs(candidate_pairs, outFile)


def inspect_pairs(candidate_pairs, outFile):

    final_pairs = set()

    # from bp1 point of view
    bp_dict = defaultdict(dict)
    bp_list = []
    for sv in candidate_pairs:
        bp1, bp2 = sv.tuple

        bp_id = bp1.id()
        bp2_id = '_'.join([bp2.chr, bp2.strand])
        if bp2_id not in bp_dict[bp_id]:
            bp_dict[bp_id] = defaultdict(list)
        bp_dict[bp_id][bp2_id].append(bp2.pos)
        bp_list.append(bp_id)

    bp_cnt = Counter(bp_list)
    min_support_bp = [k for (k, v) in bp_cnt.items() if v >= 3]
    print('min support positions bp1; %d/%d' %
          (len(min_support_bp), len(bp_cnt)))
    for bp1_id in min_support_bp:
        bp1_chr, bp1_pos, bp1_strand = bp1_id.split('_')
        for bp2_id in bp_dict[bp1_id]:
            bp2_chr, bp2_strand = bp2_id.split('_')
            bp2_pos = max(bp_dict[bp1_id][bp2_id]) if bp1_strand == '+' else min(bp_dict[bp1_id][bp2_id])
            final_pairs.add(StructuralVariant(Breakpoint(bp1_chr, int(bp1_pos), bp1_strand),
                                              Breakpoint(bp2_chr, int(bp2_pos), bp2_strand)))

    print(len(final_pairs))

    # from bp2 point of view
    bp_dict = defaultdict(dict)
    bp_list = []
    for sv in candidate_pairs:
        bp1, bp2 = sv.tuple

        bp_id = bp2.id()
        bp1_id = '_'.join([bp1.chr, bp1.strand])
        if bp1_id not in bp_dict[bp_id]:
            bp_dict[bp_id] = defaultdict(list)
        bp_dict[bp_id][bp1_id].append(bp1.pos)
        bp_list.append(bp_id)

    bp_cnt = Counter(bp_list)
    min_support_bp = [k for (k, v) in bp_cnt.items() if v >= 3]
    print('unique positions bp2; %d/%d' %
          (len(min_support_bp), len(bp_cnt)))
    for bp1_id in min_support_bp:
        bp1_chr, bp1_pos, bp1_strand = bp1_id.split('_')
        for bp2_id in bp_dict[bp1_id]:
            bp2_chr, bp2_strand = bp2_id.split('_')
            bp2_pos = max(bp_dict[bp1_id][bp2_id]) if bp1_strand == '+' else min(bp_dict[bp1_id][bp2_id])
            final_pairs.add(StructuralVariant(Breakpoint(bp1_chr, int(bp1_pos), bp1_strand),
                                              Breakpoint(bp2_chr, int(bp2_pos), bp2_strand)))

    print(len(final_pairs))

    # Write the output in pickle format
    with bz2file.BZ2File(outFile, 'wb') as f:
        pickle.dump(final_pairs, f)


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
    parser = argparse.ArgumentParser(description='Get clipped reads positions')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='candidate_pairs.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='candidate_pairs.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()

    get_candidate_breakpoint_pairs(ibam=args.bam, chrName=args.chr, outFile=args.out)

    print('Time: candidate breakpoint pairs on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
