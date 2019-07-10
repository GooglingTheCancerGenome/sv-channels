# Imports
import argparse
import logging
import os
from collections import Counter
from time import time
import pysam
from functions import *
import numpy as np
import twobitreader as twobit

config = get_config_file()

HPC_MODE = config["DEFAULT"]["HPC_MODE"]
MAX_PILEUP_BUFFER_SIZE = 8000
minMAPQ = config["DEFAULT"]["MIN_MAPQ"]


def get_snvs(ibam, chrName, outFile):
    def get_2bit_genome():
        if HPC_MODE:
            # Path on the HPC of the 2bit version of the human reference genome (hg19)
            genome = twobit.TwoBitFile('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes/hg19.2bit')
        else:
            # Path on the local machine of the 2bit version of the human reference genome (hg19)
            genome = twobit.TwoBitFile('/Users/lsantuari/Documents/Data/GiaB/reference/hg19.2bit')
        return genome

    def get_snv_number(query_seq_list, reference_base):

        reference_base = reference_base.upper()
        if len(query_seq_list) > 0 and reference_base != 'N':
            cnt = Counter(list(map(lambda x: x.upper(), query_seq_list)))
            return cnt['A'] + cnt['T'] + cnt['C'] + cnt['G'] - cnt[reference_base]
        else:
            return 0

    def correct_for_indels(read, ref_string):

        seq = read.query_alignment_sequence
        base_quals = [i for i in read.query_alignment_qualities]
        # print(base_quals)
        bq = []
        pos = 0
        bq.append(read.reference_start + i)
        for op, l in read.cigartuples:
            if op not in [0, 1, 2]:
                continue
            elif op == 1:
                ref_string = ref_string[:pos] + '-' * l + ref_string[pos:]
                base_quals = base_quals[:pos] + base_quals[pos + l:]
            elif op == 2:
                seq = seq[:pos] + '-' * l + seq[pos:]
                base_quals = base_quals[:pos] + [0] * l + base_quals[pos:]
                pos += l
        return seq, base_quals, ref_string

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    # Load the BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Extract the header
    header_dict = bamfile.header
    # Get the chromosome length from the header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # stop_pos = 10000000

    reference_sequence = get_2bit_genome()

    # snv_list = ['BQ', 'nALN', 'nSEG', 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't']
    snv_list = ['BQ', 'SNV', 'RC']
    snv_array = np.zeros(shape=(len(snv_list), stop_pos + 1), dtype=np.uint32)
    snv_dict = {v: n for n, v in enumerate(snv_list)}
    # print(snv_dict)
    # Print every n_r alignments processed
    n_r = 10 ** 6
    # Record the current time
    last_t = time()

    for pileupcolumn in bamfile.pileup(chrName, start_pos, stop_pos, stepper='all'):
        # pileupcolumn.set_min_base_quality(0)
        # print("\ncoverage at base %s = %s" %
        #       (pileupcolumn.pos, pileupcolumn.nsegments))
        if 0 < pileupcolumn.nsegments < MAX_PILEUP_BUFFER_SIZE and start_pos <= pileupcolumn.pos <= stop_pos:
            quals = pileupcolumn.get_query_qualities()
            if len(quals) > 0:
                snv_array[snv_dict['BQ'], pileupcolumn.pos] = np.median(pileupcolumn.get_query_qualities())
            # snv_array[snv_dict['nALN'], pileupcolumn.pos] = pileupcolumn.get_num_aligned()
            # snv_array[snv_dict['nSEG'], pileupcolumn.pos] = pileupcolumn.nsegments
            try:

                query_seq_list = pileupcolumn.get_query_sequences()
                snv_number = get_snv_number(query_seq_list, reference_sequence['chr' + chrName][pileupcolumn.pos])
                snv_array[snv_dict['SNV'], pileupcolumn.pos] = snv_number/pileupcolumn.nsegments \
                    if pileupcolumn.nsegments !=0 else 0

            except AssertionError as error:
                # Output expected AssertionErrors.
                logging.info(error)
                logging.info('Position {}:{} has {} nsegments'.format(chrName,
                                                                      pileupcolumn.pos,
                                                                      pileupcolumn.nsegments))
                continue

    # # Pysam iterator to fetch the reads
    # iter = bamfile.fetch(chrName, start_pos, stop_pos)
    #
    # for i, read in enumerate(iter, start=1):
    #
    #     # Every n_r alignments, write log informations
    #     if not i % n_r:
    #         # Record the current time
    #         now_t = time()
    #         # print(type(now_t))
    #         logging.info("%d alignments processed (%f alignments / s)" % (
    #             i,
    #             n_r / (now_t - last_t)))
    #         last_t = time()
    #
    #     if not read.is_unmapped and read.mapping_quality >= minMAPQ:
    #
    #         # print(read)
    #         # print(read.cigartuples)
    #         # print(read.mapping_quality)
    #         # print(read.query_alignment_qualities)
    #         # print('{} == {}'.format(
    #         #     len(read.query_alignment_qualities), read.reference_end - read.reference_start))
    #
    #         try:
    #
    #             ref_string = reference_sequence['chr' + chrName][read.reference_start:read.reference_end].upper()
    #
    #             if len(read.cigartuples) == 1 and read.cigartuples[0][0] == 0:
    #
    #                 query_string = read.query_alignment_sequence
    #                 query_base_qualities = read.query_alignment_qualities
    #
    #             else:
    #
    #                 query_string, query_base_qualities, ref_string = correct_for_indels(read, ref_string)
    #                 # print(read.query_alignment_sequence)
    #                 # print(query_string)
    #                 # print(ref_string)
    #
    #                 snv_idx = [read.reference_start + index for index, (e1, e2) in enumerate(
    #                     zip(query_string, ref_string)) if e1 != e2 and e1 != '-' and e2 != '-' ]
    #                 # for i in snv_idx:
    #                 #     pos = i - read.reference_start
    #                 #     print('{} != {}'.format(query_string[pos], ref_string[pos]))
    #
    #                 snv_array[snv_dict['SNV'], snv_idx] += 1
    #
    #             snv_array[snv_dict['BQ'], read.reference_start:read.reference_end] += \
    #                 np.array(query_base_qualities, dtype=np.uint32)
    #             snv_array[snv_dict['RC'], read.reference_start:read.reference_end] += 1
    #
    #         except:
    #
    #             print(read)
    #             print(read.cigartuples)
    #             print(read.mapping_quality)
    #             print(read.query_alignment_qualities)

            # print(cnt)
            # for k in cnt.keys():
            #     if k in snv_list and k.upper() != reference_sequence['chr' + chrName][pileupcolumn.pos]:
            #         snv_array[snv_dict[k], pileupcolumn.pos] = cnt[k]

    # # Close the BAM file
    # bamfile.close()
    #
    # snv_array = np.vstack(
    #     (
    #         np.divide(snv_array[snv_dict['SNV']], snv_array[snv_dict['RC']],
    #                   where=snv_array[snv_dict['RC']] != 0),
    #         np.divide(snv_array[snv_dict['BQ']], snv_array[snv_dict['RC']],
    #                   where=snv_array[snv_dict['RC']] != 0)
    #     )
    # )

    # Write the output
    np.save(file=outFile, arr=snv_array)
    # os.system('gzip -f ' + outFile)


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
    parser = argparse.ArgumentParser(description='Get snv info')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='snv.npy',
                        help="Specify output")
    parser.add_argument('-p', '--outputpath', type=str,
                        default='/Users/lsantuari/Documents/Processed/channel_maker_output',
                        help="Specify output path")
    parser.add_argument('-l', '--logfile', default='snv.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file

    cmd_name = 'snv'
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
    get_snvs(ibam=args.bam, chrName=args.chr, outFile=output_file)
    logging.info('Time: SNVs on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
