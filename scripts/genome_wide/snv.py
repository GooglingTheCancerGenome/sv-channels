# Imports
import argparse
import logging
import os
import pickle
import bz2file
from collections import Counter
from time import time
import pysam
import functions as fun
import numpy as np
import twobitreader as twobit
import json

with open('./genome_wide/parameters.json', 'r') as f:
    config = json.load(f)

HPC_MODE = config["DEFAULT"]["HPC_MODE"]
MAX_PILEUP_BUFFER_SIZE = 8000

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

        if len(query_seq_list) > 0:
            cnt = Counter(list(map(lambda x:x.upper(), query_seq_list)))
            return cnt['A']+cnt['T']+cnt['C']+cnt['G']-cnt[reference_base]
        else:
            return 0


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

    # Fetch reads over the entire chromosome between positions [0, chrLen]
    start_pos = 0
    stop_pos = chrLen
    # stop_pos = 10000000

    reference_sequence = get_2bit_genome()

    # snv_list = ['BQ', 'nALN', 'nSEG', 'A', 'a', 'C', 'c', 'G', 'g', 'T', 't']
    snv_list = ['BQ', 'SNV']
    snv_array = np.zeros(shape=(len(snv_list), stop_pos+1), dtype=float)
    snv_dict = {v: n for n, v in enumerate(snv_list)}
    # print(snv_dict)
    # Print every n_r alignments processed
    n_r = 10 ** 6
    # Record the current time
    last_t = time()

    for pileupcolumn in bamfile.pileup(chrName, start_pos, stop_pos, min_mapping_quality=minMAPQ):
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

            # print(cnt)
            # for k in cnt.keys():
            #     if k in snv_list and k.upper() != reference_sequence['chr' + chrName][pileupcolumn.pos]:
            #         snv_array[snv_dict[k], pileupcolumn.pos] = cnt[k]


    # Close the BAM file
    bamfile.close()

    # Write the output
    np.savez(outFile, snv_array=snv_array)
    os.system('gzip -f ' + outFile)


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
    parser.add_argument('-o', '--out', type=str, default='snv.npz',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='snv.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    # Log file
    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()
    get_snvs(ibam=args.bam, chrName=args.chr, outFile=args.out)
    logging.info('Time: SNVs on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()