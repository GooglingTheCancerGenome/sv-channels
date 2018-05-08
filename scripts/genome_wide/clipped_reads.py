# Imports

import argparse
import pysam
import bz2file
import pickle
from time import time
import logging


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


def get_clipped_reads(ibam, chrName, outFile):
    '''

    :param ibam: input BAM alignment file
    :param chrName: chromosome name to consider
    :param outFile: output file for the dictionary of clipped reads
    :return: None
    '''

    # Check if the BAM file in input exists
    assert os.path.isfile(ibam)

    # Minimum read mapping quality to consider
    minMAPQ = 30

    # Dictionary to store number of clipped reads per position
    clipped_reads = dict()
    # For left- and right-clipped reads
    for split_direction in ['left', 'right']:
        clipped_reads[split_direction] = dict()

    # Open BAM file
    bamfile = pysam.AlignmentFile(ibam, "rb")
    # Get chromosome length from the BAM header
    header_dict = bamfile.header
    chrLen = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chrName][0]

    # Consider all the chromosome: interval [0, chrLen]
    start_pos = 0
    stop_pos = chrLen

    # Fetch the reads mapped on the chromosome
    iter = bamfile.fetch(chrName, start_pos, stop_pos)

    # Log information every n_r reads
    n_r = 10**6
    #print(n_r)
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

        # Both read and mate should be mapped
        if not read.is_unmapped and not read.mate_is_unmapped and read.mapping_quality >= minMAPQ:
            # Read is left-clipped
            if is_left_clipped(read):
                # print(str(read))
                # print('Clipped at the start: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                # print('Pos:%d, clipped_pos:%d' % (read.reference_start, read.get_reference_positions()[0]))
                # print('start:'+str(read.get_reference_positions()[0])+'=='+str(read.reference_start))
                ref_pos = read.reference_start
                if ref_pos not in clipped_reads['left'].keys():
                    clipped_reads['left'][ref_pos] = 1
                else:
                    clipped_reads['left'][ref_pos] += 1
            # Read is right-clipped
            if is_right_clipped(read):
                # print('Clipped at the end: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                # print('Pos:%d, clipped_pos:%d' %(read.reference_end, read.get_reference_positions()[-1]))
                # print('end: '+str(read.get_reference_positions()[-1]) + '==' + str(read.reference_end))
                ref_pos = read.reference_end + 1
                if ref_pos not in clipped_reads['right'].keys():
                    clipped_reads['right'][ref_pos] = 1
                else:
                    clipped_reads['right'][ref_pos] += 1

    # save clipped reads dictionary
    with bz2file.BZ2File(outFile, 'wb') as f:
        pickle.dump(clipped_reads, f)


def main():

    # Default BAM file for testing
    wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    inputBAM = wd + "T0_dedup.bam"
    # Default chromosome is 17 for the artificial data

    parser = argparse.ArgumentParser(description='Create channels with number of left/right clipped reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='clipped_reads.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='clipped_reads.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    get_clipped_reads(ibam=args.bam, chrName=args.chr, outFile=args.out)
    print('Time: clipped reads on BAM %s and Chr %s: %f' % (args.bam, args.chr, (time() - t0)))


if __name__ == '__main__':
    main()
