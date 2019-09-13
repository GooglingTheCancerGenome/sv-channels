'''

This file contains the functions to compute the channel for the exact coverage. In the current implementation of
ChannelMaker, the exact coverage was discarded because considered unnecessary, given that the coverage channel is
already present. Computing the exact coverage requires too much time.

'''

#Imports
import argparse
from time import time
import twobitreader as twobit
import logging
import h5py
import numpy as np
from time import time

def get_one_hot_genome():
    '''

        :param ibam: input BAM alignment file
        :param chrName: chromosome name to consider
        :param outFile: output file with the list containing the exact coverage per chromosome position. Exact coverage
        is the number of bases mapped in a certain position that are equal to the reference base.
        :return: None. It only writes the exact coverage in the output file
    '''

    # Path on the local machine of the 2bit version of the human reference genome (hg19)
    genome = twobit.TwoBitFile('/Users/lsantuari/Documents/Data/GiaB/reference/hg19.2bit')
    # Path on the HPC of the 2bit version of the human reference genome (hg19)
    #genome = twobit.TwoBitFile('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes/hg19.2bit')

    #print('chr'+list(map(str, range(1,23))).extend('X','Y','M'))

    chrlist = list(map(str, range(1, 23)))
    chrlist.extend(['X','Y','M'])

    def vectorizeSequence(seq):
        # the order of the letters is not arbitrary.
        # Flip the matrix up-down and left-right for reverse compliment
        ltrdict = {'a': [1, 0, 0, 0], 'c': [0, 1, 0, 0], 'g': [0, 0, 1, 0], 't': [0, 0, 0, 1], 'n': [0, 0, 0, 0]}
        return np.array([ltrdict[x] for x in seq])

    starttime = time()
    with h5py.File('genomeEncoded.h5', 'w') as hf:
        for chrname in chrlist:
            fasta_name = 'chr'+chrname
            if fasta_name in genome.keys():
                # get the fasta files.
                sequence = str(genome[fasta_name])
                # one hot encode...
                data = vectorizeSequence(sequence.lower())
                print(fasta_name + " is one hot encoded!")
                # write to hdf5
                hf.create_dataset(chrname, data=data)
                print(chrname + " is written to dataset")

    endtime = time()
    print("Encoding is done in " + str(endtime))


def main():

    # Default BAM file for testing
    #wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create exact coverage channel')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-c', '--chr', type=str, default='17',
                        help="Specify chromosome")
    parser.add_argument('-o', '--out', type=str, default='exact_coverage.pbz2',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='exact_coverage.log',
                        help='File in which to write logs.')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    # Compute the exact coverage
    get_one_hot_genome()
    print(time() - t0)


if __name__ == '__main__':
    main()
