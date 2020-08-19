import pysam
import logging
import argparse
from time import time
import os


def extract_clipped_reads(ibam, obam, ofastq):

    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    if os.path.exists(ofastq):
        os.remove(ofastq)

    aln = pysam.AlignmentFile(ibam, "rb")
    out_aln = pysam.AlignmentFile(obam, "wb", template=aln)

    for read in aln.fetch():
        # write full sequence and quality of soft clipped reads in a FASTQ file
        if read.cigartuples is not None and \
                (
                        read.cigartuples[0][0] == 4 or  # soft clipped at the start
                        read.cigartuples[-1][0] == 4  # soft clipped at the end
                ):

            seq = read.seq
            qual = read.qual

            if read.is_reverse:
                seq = "".join(complement.get(base, base) for base in reversed(seq))
                qual = qual[::-1]

            with open(ofastq, "a") as out_handle:
                out_handle.write("@{}\n{}\n+\n{}\n".format(read.query_name, seq, qual))
        else:
            out_aln.write(read)


def main():

    # Parse the arguments of the script
    parser = argparse.ArgumentParser(description='Collect clipped reads and remap them as split reads')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='../../data/test.bam',
                        help="Specify input file (BAM)")
    parser.add_argument('-ob',
                        '--outbam',
                        type=str,
                        default='test_no_clipped.bam',
                        help="Specify BAM output without clipped reads")
    parser.add_argument('-of',
                        '--outfastq',
                        type=str,
                        default='clipped.fastq',
                        help="Specify FASTQ output with clipped reads")
    parser.add_argument('-l',
                        '--logfile',
                        default='split_clipped_reads.log',
                        help='Log file')

    args = parser.parse_args()

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=args.logfile,
                        filemode='w',
                        level=logging.INFO)

    t0 = time()

    extract_clipped_reads(ibam=args.bam,
                          obam=args.outbam,
                          ofastq=args.outfastq
                          )

    logging.info('Time: extract clipped reads on BAM %s: %f' % (args.bam,
                                                               (time() - t0)))


if __name__ == '__main__':
    main()
