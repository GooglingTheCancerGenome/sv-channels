import pysam
import argparse
import logging
from time import time


def slice_bam(ibam, ibed, outFile):
    '''
    This function reduces a BAM file to a set of regions, centered on positions specified in the BED file,
    for inspection.
    It takes in input a BAM alignment file and a BED file with a set of positions. For each position, the
    reads overlapping the win_len*2 interval centered on the position are fetched from the BAM file and saved
    in the output BAM file.
    :param ibam: BAM file used to extract slice from
    :param ibed: BED file of intervals to consider for the slices
    :param outFile: BAM file in output
    :return: None
    '''

    win_len = 100
    bamfile = pysam.AlignmentFile(ibam, "rb")
    header_dict = bamfile.header
    bam_outfile = pysam.AlignmentFile(outFile, "wb", template=bamfile)

    with open(ibed, 'r') as bedfile:
        for row in bedfile:
            row_list = row.strip().split('\t')
            chr, start, end = (row_list[0], int(row_list[1]), int(row_list[2]))
            chr_len = [i['LN'] for i in header_dict['SQ'] if i['SN'] == chr][0]
            for pos in [start, end]:
                for read in bamfile.fetch(chr, max([0, pos - win_len]), min([pos + win_len, chr_len])):
                    bam_outfile.write(read)

    bam_outfile.close()
    bamfile.close()
    # pysam.sort(outFile)


def main():

    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.bam"

    inputBED = wd + "chr17_somaticallele_10k_INS_DEL.bed"

    parser = argparse.ArgumentParser(description='Create channels with split read distance for left/right split reads')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-be', '--bed', type=str,
                        default=inputBED,
                        help="Specify input file (BAM)")
    parser.add_argument('-o', '--out', type=str, default='slice.bam',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', default='slice_bam.log',
                        help='Specify log file')

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        level=logging.INFO)

    t0 = time()
    slice_bam(ibam=args.bam, ibed=args.bed, outFile=args.out)
    print(time() - t0)


if __name__ == '__main__':
    main()