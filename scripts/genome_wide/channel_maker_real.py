import pysam
import bz2
import os
import cPickle as pickle
from collections import Counter

def main():

    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    inputBAM = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"

    clipped_pos_wd = "/Users/lsantuari/Documents/Data/GiaB/NA12878_NA24385_Tumor_like/ChannelMaker/T/clipped_read_pos/"

    bamfile = pysam.AlignmentFile(inputBAM, "rb")
    header_dict = bamfile.header

    chr_list = [i['SN'] for i in header_dict['SQ'] if i['SN']]

    for chr in chr_list:

        outFile = clipped_pos_wd + chr + '_clipped_reads_pos.pbz2'

        if os.path.isfile(outFile):

            print('Chromosome %s' % chr)
            with bz2.BZ2File(outFile, 'rb') as f:
                clipped_pos = pickle.load(f)

            cpos_cnt = Counter(clipped_pos.values())
            #cpos_cnt = Counter(Counter(clipped_pos).values())

            for i in range(2, 5):
                print('Number of positions with at least %d clipped reads: %d' %
                      (i+1, sum([v for k, v in cpos_cnt.items() if k > i])))


if __name__ == '__main__':
    main()