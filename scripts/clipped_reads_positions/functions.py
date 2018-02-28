import pysam
import bz2
import logging
import sys
import cPickle as pickle
from collections import Counter

# DEBUG printing
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT, stream=sys.stderr, level=logging.DEBUG)

# Global variables
window_size = 200
clipped_read_support = 3
clipped_pos_wd = "/Users/lsantuari/Documents/Data/GiaB/NA12878_NA24385_Tumor_like/ChannelMaker/T/clipped_read_pos/"

# Return if a read is clipped on the left
def is_left_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    if read.cigartuples is not None:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right or on the left
def is_clipped(read):
    if read.cigartuples is not None:
        if is_left_clipped(read) or is_right_clipped(read):
            return True
    return False


# Return the mate of a read. Get read mate from BAM file
def get_read_mate(read, bamfile):
    # print(read)
    # print('Mate is mapped')
    iter = bamfile.fetch(read.next_reference_name,
                         read.next_reference_start, read.next_reference_start + 1,
                         multiple_iterators=True)
    for mate in iter:
        if mate.query_name == read.query_name:
            if (read.is_read1 and mate.is_read2) or (read.is_read2 and mate.is_read1):
                # print('Mate is: ' + str(mate))
                return mate
    return None


# Return clipped reads positions for a specific chromosome
def get_clipped_reads_positions(chr_name, support=clipped_read_support):

    outFile = clipped_pos_wd + chr_name + '_clipped_reads_pos.pbz2'

    with bz2.BZ2File(outFile, 'rb') as f:
        clipped_pos = pickle.load(f)

    clipped_pos_list = [k for k, v in Counter(clipped_pos).items() if v > support]
    return sorted(clipped_pos_list)
