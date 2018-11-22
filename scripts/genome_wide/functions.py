# Imports
import numpy as np
import twobitreader as twobit

'''
Generic functions used in the channel scripts
'''

# Return if a read is clipped on the left
def is_left_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the left, False otherwise
    '''
    if read.cigartuples is not None:
        if read.cigartuples[0][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right
def is_right_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the right, False otherwise
    '''
    if read.cigartuples is not None:
        if read.cigartuples[-1][0] in [4, 5]:
            return True
    return False


# Return if a read is clipped on the right or on the left
def is_clipped(read):
    '''

    :param read: read object of the class pysam.AlignedSegment
    :return: True if the read is soft (4) or hard (5) clipped on the left or on the right, False otherwise
    '''
    if read.cigartuples is not None:
        if is_left_clipped(read) or is_right_clipped(read):
            return True
    return False


# Return the mate of a read. Get read mate from BAM file
def get_read_mate(read, bamfile):
    '''
    This function was used in a previous version of the code when we needed to check if the mate was clipped.
    Retrieving the mate for each read is time consuming and should be avoided when not strictly necessary.

    :param read: object of the class pysam.AlignedSegment whose mate needs to be retrieved
    :param bamfile: BAM file containing both the read and its mate
    :return: mate, object of the class pysam.AlignedSegment, if a mate for the read is found. Return None otherwise.
    '''
    # print(read)
    # The mate is located at:
    # chromosome: read.next_reference_name
    # positions: [read.next_reference_start, read.next_reference_start+1]
    # Fetch all the reads in that location and retrieve the mate
    iter = bamfile.fetch(read.next_reference_name,
                         read.next_reference_start, read.next_reference_start + 1,
                         multiple_iterators=True)
    for mate in iter:
        # A read and its mate have the same query_name
        if mate.query_name == read.query_name:
            # Check if read is first in pair (read1) and mate is second in pair (read2) or viceversa
            if (read.is_read1 and mate.is_read2) or (read.is_read2 and mate.is_read1):
                # print('Mate is: ' + str(mate))
                return mate
    return None


#Return a one-hot encoding for the chromosome region chr:start-stop
# with Ns encoded as 1 and other chromosomes encoded as 0
def get_one_hot_sequence(chrname, start, stop, HPC_MODE):

    if HPC_MODE:
        # Path on the HPC of the 2bit version of the human reference genome (hg19)
        genome = twobit.TwoBitFile('/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/genomes/hg19.2bit')
    else:
        # Path on the local machine of the 2bit version of the human reference genome (hg19)
        genome = twobit.TwoBitFile('/Users/lsantuari/Documents/Data/GiaB/reference/hg19.2bit')

    ltrdict = {'a': 1, 'c': 2, 'g': 3, 't': 4, 'n': 0}
    # N one-hot
    #ltrdict = {'a': 0, 'c': 0, 'g': 0, 't': 0, 'n': 1}
    return np.array([ltrdict[x.lower()] for x in genome['chr'+chrname][start:stop]])


# From https://github.com/joferkington/oost_paper_code/blob/master/utilities.py
def is_outlier(points, thresh=3.5):
    """
    Returns a boolean array with True if points are outliers and False
    otherwise.
    Parameters:
    -----------
        points : An numobservations by numdimensions array of observations
        thresh : The modified z-score to use as a threshold. Observations with
            a modified z-score (based on the median absolute deviation) greater
            than this value will be classified as outliers.
    Returns:
    --------
        mask : A numobservations-length boolean array.
    References:
    ----------
        Boris Iglewicz and David Hoaglin (1993), "Volume 16: How to Detect and
        Handle Outliers", The ASQC Basic References in Quality Control:
        Statistical Techniques, Edward F. Mykytka, Ph.D., Editor.
    """
    if len(points.shape) == 1:
        points = points[:,None]
    median = np.median(points, axis=0)
    diff = np.sum((points - median)**2, axis=-1)
    diff = np.sqrt(diff)
    med_abs_deviation = np.median(diff)

    modified_z_score = 0.6745 * diff / med_abs_deviation

    return modified_z_score > thresh