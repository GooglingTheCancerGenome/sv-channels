import pysam
import numpy as np
from collections import defaultdict


class Segment(object):
    def __init__(self):
        self.segment = pysam.AlignedSegment()

    def is_left_clipped(self):
        """Check if the read is clipped on the left.

        :param read:
            read object of the class `pysam.AlignedSegment`
        :return:
            True if the read is soft (4) or hard (5) clipped on the left,
            False otherwise
        """
        if self.segment.cigartuples is not None:
            if self.segment.cigartuples[0][0] in [4, 5]:
                return True
        return False

    def is_right_clipped(self):
        """Check if the read is clipped on the right.

        :param read:
            read object of the class `pysam.AlignedSegment`
        :return:
            True if the read is soft (4) or hard (5) clipped on the right,
            False otherwise
        """
        if self.segment.cigartuples is not None:
            if self.segment.cigartuples[-1][0] in [4, 5]:
                return True
        return False

    def is_clipped(self):
        """Check if the read is clipped on the left or right."""
        if self.is_left_clipped is True or self.is_right_clipped is True:
            return True
        return False


class Alignment(object):
    def __init__(self, bam_file):
        self.alignment = pysam.AlignmentFile(bam_file, "rb")

    def get_length(self):
        """Get chromosome lenghts from the BAM header."""
        chr_len = {}
        for d in self.alignment.header["SQ"]:
            chr_len[d['SN']] = d['LN']
        return chr_len

    def get_coverage(self):
        """Compute sequence coverage per position."""
        for chr, len in self.get_length().items():
            cov = np.zeros(len, dtype=int)
            for pile in self.alignment.pileup(chr, 0, len, truncate=True):
                cov[pile.pos] = pile.n
            yield (chr, cov)

    def get_clipped_reads(self, chrom, min_mapq=30):
        """Count clipped reads per position."""
        clipped_reads = dict()  # left- or right-clipped reads per position
        for direction in ["L", "R"]:
            clipped_reads[direction] = defaultdict(int)

        # Dictionary to store number of clipped reads per position for
        # INVersion:
        # reads that are clipped AND mapped on the same chromosome AND with
        # the same orientation (FF or RR)
        # Two channels: mate is mapped before or after the read
        clipped_reads_inversion = dict()

        # DUPlication:
        # 1) reads that are right-clipped AND mapped on the same chromosome
        # AND read is forward AND mate is reverse AND mate is mapped before read
        # 2) reads that are left-clipped AND mapped on the same chromosome
        # AND read is reverse AND mate is forward AND mate is mapped after read
        clipped_reads_duplication = dict()

        # Mate is mapped (B)efore or (A)fter?
        for mate_pos in ["B", "A"]:
            clipped_reads_inversion[mate_pos] = defaultdict(int)
            clipped_reads_duplication[mate_pos] = defaultdict(int)

        # Fetch the reads mapped on the chromosome
        start, end = 0, self.get_length(chrom)
        for read in self.alignment.fetch(chrom, start, end):
            print(type(read))
            return
            # both read and mate should be mapped
            # if read.is_unmapped is False and \
            #    read.mate_is_unmapped is False and \
            #    read.mapping_quality >= min_mapq:
            #     if read.is_left_clipped is True:
            #         ref_pos = read.reference_start
            #         clipped_reads["L"][ref_pos] += 1
            #     if read.is_right_clipped is True:
            #         ref_pos = read.reference_end + 1
            #         clipped_reads["R"][ref_pos] += 1
            #     # INVersion
            #     if read.is_clipped is True and \
            #        read.reference_name == read.next_reference_name:
            #         if read.is_left_clipped is True:
            #             ref_pos = read.reference_start
            #             # DUPlication, channel 2
            #             if read.is_reverse is True and \
            #                read.mate_is_reverse is False and \
            #                read.reference_start < read.next_reference_start:
            #                 clipped_reads_duplication["A"][ref_pos] += 1
            #         elif read.is_right_clipped is True:
            #             ref_pos = read.reference_end + 1
            #             # DUPlication, channel 1
            #             if read.is_reverse is False and \
            #                read.mate_is_reverse is True and \
            #                read.reference_start > read.next_reference_start:
            #                 clipped_reads_duplication["A"][ref_pos] += 1
            #         if read.is_reverse == read.mate_is_reverse:
            #             if read.reference_start > read.next_reference_start:
            #                 clipped_reads_inversion["B"][ref_pos] += 1
            #             else:
            #                 clipped_reads_inversion["A"][ref_pos] += 1

        return (clipped_reads, clipped_reads_inversion, clipped_reads_duplication)
