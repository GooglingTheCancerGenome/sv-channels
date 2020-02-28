import numpy as np

from pysam import AlignmentFile

class Channel(object):
    def __init__(self, seqid, len, dtype=np.uint):
        self._seqid = seqid
        self._vec = np.zeros(len, dtype)

    def set_value(self, idx, value):
        self._vec[idx] = value

    def get_vec(self):
        return self._vec

    def get_value(self, idx):
        return self._vec[idx]

    def get_seqid(self):
        return self._seqid

    def save(self, outfile):
        np.save(outfile, self.get_vec())

class Alignment(AlignmentFile):
    def __init__(self, bam_file):
        super().__init__()

    def get_seqlen(self):
        """Get sequence lenghts from BAM header."""
        seqlen = dict()
        for d in self.header['SQ']:
            seqlen[d['SN']] = d['LN']
        return seqlen

    def get_seqcov(self):
        """Compute sequence coverage per position."""
        for seqid, len in self.get_seqlen().items():
            ch = Channel(seqid, len)
            for pile in self.pileup(seqid, 0, len):
                ch.set_value(pile.pos, pile.n)
            yield ch

    def get_seqcov2(self, min_mapq=10):
        """Compute sequence coverage per position."""
        for seqid, len in self.get_seqlen().items():
            ch = Channel(seqid, len)
            for read in self.fetch(seqid, 0, len):
                if read.reference_name == read.next_reference_name and \
                    read.is_reverse != read.mate_is_reverse and \
                    read.mapping_quality >= min_mapq:
                    idx = slice(read.reference_start, read.reference_end - 1)
                    ch.set_value(idx, ch.get_value(idx) + 1)
            yield ch
