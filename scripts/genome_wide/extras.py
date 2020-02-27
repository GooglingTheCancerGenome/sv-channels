import numpy as np

from pysam import AlignmentFile

class Alignment(AlignmentFile):
    def __init__(self, bam_file):
        AlignmentFile.__init__(self)

    def get_seqlen(self):
        """Get sequence lenghts from BAM header."""
        seqlen = dict()
        for d in self.header['SQ']:
            seqlen[d['SN']] = d['LN']
        return seqlen

    def get_seqcov(self):
        """Compute sequence coverage per position."""
        for seqid, len in self.get_seqlen().items():
            cov = np.zeros(len, dtype=int)
            for pile in self.pileup(seqid, 0, len):
                cov[pile.pos] = pile.n
            yield (seqid, cov)

    def get_seqcov2(self):
        pass
