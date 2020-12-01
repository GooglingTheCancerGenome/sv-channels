import re

from pysam import VariantRecord

__bpRE__ = None
__symbolicRE__ = None


def setupREs():
    '''
    Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
    URL: https://github.com/ljdursi/mergevcf
    '''
    global __symbolicRE__
    global __bpRE__
    if __symbolicRE__ is None or __bpRE__ is None:
        __symbolicRE__ = re.compile(r'.*<([A-Z:]+)>.*')
        __bpRE__ = re.compile(
            r'([ACGTNactgn\.]*)([\[\]])([a-zA-Z0-9\.\_]+:\d+)([\[\]])([ACGTNacgtn\.]*)'
        )


class SVRecord:
    def __init__(self, record, svcaller):
        ci_slop = 0
        # For CHM[1|13] SVs
        svtype_dict = {
            'deletion': 'DEL',
            'insertion': 'INS',
            'inversion': 'INV'
        }

        if isinstance(record) != VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')

        if record.info['SVTYPE'] == 'BND':
            ct, chr2, pos2, indellen = self.get_bnd_info(record)
        else:
            ct = None
            chr2 = record.chrom
            pos2 = record.stop
            if 'SVLEN' in record.info.keys():
                self.indellen = record.info['SVLEN']
            else:
                self.indellen = abs(record.stop - record.pos)

        self.id = record.id
        self.chrom = record.chrom.replace('chr', '')
        self.start = record.pos
        self.chrom2 = chr2.replace('chr', '')
        self.end = pos2
        self.alt = record.alts[0]
        self.insLen = len(self.alt)
        self.svLen = self.end - self.start + 1 if self.chrom == self.chrom2 else None
        self.ct = ct

        # CIPOS
        if 'CIPOS' in record.info.keys():
            if 'CIPOS95' in record.info.keys():
                self.cipos = record.info['CIPOS95']
            else:
                self.cipos = record.info['CIPOS']
        else:
            self.cipos = (-ci_slop, ci_slop)

        # CIEND
        if 'CIEND' in record.info.keys():
            if 'CIEND95' in record.info.keys():
                self.ciend = record.info['CIEND95']
            else:
                self.ciend = record.info['CIEND']
        elif 'CIRPOS' in record.info.keys():
            self.ciend = record.info['CIRPOS']
        else:
            self.ciend = (-ci_slop, ci_slop)
        self.filter = record.filter

        # set SVTYPE
        if svcaller is None:
            self.svtype = record.info['SVTYPE']
        else:
            if self.chrom != self.chrom2:
                self.svtype = 'BP'
            elif self.insLen >= abs(self.svLen) * 0.7:
                self.svtype = 'INS'
            elif self.ct in ['5to5', '3to3']:
                self.svtype = 'INV'
            elif (self.start < self.end) != (ct in ['5to3', '3to5']):
                self.svtype = 'DEL'
            else:
                self.svtype = record.info['SVTYPE']

    @staticmethod
    def stdchrom(chrom):
        if chrom[0] == 'c':
            return chrom[3:]
        return chrom

    def locFromBkpt(self, ref, pre, delim1, pair, delim2, post):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''
        chpos = pair.split(':')
        chr2 = self.stdchrom(chpos[0])
        pos2 = int(chpos[1])
        assert delim1 == delim2  # '['..'[' or ']'...']'
        joinedAfter = True
        extendRight = True
        connectSeq = ""

        if len(pre) > 0:
            connectSeq = pre
            joinedAfter = True
            assert len(post) == 0
        elif len(post) > 0:
            connectSeq = post
            joinedAfter = False

        if delim1 == "]":
            extendRight = False
        else:
            extendRight = True
        indellen = len(connectSeq) - len(ref)

        if joinedAfter:
            if extendRight:
                ct = '3to5'
            else:
                ct = '3to3'
        else:
            if extendRight:
                ct = '5to5'
            else:
                ct = '5to3'
        return ct, chr2, pos2, indellen

    def get_bnd_info(self, record):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''
        setupREs()
        altstr = str(record.alts[0])
        resultBP = re.match(__bpRE__, altstr)
        if resultBP:
            ct, chr2, pos2, indellen = self.locFromBkpt(
                str(record.ref), resultBP.group(1), resultBP.group(2),
                resultBP.group(3), resultBP.group(4), resultBP.group(5))
            res = (ct, chr2, pos2, indellen)
        else:
            res = None
        return res
