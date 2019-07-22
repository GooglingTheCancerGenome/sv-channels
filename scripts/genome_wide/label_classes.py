from pysam import VariantRecord

# Classes

class SVRecord:

    def __init__(self, record):

        ci_slop = 0

        # For CHM[1|13] SVs
        svtype_dict = {'deletion': 'DEL', 'insertion': 'INS', 'inversion': 'INV'}

        if type(record) != pysam.VariantRecord:
            raise TypeError('VCF record is not of type pysam.VariantRecord')
        # logging.info(record)

        if record.info['SVTYPE'] == 'BND':
            ct, chr2, pos2, indellen = self.get_bnd_info(record)
        else:
            # logging.info(record.info['SVTYPE'])
            ct = None
            chr2 = record.chrom
            pos2 = record.stop
            if 'SVLEN' in record.info.keys():
                indellen = record.info['SVLEN']
            else:
                indellen = abs(record.stop - record.pos)

        # logging.info(record.info.keys())

        self.id = record.id
        self.chrom = record.chrom.replace('chr', '')
        self.start = record.pos
        self.chrom2 = chr2.replace('chr', '')
        self.end = pos2
        self.alt = record.alts[0]

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
        self.svtype = record.info['SVTYPE']

    @staticmethod
    def stdchrom(chrom):

        if chrom[0] == 'c':
            return chrom[3:]
        else:
            return chrom

        # Modified from the function ctAndLocFromBkpt in mergevcf

    def locFromBkpt(self, ref, pre, delim1, pair, delim2, post):
        '''
        Function of the mergevcf tool by Jonathan Dursi (Simpson Lab)
        URL: https://github.com/ljdursi/mergevcf
        :param record: pysam.VariantRecord
        :return: tuple with connection (3' to 5', 3' to 3', 5' to 5' or 5' to 3'), chromosome and position of the
        second SV endpoint, length of the indel
        '''

        chpos = pair.split(':')
        # logging.info(chpos[0])
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
            ct, chr2, pos2, indellen = self.locFromBkpt(str(record.ref), resultBP.group(1),
                                                        resultBP.group(2), resultBP.group(3), resultBP.group(4),
                                                        resultBP.group(5))
        return (ct, chr2, pos2, indellen)
