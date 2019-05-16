# Create the intersection between the Mills BEDPE call set and the deletions (DEL) of the nanoSV VCF file
# A match is found if the Mills confidence intervals for start and end of a deletion fully contain the
# confidence interval for the start and end of one or more nanoSV deletion respectively

from pysam import VariantFile
from intervaltree import IntervalTree
from collections import defaultdict
import os


def load_vcf_as_interval_tree(nanosv_vcf):
    vcf_in = VariantFile(nanosv_vcf)

    trees = defaultdict(dict)
    trees['pos'] = defaultdict(IntervalTree)
    trees['end'] = defaultdict(IntervalTree)

    for record in vcf_in.fetch():
        if record.info['SVTYPE'] == 'DEL':
            if 'CIPOS' in record.info and 'CIEND' in record.info:
                bp_id = '_'.join([record.chrom, str(record.pos),
                                  str(record.stop)])

                bp1_left = record.pos + record.info['CIPOS'][0]
                bp1_right = record.pos + record.info['CIPOS'][1] + 1
                trees['pos'][record.chrom][bp1_left:bp1_right] = bp_id

                bp2_left = record.stop + record.info['CIEND'][0]
                bp2_right = record.stop + record.info['CIEND'][1] + 1
                trees['end'][record.chrom][bp2_left:bp2_right] = bp_id

    return trees


mills_bedpe = '/Users/lsantuari/Documents/External_GitHub/' + \
              'sv_benchmark/input.na12878/lumpy-Mills2012-call-set.chr.bedpe'
context = 'trio/NA12878'
nanosv_vcf = '/Users/lsantuari/Documents/Data/germline/' + context + \
             '/SV/Filtered/last_nanosv.sorted.sym.vcf'

mills_nanosv = '/Users/lsantuari/Documents/External_GitHub/' + \
               'sv_benchmark/input.na12878/lumpy-Mills2011-call-set.nanosv.bedpe'

nanosv_vcf_in = VariantFile(nanosv_vcf)

nanosv_trees = load_vcf_as_interval_tree(nanosv_vcf)
# print(nanosv_trees)

output_lines = []

for line in open(mills_bedpe):

    chrom, bp1_begin, bp1_end, chrom2, bp2_begin, bp2_end = line.rstrip().split("\t")[:6]

    bp1_begin = int(bp1_begin)
    bp1_end = int(bp1_end)

    bp2_begin = int(bp2_begin)
    bp2_end = int(bp2_end)

    assert chrom == chrom2

    # print('%s-%s:%s-%s' % (bp1_begin, bp1_end, bp2_begin, bp2_end))
    bp1_match = nanosv_trees['pos'][chrom][bp1_begin:bp1_end]
    bp2_match = nanosv_trees['end'][chrom][bp2_begin:bp2_end]

    set_bp1 = set([l for s, e, l in bp1_match if s >= bp1_begin and e <= bp1_end])
    len_set_bp1 = len(set_bp1)

    set_bp2 = set([l for s, e, l in bp2_match if s >= bp2_begin and e <= bp2_end])
    len_set_bp2 = len(set_bp2)

    if len(set_bp1 & set_bp2) == len_set_bp1 and len_set_bp1 == len_set_bp2 and len_set_bp1 > 0:

        # print('BP1 => %s' % bp1_match)
        # print('BP2 => %s' % bp2_match)

        output_lines.append(line)

print('%d entries' % len(output_lines))

with open(mills_nanosv, "w") as text_file:
    for line in output_lines:
        chrom, bp1_begin, bp1_end, chrom2, bp2_begin, bp2_end = line.rstrip().split("\t")[:6]
        #if chrom in ['1', '2', '3']:
        text_file.write(line)

filename, file_extension = os.path.splitext(mills_nanosv)
# requires SURVIVOR and bedtools in PATH
os.system(' '.join(["sortBed", "-i", filename + ".bedpe", ">", filename + ".sorted.bedpe"]))
os.system(' '.join(["SURVIVOR", "bedpetovcf", filename+".sorted.bedpe", filename+".sorted.vcf"]))