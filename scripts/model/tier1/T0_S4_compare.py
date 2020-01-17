import argparse
import sys, os

import numpy as np
import pandas as pd

# add scripts folder in path
sys.path.append(
    os.path.dirname(
        os.path.dirname(
            os.path.dirname(
                os.path.realpath(__file__)
            )
        )
    )
)
from genome_wide.label_classes import SVRecord

from T0_S1_generate_training_data import get_truth_set_trees


def read_vcf(invcf):
    # Check file existence
    assert os.path.isfile(invcf), invcf + ' not found!'
    # Dictionary with chromosome keys to store SVs
    sv_list = []

    vcf_in = VariantFile(invcf, 'r')
    for rec in vcf_in.fetch():

        var = SVRecord(rec, 'gridss', 'NA24385')

        chrom1 = var.chrom
        pos1_start = var.start + var.cipos[0]
        pos1_end = var.start + var.cipos[1] + 1

        chrom2 = var.chrom2
        pos2_start = var.end + var.ciend[0]
        pos2_end = var.end + var.ciend[1] + 1
        svtype = var.svtype

        if svtype == "DEL":
            sv_list.append((
                chrom1, pos1_start, pos1_end,
                chrom2, pos2_start, pos2_end,
                svtype
            ))

    print('{} SVs'.format(len(sv_list)))

    return sv_list


def compare(args):

    def write_bed(hit_ids):

        lines = []
        half_interval = 1000/2

        for i in hit_ids:
            c1, s1_0, s1_1, c2, s2_0, s2_1 = i.split('_')
            s1_0 = int(s1_0)
            s2_0 = int(s2_0)
            lines.append('\t'.join([str(c1), str(int(s1_0) - half_interval), str(s1_0 + half_interval)]) + '\n')
            lines.append('\t'.join([str(c2), str(s2_0 - half_interval), str(s2_0 + half_interval)]) + '\n')

        f = open(args.outputbed, 'w')
        try:
            # use set to make lines unique
            for l in lines:
                f.write(l)
        finally:
            f.close()

    def load_predicted_positions(chrlist, inputdir):

        bp1_list = {}
        bp2_list = {}

        for c in chrlist:

            infile = os.path.join(inputdir, '_'.join([c, 'predictions.npz']))
            assert os.path.exists(infile), infile + ' not found'

            with np.load(infile) as data:
                bp1_list[c] = data['start']
                bp2_list[c] = data['end']

        return bp1_list, bp2_list

    def get_hit_ids(bp_hits):

        bp_hits_idlist = set()

        for c in bp_hits.keys():
            for e in bp_hits[c]:
                if len(e) > 0:
                    for sub_e in e:
                        s, e, d = sub_e
                        bp_hits_idlist.add(d)
        return bp_hits_idlist

    # load truth set
    sv_list, sv_bp1_tree, sv_bp2_tree = get_truth_set_trees(args.truthset)

    # window half length
    win_hlen = int(args.win / 2)

    hit_ids_bp1 = set()
    hit_ids_bp2 = set()

    for inputdir in args.inputdirlist:

        # load results of the genome scanning
        bp1_list, bp2_list = load_predicted_positions(args.chrlist, inputdir)

        bp1_hits = {}
        bp2_hits = {}

        # Look for matches
        for c in args.chrlist:
            if c in sv_bp1_tree.keys() and c in bp1_list.keys():
                bp1_hits[c] = [sv_bp1_tree[c][p - win_hlen:p + win_hlen] for p in bp1_list[c]]
            if c in sv_bp2_tree.keys() and c in bp2_list.keys():
                bp2_hits[c] = [sv_bp2_tree[c][p - win_hlen:p + win_hlen] for p in bp2_list[c]]

        hit_ids_bp1 = hit_ids_bp1 | get_hit_ids(bp1_hits)
        hit_ids_bp2 = hit_ids_bp2 | get_hit_ids(bp2_hits)

    # writing BED output
    write_bed(hit_ids_bp1 | hit_ids_bp2)

    # Collect stats

    hit_ids_bp1_len = len(hit_ids_bp1)
    hit_ids_bp2_len = len(hit_ids_bp2)
    bp1_bp2_union_len = len(hit_ids_bp1 | hit_ids_bp2)

    n_svs = len([l for l in sv_list if l[0] in args.chrlist])

    # write CSV output file
    df = pd.DataFrame({'bp1_num': [hit_ids_bp1_len],
                       'bp1_sv_cov': [hit_ids_bp1_len / n_svs * 100],
                       'bp2_num': [hit_ids_bp2_len],
                       'bp2_sv_cov': [hit_ids_bp2_len / n_svs * 100],
                       'bp1_bp2__sv_cov': [bp1_bp2_union_len / n_svs * 100],
                       })
    df.to_csv(path_or_buf=args.output, index=False)


def main():

    truth_sets = {
        'NA24385': os.path.join(
            '/Users/lsantuari/Documents/Data/germline/NA24385/NIST_SVs_Integration_v0.6/processed',
            'HG002_SVs_Tier1_v0.6.PASS.vcf.gz')
    }

    parser = argparse.ArgumentParser(
        description='Train model',
        usage='''T0_S2_training.py [<args>]
            ''')
    parser.add_argument('-truthset', type=str,
                        default=truth_sets['NA24385'],
                        help="Positive set file(s)")
    parser.add_argument('-chrlist', nargs='+', default=['17'],
                        help="List of chromosomes to consider")
    parser.add_argument('-win', type=int, default=200,
                        help="Window length")
    parser.add_argument('-inputdirlist', nargs='+',
                        default=[''],
                        help="Positive set file(s)")
    parser.add_argument('-output', type=str,
                        default='results.csv',
                        help="Specify output")
    parser.add_argument('-outputbed', type=str,
                        default='regions_of_interest.bed',
                        help="Specify output file for regions of interest for GRIDSS targeted approach")

    args = parser.parse_args()
    compare(args)


if __name__ == '__main__':
    main()
