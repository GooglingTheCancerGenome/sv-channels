import argparse
import gzip
import json
import logging
import os
import sys

from collections import Counter, defaultdict
from time import time

import pysam
from intervaltree import IntervalTree

from functions import *
from label_classes import SVRecord


def read_vcf(invcf):

    # Check file existence
    assert os.path.isfile(invcf), invcf + ' not found!'
    # Dictionary with chromosome keys to store SVs
    sv_list = []

    vcf_in = pysam.VariantFile(invcf, 'r')
    for rec in vcf_in.fetch():

        var = SVRecord(rec, 'gridss')

        chrom1 = var.chrom
        pos1_start = var.start + var.cipos[0]
        pos1_end = var.start + var.cipos[1] + 1

        chrom2 = var.chrom2
        pos2_start = var.end + var.ciend[0]
        pos2_end = var.end + var.ciend[1] + 1
        svtype = var.svtype

        sv_list.append((chrom1, pos1_start, pos1_end, chrom2, pos2_start,
                        pos2_end, svtype))

    logging.info('{} SVs'.format(len(sv_list)))

    return sv_list


def read_bedpe(inbedpe, svtype_to_select):
    # Check file existence
    assert os.path.isfile(inbedpe), inbedpe + ' not found!'
    # Dictionary with chromosome keys to store SVs
    sv_list = []

    with (open(inbedpe, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])
            svtype = columns[-1]

            if svtype == "TYPE:DELETION":
                svtype = "DEL"

            if svtype_to_select == svtype:
                if svtype in ['DEL', 'INV', 'DUP', 'CTX']:
                    sv_list.append((chrom1, pos1_start, pos1_end, chrom2,
                                    pos2_start, pos2_end, svtype))
                elif svtype == "INS":
                    sv_list.append((chrom1, pos1_start, pos1_end, chrom1,
                                    pos1_start + 1, pos1_end + 1, svtype))

    logging.info('{} SVs'.format(len(sv_list)))

    return sv_list


def filter_bedpe(inbedpe, sv_id_list, outDir):

    # Check file existence
    assert os.path.isfile(inbedpe), inbedpe + ' not found!'
    # Dictionary with chromosome keys to store SVs
    logging.info('{} SVs to filter out'.format(len(sv_id_list)))
    lines_to_keep = []

    with (open(inbedpe, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])
            svtype = columns[-1]
            svtype = "DEL" if svtype == "TYPE:DELETION" else svtype

            sv_id = '_'.join(
                (svtype, chrom1, str(pos1_start), chrom2, str(pos2_start)))

            if svtype in ['DEL', 'INS', 'INV', 'DUP', 'CTX'] and sv_id not in sv_id_list:
                lines_to_keep.append(line)

        fileout = os.path.join(outDir, 'uncaptured_SVs.bedpe')
        logging.info('Writing {}'.format(fileout))

        with (open(fileout, 'w')) as fout:
            for line in lines_to_keep:
                fout.write(line)

        logging.info('{} SVs written'.format(len(lines_to_keep)))


def read_svcaller_bedpe(inbedpe):

    # Check file existence
    assert os.path.isfile(inbedpe), inbedpe + ' not found!'
    # Dictionary with chromosome keys to store SVs
    cr_pos = []

    with (open(inbedpe, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(columns[0]), int(
                columns[1]), int(columns[2])
            chrom2, pos2_start, pos2_end = str(columns[3]), int(
                columns[4]), int(columns[5])

            cr_pos.append((chrom1, pos1_start, chrom2, pos2_start))

    logging.info('{} candidate positions'.format(len(cr_pos)))

    return cr_pos


def make_gtrees_from_svlist(sv_list):

    logging.info('Building SV GenomicTrees...')
    # Tree with windows for candidate positions
    trees_start = defaultdict(IntervalTree)
    trees_end = defaultdict(IntervalTree)

    # Populate tree
    for sv in sv_list:
        chrom1, pos1_start, pos1_end, chrom2, pos2_start, pos2_end, svtype = sv
        sv_id = '_'.join(
            (svtype, chrom1, str(pos1_start), chrom2, str(pos2_start)))

        trees_start[chrom1][pos1_start:pos1_end] = (svtype, sv_id)
        trees_end[chrom2][pos2_start:pos2_end] = (svtype, sv_id)

    # print('Tree start')
    # for k in trees_start.keys():
    #     print('{} : {}'.format( k, len(trees_start[k])))
    # print('Tree end')
    # for k in trees_end.keys():
    #     print('{} : {}'.format( k, len(trees_end[k])))

    return trees_start, trees_end


def search_tree_with_cpos(cpos, trees_start, trees_end, win_hlen):

    logging.info('Searching SV GenomicTrees with candidate positions...')

    lookup_start = []
    lookup_end = []

    # Log info every n_r times
    n_r = 10**6
    last_t = time()

    for i, p in enumerate(cpos, start=1):

        if not i % n_r:
            now_t = time()
            # print(type(now_t))
            logging.info(
                "%d candidate positions processed (%f positions / s)" %
                (i, n_r / (now_t - last_t)))
            last_t = time()

        chrom1, pos1, chrom2, pos2, strand_info = p
        lookup_start.append(trees_start[chrom1].envelop(pos1 - win_hlen,
                                                        pos1 + win_hlen + 1))
        lookup_end.append(trees_end[chrom2].envelop(pos2 - win_hlen,
                                                    pos2 + win_hlen + 1))

    return lookup_start, lookup_end


def overlap(svtype, sv_list, cpos_list, win_hlen, ground_truth, outDir):
    '''

    :param sv_list: list, list of SVs
    :param cr_pos: list, list of clipped read positions
    :return: list, list of clipped read positions whose window completely overlap either the CIPOS interval
    or the CIEND interval
    '''

    trees_start, trees_end = make_gtrees_from_svlist(sv_list)
    lookup_start, lookup_end = search_tree_with_cpos(cpos_list, trees_start,
                                                     trees_end, win_hlen)

    # print([l for l in lookup_start if len(l) > 0])
    # print([l for l in lookup_end if len(l) > 0])

    labels = dict()

    sv_covered = set()

    for p, lu_start, lu_end in zip(cpos_list, lookup_start, lookup_end):

        chrom1, pos1, chrom2, pos2, strand_info = p
        pos_id = '_'.join((chrom1, str(pos1), chrom2, str(pos2), strand_info))

        l1 = len(lu_start)
        l2 = len(lu_end)

        if l1 == 1 and l1 == l2:

            # print(lu_start)
            # print(lu_end)
            lu_start_elem_start, lu_start_elem_end, lu_start_elem_data = lu_start.pop(
            )
            lu_end_elem_start, lu_end_elem_end, lu_end_elem_data = lu_end.pop()

            lu_start_elem_svtype, lu_start_elem_svid = lu_start_elem_data
            lu_end_elem_svtype, lu_end_elem_svid = lu_end_elem_data

            # if lu_start_elem_svtype != 'TRA':

            if pos1 - win_hlen <= lu_start_elem_start and lu_start_elem_end <= pos1 + win_hlen and \
                    pos2 - win_hlen <= lu_end_elem_start and lu_end_elem_end <= pos2 + win_hlen and \
                    lu_start_elem_svid == lu_end_elem_svid:
                # logging.info(
                #     'Chr1:{}\tpos1:{}-{}\tChr2:{}\tpos2:{}-{}'.format(
                #         chrom1, pos1 - win_hlen, pos1 + win_hlen, chrom2, pos2 - win_hlen, pos2 + win_hlen
                #     )
                # )
                # logging.info(
                #     'LookUp_start:{}-{}_{}\tLookUp_end:{}-{}_{}'.format(
                #         lu_start_elem_start, lu_start_elem_end, lu_start_elem_data,
                #         lu_end_elem_start, lu_end_elem_end, lu_end_elem_data
                #     )
                # )
                # if pos1 in np.arange(lu_start_elem_start-2, lu_start_elem_end+2) and \
                #         pos2 in np.arange(lu_end_elem_start-2, lu_end_elem_end+2):
                sv_covered.add(lu_start_elem_svid)
                labels[pos_id] = lu_start_elem_svtype
                # else:
                #     sv_covered.add(lu_start_elem_svid)
                #     labels[pos_id] = 'UK_overlap_not_matching'
            else:
                labels[pos_id] = 'UK_single_partial'

            # else:
            #
            #     sv_covered.add(lu_start_elem_svid)
            #     labels[pos_id] = lu_start_elem_svtype

        elif l1 > 1 or l2 > 1:

            lu_start_set = set()
            lu_end_set = set()

            for s in lu_start:
                lu_start_elem_start, lu_start_elem_end, lu_start_elem_data = s
                lu_start_elem_svtype, lu_start_elem_svid = lu_start_elem_data
                lu_start_set.add(lu_start_elem_svid)

            for s in lu_end:
                lu_end_elem_start, lu_end_elem_end, lu_end_elem_data = s
                lu_end_elem_svtype, lu_end_elem_svid = lu_end_elem_data
                lu_end_set.add(lu_end_elem_svid)

            sv_covered = sv_covered | (lu_start_set & lu_end_set)

            # if svtype in ['DEL', 'INV', 'DUP', 'TRA']:
            #     labels[pos_id] = 'UK_multiple_on_either_windows'
            # elif svtype == 'INS':
            #     labels[pos_id] = svtype

            labels[pos_id] = svtype

        elif l1 == 0 and l1 == l2:
            # logging.info('CPOS->Partial: %s\t%d\t%d' % (elem, start, end))
            labels[pos_id] = 'no'+svtype

        elif (l1 == 1 and l2 > 1) or (l2 == 1 and l1 > 1) or (l2 > 1 and l1 > 1):

            lu_start_set = set()
            lu_end_set = set()

            for s in lu_start:
                lu_start_elem_start, lu_start_elem_end, lu_start_elem_data = s
                lu_start_elem_svtype, lu_start_elem_svid = lu_start_elem_data
                lu_start_set.add(lu_start_elem_svid)

            for s in lu_end:
                lu_end_elem_start, lu_end_elem_end, lu_end_elem_data = s
                lu_end_elem_svtype, lu_end_elem_svid = lu_end_elem_data
                lu_end_set.add(lu_end_elem_svid)

            sv_covered = sv_covered | (lu_start_set & lu_end_set)
            # labels[pos_id] = 'UK_single_and_multiple'
            labels[pos_id] = svtype

        else:
            # (l1 == 1 and l2 == 0) or (l1 == 0 and l2 == 1)
            labels[pos_id] = 'no'+svtype

    logging.info(Counter(labels.values()))
    sv_coverage = int(len(sv_covered) / len(sv_list) * 100)
    logging.info('SV coverage: {}/{}={}%'.format(len(sv_covered), len(sv_list),
                                                 sv_coverage))

    filename, file_extension = os.path.splitext(ground_truth)

    if file_extension == '.bedpe':
        # print(sv_covered)
        filter_bedpe(ground_truth, sv_covered, outDir)
    elif file_extension == '.sur':
        filter_survivor_output(ground_truth, sv_covered, outDir)

    return labels


# Get labels
def get_labels(chrlist, chr_dict, win_len, svtype, ground_truth, sv_positions,
               channelDataDir, outFile, outDir):

    # windows half length
    win_hlen = int(int(win_len) / 2)

    sv_caller_file = os.path.join(channelDataDir, sv_positions + '.bedpe')
    sv_caller_name = os.path.basename(sv_positions)

    if os.path.exists(sv_caller_file):

        cpos_list = read_svcaller_bedpe(sv_caller_file)

    elif sv_caller_name == 'split_reads':

        cpos_list = load_all_clipped_read_positions(win_hlen, svtype, chr_dict,
                                                    channelDataDir)
    else:
        sys.exit('I cannot find {} nor {}'.format(sv_caller_file, sv_caller_name))

    # Keep only positions that can be used to create windows
    cpos_list = [
        (chrom1, pos1, chrom2, pos2, strand_info)
        for chrom1, pos1, chrom2, pos2, strand_info, in cpos_list if chrom1 in chrlist
        and chrom2 in chrlist and win_hlen <= pos1 <= chr_dict[chrom1] -
        win_hlen and win_hlen <= pos2 <= chr_dict[chrom2] - win_hlen
    ]

    filename, file_extension = os.path.splitext(ground_truth)
    if file_extension == '.bedpe':
        sv_list = read_bedpe(ground_truth, svtype)
    elif file_extension == '.vcf' or file_extension == '.gz':
        sv_list = read_vcf(ground_truth)

    # Get overlap of candidate positions with all SV breakpoints (all 4 SV callers)
    # crpos_all_sv = get_crpos_overlap_with_sv_callsets(sv_dict, cr_pos_dict)

    # filename, file_extension = os.path.splitext(ground_truth)
    # trees_start, trees_end = make_gtrees_from_truth_set(sv_list, file_extension.upper())
    # print(sv_list)
    labels = overlap(svtype, sv_list, cpos_list, win_hlen, ground_truth, outDir)

    with gzip.GzipFile(outFile, 'wb') as fout:
        fout.write(json.dumps(labels).encode('utf-8'))


def main():
    '''
    Label windows according to truth set
    :return: None
    '''

    parser = argparse.ArgumentParser(description='Create labels')
    parser.add_argument('-b',
                        '--bed',
                        type=str,
                        default='../../data/seqs.bed',
                        help="Specify chromosome regions to consider (BED)")
    parser.add_argument('-l',
                        '--logfile',
                        type=str,
                        default='labels.log',
                        help="Specify log file")
    parser.add_argument('-c',
                        '--chrlist',
                        type=str,
                        default='12,22',
                        help="Comma separated list of chromosomes to consider")
    parser.add_argument('-w',
                        '--window',
                        type=str,
                        default=250,
                        help="Specify window size")
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='DEL',
                        help="Specify SV type")
    parser.add_argument('-sv',
                        '--sv_positions',
                        type=str,
                        #default=os.path.join('..', '..', 'data', 'gridss'),
                        default=os.path.join('..', '..', 'data', 'split_reads'),
                        help="Specify Manta/GRIDSS BEDPE file")
    parser.add_argument('-gt',
                        '--ground_truth',
                        type=str,
                        default='../../data/test.bedpe',
                        help="Specify ground truth VCF/BEDPE file")
    parser.add_argument('-o',
                        '--out',
                        type=str,
                        default='labels.json.gz',
                        help="Specify output")
    parser.add_argument('-p',
                        '--outputpath',
                        type=str,
                        default='',
                        help="Specify output path")

    args = parser.parse_args()

    sv_caller_name = os.path.basename(args.sv_positions)
    output_dir = os.path.join(args.outputpath, 'labels',
                              'win' + str(args.window),
                              args.svtype,
                              sv_caller_name
                              )
    os.makedirs(output_dir, exist_ok=True)
    # Log file
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, args.out)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)

    t0 = time()

    # Get dictionary of chromosome lengths
    chr_dict = chr_dict_from_bed(args.bed)

    get_labels(chrlist=args.chrlist.split(','),
               chr_dict=chr_dict,
               win_len=args.window,
               svtype=args.svtype,
               ground_truth=args.ground_truth,
               sv_positions=args.sv_positions,
               channelDataDir=args.outputpath,
               outFile=output_file,
               outDir=output_dir)

    logging.info('Elapsed time making labels = %f' % (time() - t0))


if __name__ == '__main__':
    main()
