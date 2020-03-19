# Imports
import argparse
import gzip
import json
import logging
import os
from collections import Counter, defaultdict
from time import time

from intervaltree import IntervalTree
from pysam import VariantFile

from functions import *
from label_classes import SVRecord

with open('parameters.json', 'r') as f:
    config = json.load(f)

HPC_MODE = config["DEFAULT"]["HPC_MODE"]


def get_chr_list():

    chrlist = list(map(str, range(1, 23)))
    chrlist.extend(['X'])
    #chrlist = ['17']

    return chrlist


def read_vcf(invcf):

    # Check file existence
    assert os.path.isfile(invcf), invcf + ' not found!'
    # Dictionary with chromosome keys to store SVs
    sv_list = []

    vcf_in = VariantFile(invcf, 'r')
    for rec in vcf_in.fetch():

        var = SVRecord(rec, 'gridss')

        chrom1 = var.chrom
        pos1_start = var.start + var.cipos[0]
        pos1_end = var.start + var.cipos[1] + 1

        chrom2 = var.chrom2
        pos2_start = var.end + var.ciend[0]
        pos2_end = var.end + var.ciend[1] + 1
        svtype = var.svtype

        if svtype == "DEL":
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
                if svtype == "DEL":
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

            if svtype == "DEL" and sv_id not in sv_id_list:
                lines_to_keep.append(line)

        fileout = os.path.join(outDir, 'uncaptured_SVs.bedpe')
        logging.info('Writing {}'.format(fileout))

        with (open(fileout, 'w')) as fout:
            for line in lines_to_keep:
                fout.write(line)

        logging.info('{} SVs written'.format(len(lines_to_keep)))


def read_survivor_simsv_output(insur, svtype):
    # Check file existence
    assert os.path.isfile(insur), insur + ' not found!'
    # Dictionary with chromosome keys to store SVs
    sv_list = []

    with (open(insur, 'r')) as bed:
        for line in bed:

            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(
                columns[0]), int(columns[1]) - 1, int(columns[1])
            chrom2, pos2_start, pos2_end = str(
                columns[2]), int(columns[3]) - 1, int(columns[3])
            type = columns[4]

            if type == svtype:

                if type == 'DEL':
                    sv_list.append((chrom1, pos1_start, pos1_end, chrom2,
                                    pos2_start, pos2_end, svtype))

                elif type == 'INS':
                    sv_list.append((chrom1, pos1_start, pos1_end, chrom1,
                                    pos1_start, pos1_end, svtype))

    logging.info('{} SVs'.format(len(sv_list)))

    return sv_list


def filter_survivor_output(insur, sv_id_list, outDir):
    # Check file existence
    assert os.path.isfile(insur), insur + ' not found!'

    logging.info('{} SVs to filter out'.format(len(sv_id_list)))
    lines_to_keep = []

    with (open(insur, 'r')) as bed:
        for line in bed:
            columns = line.rstrip().split("\t")
            chrom1, pos1_start, pos1_end = str(
                columns[0]), int(columns[1]) - 1, int(columns[1])
            chrom2, pos2_start, pos2_end = str(
                columns[2]), int(columns[3]) - 1, int(columns[3])
            svtype = columns[4]

            sv_id = '_'.join(
                (svtype, chrom1, str(pos1_start), chrom2, str(pos2_start)))

            if svtype == "DEL" and sv_id not in sv_id_list:
                lines_to_keep.append(line)

    fileout = os.path.join(outDir, 'uncaptured_SVs.sur')
    logging.info('Writing {}'.format(fileout))

    with (open(fileout, 'w')) as fout:
        for line in lines_to_keep:
            fout.write(line)

    logging.info('{} SVs written'.format(len(lines_to_keep)))


def overlap(svtype, sv_list, cpos_list, win_hlen, ground_truth, outDir):
    '''

    :param sv_list: list, list of SVs
    :param cr_pos: list, list of clipped read positions
    :return: list, list of clipped read positions whose window completely overlap either the CIPOS interval
    or the CIEND interval
    '''
    def make_gtrees_from_svlist(sv_list):

        logging.info('Building SV GenomicTrees...')
        # Tree with windows for candidate positions
        trees_start = defaultdict(IntervalTree)
        trees_end = defaultdict(IntervalTree)

        # Populate tree
        for sv in sv_list:
            # print(sv)

            chrom1, pos1_start, pos1_end, chrom2, pos2_start, pos2_end, svtype = sv
            sv_id = '_'.join(
                (svtype, chrom1, str(pos1_start), chrom2, str(pos2_start)))
            # print(sv_id)

            trees_start[chrom1][pos1_start:pos1_end] = (svtype, sv_id)
            trees_end[chrom2][pos2_start:pos2_end] = (svtype, sv_id)

        # print('Tree start')
        # for k in trees_start.keys():
        #     print('{} : {}'.format( k, len(trees_start[k])))
        # print('Tree end')
        # for k in trees_end.keys():
        #     print('{} : {}'.format( k, len(trees_end[k])))

        return trees_start, trees_end

    def search_tree_with_cpos(cpos, trees_start, trees_end):

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

            chrom1, pos1, chrom2, pos2 = p

            lookup_start.append(trees_start[chrom1][pos1 - win_hlen:pos1 +
                                                    win_hlen + 1])
            lookup_end.append(trees_end[chrom2][pos2 - win_hlen:pos2 +
                                                win_hlen + 1])

        return lookup_start, lookup_end

    trees_start, trees_end = make_gtrees_from_svlist(sv_list)
    lookup_start, lookup_end = search_tree_with_cpos(cpos_list, trees_start,
                                                     trees_end)

    #print([l for l in lookup_start if len(l) > 0])
    #print([l for l in lookup_end if len(l) > 0])

    labels = dict()

    sv_covered = set()

    for p, lu_start, lu_end in zip(cpos_list, lookup_start, lookup_end):

        chrom1, pos1, chrom2, pos2 = p
        pos_id = '_'.join((chrom1, str(pos1), chrom2, str(pos2)))

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

            if svtype == 'DEL':
                labels[pos_id] = 'UK_multiple_on_either_windows'
            elif svtype == 'INS':
                labels[pos_id] = svtype

        elif l1 == 0 and l1 == l2:
            # logging.info('CPOS->Partial: %s\t%d\t%d' % (elem, start, end))
            labels[pos_id] = 'no' + svtype

        elif (l1 == 1 and l2 > 1) or (l2 == 1 and l1 > 1):

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
            labels[pos_id] = 'UK_single_and_multiple'

        else:
            # (l1 == 1 and l2 == 0) or (l1 == 0 and l2 == 1)
            labels[pos_id] = 'no' + svtype

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
def get_labels(ibam, chrlist, win_len, svtype, ground_truth, channelDataDir,
               outFile, outDir):

    # windows half length
    win_hlen = int(int(win_len) / 2)
    # get chromosome lengths
    chr_dict = get_chr_len_dict(ibam)

    cpos_list = load_all_clipped_read_positions(win_hlen, svtype, chr_dict,
                                                channelDataDir)

    # Keep only positions that can be used to create windows
    chr_len_dict = get_chr_len_dict(ibam)

    cpos_list = [
        (chrom1, pos1, chrom2, pos2)
        for chrom1, pos1, chrom2, pos2 in cpos_list if chrom1 in chrlist
        and chrom2 in chrlist and win_hlen <= pos1 <= chr_len_dict[chrom1] -
        win_hlen and win_hlen <= pos2 <= chr_len_dict[chrom2] - win_hlen
    ]

    filename, file_extension = os.path.splitext(ground_truth)
    if file_extension == '.bedpe':
        sv_list = read_bedpe(ground_truth, svtype)
    elif file_extension == '.sur':
        sv_list = read_survivor_simsv_output(ground_truth, svtype)
    elif file_extension == '.vcf' or file_extension == '.gz':
        sv_list = read_vcf(ground_truth)

    # Get overlap of candidate positions with all SV breakpoints (all 4 SV callers)
    # crpos_all_sv = get_crpos_overlap_with_sv_callsets(sv_dict, cr_pos_dict)

    # filename, file_extension = os.path.splitext(ground_truth)
    # trees_start, trees_end = make_gtrees_from_truth_set(sv_list, file_extension.upper())

    labels = overlap(svtype, sv_list, cpos_list, win_hlen, ground_truth,
                     outDir)

    with gzip.GzipFile(outFile, 'wb') as fout:
        fout.write(json.dumps(labels).encode('utf-8'))


def main():
    '''
    Label windows according to truth set
    :return: None
    '''

    parser = argparse.ArgumentParser(description='Create labels')
    parser.add_argument('-b',
                        '--bam',
                        type=str,
                        default='',
                        help="Specify input file (BAM)")
    parser.add_argument('-c',
                        '--chrlist',
                        nargs='+',
                        default=['17'],
                        help="List of chromosomes to consider")
    parser.add_argument('-l',
                        '--logfile',
                        type=str,
                        default='labels.log',
                        help="Specify log file")
    # parser.add_argument('-s', '--sample', type=str, default='NA24385',
    #                     help="Specify sample")
    parser.add_argument('-w',
                        '--window',
                        type=str,
                        default=200,
                        help="Specify window size")
    parser.add_argument('-s',
                        '--svtype',
                        type=str,
                        default='INS',
                        help="Specify SV type")
    parser.add_argument('-gt',
                        '--ground_truth',
                        type=str,
                        default='',
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

    # Log file
    output_dir = os.path.join(args.outputpath, 'labels',
                              'win' + str(args.window))
    os.makedirs(output_dir, exist_ok=True)
    logfilename = os.path.join(output_dir, args.logfile)
    output_file = os.path.join(output_dir, args.out)

    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(format=FORMAT,
                        filename=logfilename,
                        filemode='w',
                        level=logging.INFO)

    t0 = time()

    get_labels(
        ibam=args.bam,
        chrlist=args.chrlist,
        win_len=args.window,
        svtype=args.svtype,
        ground_truth=args.ground_truth,
        channelDataDir=args.outputpath,
        outFile=output_file,
        outDir=output_dir,
    )

    logging.info('Elapsed time making labels = %f' % (time() - t0))


if __name__ == '__main__':
    main()
