from label_classes import *
from time import time
import argparse
import logging
from collections import Counter, defaultdict


def create_labels(sampleName, ibam):

    def make_tree_from_vcf(sv_list):

        t = IntervalTree()

        for var in sv_list:

            if type(var) == StructuralVariant:

                id_start = var.svtype + '_start'
                id_end = var.svtype + '_end'

                assert var.start <= var.end, "Start: "+str(var.start)+" End: "+str(var.end)

                t[var.start + var.cipos[0]:var.start + var.cipos[1] + 1] = id_start
                t[var.end + var.ciend[0]:var.end + var.ciend[1] + 1] = id_end

            elif type(var) == tuple:

                start, end, lab = var
                assert start <= end, "Start: " + str(start) + " End: " + str(end)
                t[start:end + 1] = lab

        return t

    def get_bed_sv():

        inbed_path = hpc_path if HPC_MODE else \
            '/Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/'
        inbed = os.path.join(inbed_path, 'Mills2011_nanosv_full_inclusion.unique.bed')

        bed_sv = read_bed_sv(inbed)

        inbed_path = '/Users/lsantuari/Documents/Processed/NA12878/Long_read_validation/'
        inbed = os.path.join(inbed_path, 'lumpy-Mills2012-DEL.bedpe')
        bedpe_sv = read_bedpe_sv(inbed)

        total_sv = 0
        for k in bedpe_sv.keys():
            total_sv += len(bedpe_sv[k])
        logging.info('Total bedpe: %d' % total_sv)

        sv_ci = defaultdict(list)
        # logging.info(bed_sv.keys())
        for chrom in bed_sv.keys():

            bp1_start_list = [start for start, end, lab in bed_sv[chrom] if lab == 'DEL_start']
            # logging.info(bp1_start_list)
            bp1_end_list = [end for start, end, lab in bed_sv[chrom] if lab == 'DEL_start']
            # logging.info(bp1_end_list)

            bp2_start_list = [start for start, end, lab in bed_sv[chrom] if lab == 'DEL_end']
            # logging.info(bp2_start_list)
            bp2_end_list = [end for start, end, lab in bed_sv[chrom] if lab == 'DEL_end']
            # logging.info(bp2_end_list)

            for s in bedpe_sv[chrom]:
                # logging.info(s)
                bp1_start, bp1_end, bp1_lab, bp2_start, bp2_end, bp2_lab = s
                if bp1_start in bp1_start_list and bp1_end in bp1_end_list and \
                        bp2_start in bp2_start_list and bp2_end in bp2_end_list:
                    sv_ci[chrom].append((bp1_start, bp1_end, bp1_lab, bp2_start, bp2_end, bp2_lab))
                else:
                    pass
                    # logging.info('Not found!')

        return sv_ci

    def get_sv_dict():

        hpc_path = os.path.join('/hpc/cog_bioinf/ridder/users/lsantuari/Processed/Data_for_labels/', sampleName)

        sv_dict = dict()
        # sv_dict['nanosv'] = read_nanosv_vcf(sampleName)
        sv_dict['Mills2011_nanosv'] = get_bed_sv()

        total_sv = 0
        for k in sv_dict['Mills2011_nanosv'].keys():
            total_sv += len(sv_dict['Mills2011_nanosv'][k])
            # logging.info('Mills2011_nanosv %s length: %d'% (k, len(sv_dict['Mills2011_nanosv'][k])))
        logging.info('Total bedpe with unique ci: %d' % total_sv)

        return sv_dict

    sv_dict = get_sv_dict()
    labels = dict()
    labels['id'] = defaultdict(list)

    chr_len = get_chr_len_dict(ibam)
    candidate_pairs = load_candidate_pairs(sampleName, ibam)

    candidate_pairs_chr = dict()

    for chrom in chr_len.keys():

        candidate_pairs_chr[chrom] = [sv for sv in candidate_pairs[chrom]
                                      if sv.tuple[0].chr == sv.tuple[1].chr
                                      and sv.tuple[0].chr == chrom]
        logging.info('Chr%s: %d pairs, %d on single chr' %
              (chrom, len(candidate_pairs[chrom]), len(candidate_pairs_chr[chrom])))

        for pair in candidate_pairs_chr[chrom]:
            pair_id = pair.id()
            labels['id'][chrom].append(pair_id)

    for sv_dict_key in sv_dict.keys():

        logging.info(f'running {sv_dict_key}')

        labels[sv_dict_key] = defaultdict(list)

        sv_list = sv_dict[sv_dict_key]

        if type(sv_list) is list:
            logging.info('VCF mode')
            # Select deletions (DELs)
            logging.info('%d SVs (all)' % len(sv_list))

            sv_list = [sv for sv in sv_list if sv.svtype == 'DEL']
            logging.info('%d SVs' % len(sv_list))

            # list of chromosomes
            chr_list = set([var.chrom for var in sv_list])

        else:
            logging.info('BED mode')
            chr_list = sv_list.keys()

        for chrName in chr_list:

            logging.info(f'running Chr{chrName}')

            # Load CR positions, once
            pairs = candidate_pairs_chr[chrName]
            # logging.info('%d/%d' % (len(set([sv.id() for sv in pairs])), len(pairs)))

            # VCF file SVs
            if type(sv_list) is list:

                sv_list_chr = [var for var in sv_list if var.chrom == chrName and var.chrom == var.chrom2]
                # tree_vcf = make_tree_from_vcf(sv_list_chr)

                partial_overlap, full_overlap = get_pairs_with_overlap(sv_list_chr, pairs, 'VCF')

                # logging.info(full_overlap)

            else:

                sv_list_chr = sv_list[chrName]

                partial_overlap, full_overlap = get_pairs_with_overlap(sv_list_chr, pairs, 'BED')
                # logging.info(full_overlap)
                for pair in candidate_pairs_chr[chrName]:
                    pair_id = pair.id()
                    if pair_id in full_overlap:
                        labels[sv_dict_key][chrName].append('DEL')
                    elif pair_id in partial_overlap:
                        labels[sv_dict_key][chrName].append('UK')
                    else:
                        labels[sv_dict_key][chrName].append('noSV')
                logging.info(Counter(labels[sv_dict_key][chrName]))

                # logging.info(full_overlap)

    if not HPC_MODE:
        channel_dir = '/Users/lsantuari/Documents/Data/HPC/DeepSV/GroundTruth'
    else:
        channel_dir = '.'

    output_dir = '/'.join((channel_dir, sampleName, 'label_pairs'))
    create_dir(output_dir)

    data_file = '/'.join((output_dir, 'labels.pickle'))
    # logging.info(output_dir)
    pickle.dump(labels, open(data_file, "wb"))
    os.system('gzip ' + data_file)


def count_pairs(sampleName, ibam):

    chr_len = get_chr_len_dict(ibam)

    candidate_pairs = load_candidate_pairs(sampleName, ibam)

    candidate_pairs_chr = dict()

    chr_list = list(map(str, np.arange(4, 23)))
    chr_list.append('X')

    total_pairs = 0

    for chrom in chr_list:

        candidate_pairs_chr[chrom] = [sv for sv in candidate_pairs[chrom]
                                      if sv.tuple[0].chr == sv.tuple[1].chr
                                      and sv.tuple[0].chr == chrom]

        total_pairs += len(candidate_pairs_chr[chrom])

    print('Number of pairs for chrlist: %d' % total_pairs)


def main():

    '''
    Main function for parsing the input arguments and calling the channel_maker function
    :return: None
    '''

    # Default BAM file for testing
    # On the HPC
    # wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/' + \
    # 'DeepSV/artificial_data/run_test_INDEL/samples/T0/BAM/T0/mapping'
    # inputBAM = wd + "T0_dedup.bam"
    # Locally
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/run_test_INDEL/BAM/'
    inputBAM = wd + "T1_dedup.bam"

    parser = argparse.ArgumentParser(description='Create channels from saved data')
    parser.add_argument('-b', '--bam', type=str,
                        default=inputBAM,
                        help="Specify input file (BAM)")
    parser.add_argument('-o', '--out', type=str, default='labels.npy.gz',
                        help="Specify output")
    parser.add_argument('-l', '--logfile', type=str, default='labels.log',
                        help="Specify log file")
    parser.add_argument('-s', '--sample', type=str, default='NA12878',
                        help="Specify sample")

    args = parser.parse_args()

    logfilename = args.logfile
    FORMAT = '%(asctime)s %(message)s'
    logging.basicConfig(
        format=FORMAT,
        filename=logfilename,
        filemode='w',
        level=logging.INFO)

    t0 = time()

    #for sampleName in ['NA12878', 'Patient1', 'Patient2']:
    for sampleName in ['NA12878']:
        create_labels(sampleName, inputBAM)
        # count_pairs(sampleName, inputBAM)

    logging.info('Elapsed time making labels = %f' % (time() - t0))


if __name__ == '__main__':
    main()
