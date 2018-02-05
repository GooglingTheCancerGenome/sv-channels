#
# This script takes FASTA/BED/BAM/SAM files as input and generates N feature vectors.
#
# Dependencies: sambamba or samtools
#

from ChannelMakerNew_lib import *

# argparser look up for flag args for future pipeline

''' >>> Some parameters <<< '''

# Is it generating the training set?
TRAINING_MODE = True

window_to_each_side = 100  # 5 #10 #100 #500                           #SHOUD BE 100
logging.debug('window_to_each_side: %d', window_to_each_side)
full_window_length = window_to_each_side * 2  # +1                  #SHOULD BE 200
logging.debug('full_window_length: %d', full_window_length)
# max_coverage_depth_from_30x = 80 #1000 #200 #hardcoded number based on 6*30; expect at most 100 reads at one location; matrix below is just for testing really.
# '''
# For future use, need to have max_coverage_depth_from_30x as a non-hardcoded number
# read in max coverage dynamically with sambamba depth for a given window, max depth for window
# Perhaps use: "sambamba depth region -L chr:start-end bamfile"
# '''
# print 'max_coverage_depth_from_30x:', max_coverage_depth_from_30x
illumina_rd_length_fixed = 150
logging.debug('illumina_rd_length_fixed: %d', illumina_rd_length_fixed)
chromosome_number = 17  # Perhaps un-hardcode in future
logging.debug('Using chromosome %d', chromosome_number)

''' >>> Some parameters <<<'''

# data_chr17_fasta=open("chr17_human.fasta",'r')
data_chr17_fasta = open(data_chr17_fasta_file, 'r')
current_line_on_reference = one_liner(data_chr17_fasta)
data_chr17_fasta.close()

# Truth_set_file = "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed" #contains DELs and INs 4.815 DELs and 5.185 INS, INS we don't consider for the moment
# lines in this file look like:  17	982	17	1910	INS

from time import time

# from subprocess import call


''' Make these also into list comprehensions like the get functions that have '''
channel_matrix_list_Germline_categ = []
channel_matrix_list_Somatic_categ = []
channel_matrix_list_NoSV_categ = []

label_list_Somatic = []
label_list_Germline = []
label_list_NOSV = []

#Is the window centered on the breakpoint junction (BPJ)?
bpj_flag = []

Gclass = {'GS': GS_bam_30xsubsample_file, 'G1': G1_bam_30x_file}
Sclass = {'S2': S2_bam_file, 'S3N': S3N_bam_file}
''' NOSVclass = {'':,'':} '''
Gclass_keys = Gclass.keys()
Sclass_keys = Sclass.keys()
'''NOSVclass_keys = NOSVclass.keys()'''


# Output file with information about the program execution
# info_file = open('ChannelMaker_run_info.txt', 'w')
# info_file.write('\t'.join(['OUTCOORD', 'COORD', 'INDEX']) + '\n')


def get_ch_mtx(coord, bam_class, win_left, win_right, current_genome_reference, GC_content):
    '''
    Returns a channel matrix for a location and for a class of BAM files
    :param coord:
    :param bam_class:
    :param win_left:
    :param win_right:
    :param GC_content:
    :return:
    '''

    vstack_pair = []

    for element in bam_class.keys():
        # print(element)
        all_reads = all_reads_in_window(bam_class[element], chromosome_number,
                                        win_left, win_right,
                                        element, 0)
        n_reads = number_of_reads_in_window_compute(all_reads)
        # print('N reads:%d' % n_reads)
        matrix_str_updated, matrix_int_left_updated, matrix_int_right_updated = matrix_read_updater_for_str_int(
            all_reads, coord, window_to_each_side, n_reads, full_window_length)

        vstack_ch = channels_12_vstacker(matrix_str_updated, matrix_int_left_updated,
                                         matrix_int_right_updated, current_genome_reference)
        vstack_pair.append(vstack_ch)

    vstack_with_GC = vstack_12_channel_pairer_plus_GC_chanel_fcn(
        vstack_pair[0], vstack_pair[1], GC_content)

    return vstack_with_GC


def generate_training_set():
    '''

    :return:
    '''

    t0 = time()
    counter = 0
    # positions with no clipped reads
    no_clip_coord = []
    ''' THESE TRUTH COORDINATES ARE ONLY FOR G1,GS AND S2,S3N FILES AND NOT YET FOR NOSV1,NOSV2 FILES '''
    '''start_end_SV_DEL = locations_DEL_INS(Truth_set_file)[0]'''
    # start_SV_DEL, end_SV_DEL, start_SV_INS = locations_DEL_INS(Truth_set_file)

    start_SV_DEL, end_SV_DEL, start_SV_INS = locations_DEL_INS(Truth_set_file)

    # print '##DEBUG',start_end_SV_DEL, start_end_SV_DEL_INS
    # print(start_end_SV_DEL_INS)

    '''*******CHANGE TO THIS: for coord in start_end_SV_DEL: change coord dependency in name after have tried making 10 windows! change window length to each side to 100*******'''
    for outzipped in zip(start_SV_DEL + end_SV_DEL + start_SV_INS,
                         ['del_start'] * len(start_SV_DEL) +
                         ['del_end'] * len(end_SV_DEL) +
                         ['ins_start'] * len(start_SV_INS)):  # [0:100] #changed 'i' to 'coord'

        outcoord = outzipped[0]
        # print(outzipped[1])

        logging.debug('coord: %d', outcoord)

        clipped_pos_list = [el
                            # for sample in Gclass.values() + Sclass.values()
                            for el in get_clipped_positions(Gclass.values()[0], str(chromosome_number),
                                                            outcoord - window_to_each_side,
                                                            outcoord + window_to_each_side)
                            ]

        # List of clipped positions should have at least the PBJ position
        if len(clipped_pos_list) == 0:
            print outcoord
        # print(clipped_pos_list)

        '''
        for sample in Gclass.values() + Sclass.values():
            get_del_reads_per_pos(sample, str(chromosome_number),
                            outcoord - window_to_each_side,
                            outcoord + window_to_each_side)
        '''

        for coord in sorted(list(set((clipped_pos_list)))):

            logging.debug('inner coord: %d', coord)

            #print(outcoord == coord)
            bpj_flag.append(outcoord == coord)

            # info_file.write('\t'.join([str(outcoord), str(coord), str(counter)]) + '\n')

            window_arange, left, right = make_window(coord, window_to_each_side)
            current_genome_reference = current_reference_fcn(current_line_on_reference, left, right)
            # print 'current_genome_reference:', current_genome_reference
            GC_content = GC_content_dict_fcn(current_genome_reference)
            # print 'GC_content:', GC_content

            '''The _G indicates that have a germline pure plus 50 percent mix germline and somatic and these files start with a 'G', these are used to form the somatic category '''

            channel_matrix_list_Somatic_categ.append(
                get_ch_mtx(coord, Gclass, left, right, current_genome_reference, GC_content)
            )
            label_list_Somatic.append('somatic_' + outzipped[1])

            '''print 'counter == 9630:', counter == 9630 #TESTS THAT ARE USING all DELs that have simulated only and not the INS!'''

            ''' The _S here indicates that have somatic files starting with 'S' but these are used to form the germline category!'''

            channel_matrix_list_Germline_categ.append(
                get_ch_mtx(coord, Sclass, left, right, current_genome_reference, GC_content)
            )
            label_list_Germline.append('germline_' + outzipped[1])

            # for i in range(len(channel_matrix_list_Somatic_categ[0])):
            #    print(channel_matrix_list_Somatic_categ[0][i])

            # Plot channels
            # plot_channels_mtx(channel_matrix_list_Somatic_categ[-1], 'Somatic: outCoord:'+str(outcoord)+
            #                  ' coord:'+str(coord))
            # plot_channels_mtx(channel_matrix_list_Germline_categ[-1], 'Germline: outCoord:'+str(outcoord)+
            #                  ' coord:'+str(coord))

            # print(channel_matrix_list_Germline_categ[0])

            counter += 1

    np.save('somatic_cube_data_file', channel_matrix_list_Somatic_categ)
    np.save('germline_cube_data_file', channel_matrix_list_Germline_categ)
    np.save('somatic_label_array_file', label_list_Somatic)
    np.save('germline_label_array_file', label_list_Germline)
    np.save('BPJ_flag', bpj_flag)

    print 'done in ', time() - t0
    # print 'counter/2:', counter / 2
    # print 'counter/2 == 9630:', counter / 2 == 9630


def bam_to_channels():
    '''
    Converts the Tumor and Normal BAM files into Channel data
    :return:
    '''

    t0 = time()
    counter = 0

    # Input files
    wd = '/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/GiaB/Synthetic_tumor/BAM/'

    # Requires BAM files with Clipped reads only:
    # sambamba view -F "cigar =~ /^\d[S|H]/ or cigar =~ /[S|H]$/" -f bam CPCT11111111T_dedup.realigned.bam
    inbam = wd + 'Tumor/' + 'CPCT11111111T_dedup.realigned.cr.bam'

    tumor_bam = wd + 'Tumor/' + 'CPCT11111111T_dedup.realigned.bam'
    normal_bam = wd + 'Reference/' + 'CPCT11111111R_dedup.realigned.bam'
    tn_dict = {'Tumor': tumor_bam, 'Normal': normal_bam}

    clipped_pos = get_clipped_positions_from_CR_BAM(inbam)

    for chr in clipped_pos.keys():

        channel_matrix_list_TumorNormal_categ = []
        label_list_TumorNormal = []

        for coord in sorted(clipped_pos[chr]):
            logging.debug('coord: %d', coord)

            window_arange, left, right = make_window(coord, window_to_each_side)
            current_genome_reference = current_reference_fcn(current_line_on_reference, left, right)
            # print 'current_genome_reference:', current_genome_reference
            GC_content = GC_content_dict_fcn(current_genome_reference)
            # print 'GC_content:', GC_content

            channel_matrix_list_TumorNormal_categ.append(
                get_ch_mtx(coord, tn_dict, left, right, current_genome_reference, GC_content)
            )
            label_list_TumorNormal.append('tumor_normal')

            counter += 1

        np.save('tumornormal_cube_data_file_chr_' + chr, channel_matrix_list_TumorNormal_categ)
        np.save('tumornormal_label_array_file_chr_' + chr, label_list_TumorNormal)

    print 'done in ', time() - t0
    print 'counter:', counter


def main():
    # load_channels()

    if TRAINING_MODE:
        generate_training_set()
    else:
        bam_to_channels()

    # info_file.close()


if __name__ == '__main__':
    main()
