#
# This script takes FASTA/BED/BAM/SAM files as input and generates N feature vectors.
#
# Dependencies: sambamba or samtools
#
import re  # regex
from numpy import *
import numpy as np
import os
import sys

from ChannelMakerNew_lib import *
from subprocess import call

# argparser look up for flag args for future pipeline

''' >>> Some parameters <<< '''
window_to_each_side = 100  # 5 #10 #100 #500                           #SHOUD BE 100
print 'window_to_each_side:', window_to_each_side
full_window_length = window_to_each_side * 2  # +1                  #SHOULD BE 200
print 'full_window_length:', full_window_length
# max_coverage_depth_from_30x = 80 #1000 #200 #hardcoded number based on 6*30; expect at most 100 reads at one location; matrix below is just for testing really.
# '''
# For future use, need to have max_coverage_depth_from_30x as a non-hardcoded number
# read in max coverage dynamically with sambamba depth for a given window, max depth for window
# Perhaps use: "sambamba depth region -L chr:start-end bamfile"
# '''
# print 'max_coverage_depth_from_30x:', max_coverage_depth_from_30x
illumina_rd_length_fixed = 150
print 'illumina_rd_length_fixed:', illumina_rd_length_fixed
chromosome_number = 17  # Perhaps un-hardcode in future
print chromosome_number
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

Gclass = {'GS': GS_bam_30xsubsample_file, 'G1': G1_bam_30x_file}
Sclass = {'S2': S2_bam_file, 'S3N': S3N_bam_file}
''' NOSVclass = {'':,'':} '''
Gclass_keys = Gclass.keys()
Sclass_keys = Sclass.keys()
'''NOSVclass_keys = NOSVclass.keys()'''


def main():
    t0 = time()
    counter = 0
    ''' THESE TRUTH COORDINATES ARE ONLY FOR G1,GS AND S2,S3N FILES AND NOT YET FOR NOSV1,NOSV2 FILES '''
    '''start_end_SV_DEL = locations_DEL_INS(Truth_set_file)[0]'''
    start_end_SV_DEL, start_end_SV_DEL_INS = locations_DEL_INS(Truth_set_file)
    # print '##DEBUG',start_end_SV_DEL, start_end_SV_DEL_INS
    '''*******CHANGE TO THIS: for coord in start_end_SV_DEL: change coord dependency in name after have tried making 10 windows! change window length to each side to 100*******'''
    for outcoord in start_end_SV_DEL:  # [0:100] #changed 'i' to 'coord'
        print 'coord:', outcoord

        clipped_pos_list = [ el
                             for sample in Gclass.values() + Sclass.values()
                             for el in get_clipped_positions(sample, str(chromosome_number),
                                                  outcoord - window_to_each_side, outcoord + window_to_each_side - 1 )
                              ]
        for coord in sorted(list(set((clipped_pos_list)))):
            print 'inner coord:', coord

            window_arange, left, right = make_window(coord, window_to_each_side)
            current_genome_reference = current_reference_fcn(current_line_on_reference, left, right)
            # print 'current_genome_reference:', current_genome_reference
            GC_content = GC_content_dict_fcn(current_genome_reference)
            # print 'GC_content:', GC_content

            '''The _G indicates that have a germline pure plus 50 percent mix germline and somatic and these files start with a 'G', these are used to form the somatic category '''
            vstack_12_channels_pair_Gclass_list = []
            for element1 in Gclass.keys():  # should be 2 elements
                '''counter = 0 #counts number of windows #reset counter for every new break point
                print 'counter:', counter
                if counter == 1: #then do 10 #then do all #so 10 main windows only w/o embedded windows for the moment
                    break
                '''
                # print 'element1:', element1
                # print 'Gclass[element1]:', Gclass[element1]
                all_reads_in_window_file_name_G = all_reads_in_window(Gclass[element1], chromosome_number, left, right,
                                                                      element1, 0)
                # print 'all_reads_in_window_file_name_G:', all_reads_in_window_file_name_G
                number_of_reads_in_window_total_G = number_of_reads_in_window_compute(all_reads_in_window_file_name_G)
                # print 'number_of_reads_in_window_total_G:', number_of_reads_in_window_total_G
                matrix_str_updated_G, matrix_int_left_updated_G, matrix_int_right_updated_G = matrix_read_updater_for_str_int(
                    all_reads_in_window_file_name_G, coord, window_to_each_side, number_of_reads_in_window_total_G,
                    full_window_length)
                vstack_12_channels_G = channels_12_vstacker(matrix_str_updated_G, matrix_int_left_updated_G,
                                                            matrix_int_right_updated_G, current_genome_reference)
                vstack_12_channels_pair_Gclass_list.append(vstack_12_channels_G)
                counter += 1

            # print 'vstack_12_channels_pair_Gclass_list:', vstack_12_channels_pair_Gclass_list
            # print 'len(vstack_12_channels_pair_Gclass_list):', len(vstack_12_channels_pair_Gclass_list)
            vstack_12_channel_pairer_plus_GC_chanel_G = vstack_12_channel_pairer_plus_GC_chanel_fcn(
                vstack_12_channels_pair_Gclass_list[0], vstack_12_channels_pair_Gclass_list[1], GC_content)
            # print 'vstack_12_channel_pairer_plus_GC_chanel_G:', vstack_12_channel_pairer_plus_GC_chanel_G
            channel_matrix_list_Somatic_categ.append(vstack_12_channel_pairer_plus_GC_chanel_G)
            ##############print 'channel_matrix_list_Somatic_categ:', channel_matrix_list_Somatic_categ
            label_list_Somatic.append('somatic')
            '''print 'counter == 9630:', counter == 9630 #TESTS THAT ARE USING all DELs that have simulated only and not the INS!'''

            ''' The _S here indicates that have somatic files starting with 'S' but these are used to form the germline category!'''
            vstack_12_channels_pair_Sclass_list = []  # should be 2 elements
            for element2 in Sclass.keys():
                all_reads_in_window_file_name_S = all_reads_in_window(Sclass[element2], chromosome_number, left, right,
                                                                      element2, 0)
                number_of_reads_in_window_total_S = number_of_reads_in_window_compute(all_reads_in_window_file_name_S)
                matrix_str_updated_S, matrix_int_left_updated_S, matrix_int_right_updated_S = matrix_read_updater_for_str_int(
                    all_reads_in_window_file_name_S, coord, window_to_each_side, number_of_reads_in_window_total_S,
                    full_window_length)
                vstack_12_channels_S = channels_12_vstacker(matrix_str_updated_S, matrix_int_left_updated_S,
                                                            matrix_int_right_updated_S, current_genome_reference)
                vstack_12_channels_pair_Sclass_list.append(vstack_12_channels_S)
            # print 'vstack_12_channels_pair_Sclass_list:', vstack_12_channels_pair_Sclass_list
            # print 'len(vstack_12_channels_pair_Sclass_list):', len(vstack_12_channels_pair_Sclass_list)
            vstack_12_channel_pairer_plus_GC_chanel_S = vstack_12_channel_pairer_plus_GC_chanel_fcn(
                vstack_12_channels_pair_Sclass_list[0], vstack_12_channels_pair_Sclass_list[1], GC_content)
            # print 'vstack_12_channel_pairer_plus_GC_chanel_S:', vstack_12_channel_pairer_plus_GC_chanel_S
            channel_matrix_list_Germline_categ.append(vstack_12_channel_pairer_plus_GC_chanel_S)
            label_list_Germline.append('germline')

    np.save('somatic_cube_data_file', channel_matrix_list_Somatic_categ)
    np.save('germline_cube_data_file', channel_matrix_list_Germline_categ)
    np.save('somatic_label_array_file', label_list_Somatic)
    np.save('germline_label_array_file', label_list_Germline)

    print 'done in ', time() - t0
    print 'counter/2:', counter / 2
    print 'counter/2 == 9630:', counter / 2 == 9630


if __name__ == '__main__':
    main()
