#
# This script takes FASTA/BED/BAM/SAM files as input and generates N feature vectors.
#
# Dependencies: sambamba or samtools
#
import re #regex
from numpy import *
import numpy as np
import os
import sys

from ChannelMakerNew_lib import *
from subprocess import call

#argparser look up for flag args for future pipeline

''' >>> Some parameters <<< '''
window_to_each_side = 5 #10 #100 #500                           #SHOUD BE 100
print 'window_to_each_side:', window_to_each_side
full_window_length = window_to_each_side*2 #+1                  #SHOULD BE 200
print 'full_window_length:', full_window_length
#max_coverage_depth_from_30x = 80 #1000 #200 #hardcoded number based on 6*30; expect at most 100 reads at one location; matrix below is just for testing really.
#'''
#For future use, need to have max_coverage_depth_from_30x as a non-hardcoded number
#read in max coverage dynamically with sambamba depth for a given window, max depth for window
#Perhaps use: "sambamba depth region -L chr:start-end bamfile"
#'''
#print 'max_coverage_depth_from_30x:', max_coverage_depth_from_30x
illumina_rd_length_fixed = 150
print 'illumina_rd_length_fixed:', illumina_rd_length_fixed
chromosome_number = 17 #Perhaps un-hardcode in future
print chromosome_number
''' >>> Some parameters <<<'''

#data_chr17_fasta=open("chr17_human.fasta",'r')
data_chr17_fasta = open(data_chr17_fasta_file,'r')
current_line_on_reference = one_liner(data_chr17_fasta)
data_chr17_fasta.close()

#Truth_set_file = "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed" #contains DELs and INs 4.815 DELs and 5.185 INS, INS we don't consider for the moment
#lines in this file look like:  17	982	17	1910	INS

'''on local machine'''
cmd = "sambamba"

''' #uncomment for use on hpc
#cmd = "module load sambamcram/sambamba"
cmd_hpc = "module load sambamcram/sambamba" #on the HPC cluster #sambamcram/samtools"
call(cmd_hpc, shell=True)
'''

#from time import time
#from subprocess import call

''' Remember that need to combine channel info into 3 pairs! '''
''' Make these also into list comprehensions '''
channel_matrix_list_Germline = []
channel_matrix_list_Somatic = []
channel_matrix_list_NoSV = []

def main():
    ''' THESE TRUTH COORDINATES ARE ONLY FOR G1,GS AND S2,S3N FILES AND NOT YET FOR NOSV1,NOSV2 FILES '''
    start_end_SV_DEL = locations_DEL_INS(Truth_set_file)[0]
    #######start_end_SV_DEL, start_end_SV_DEL_INS = locations_DEL_INS(Truth_set_file)
    #print '##DEBUG',start_end_SV_DEL, start_end_SV_DEL_INS
    for coord in start_end_SV_DEL: #changed 'i' to 'coord'
        counter = 0 #counts number of windows #reset counter for every new break point
        window_arange, left, right = make_window(i,window_to_each_side)
        current_genome_reference = current_reference_fcn(current_line_on_reference,left,right)
        GC_content = GC_content_dict(current_genome_reference)
        names = ['G1','GS','S2','S3N','NoSV1','NoSV2'] #the three pairs of files that will have
        ''' Need to convert what's below this line into a meta-function to pair GS with G1 and S3N with S2 and NoSV1 with NoSV2 '''
        ''' Make another loop here where will go through the G1, GS, S2, S3N file names: The NoSV category will be treated differently because will zoom in on clipped regions specifically so will have different truth set file for those '''
        ''' indent; for k in [G1_sam_30x_file, GS_sam_30xsubsample_file, S2_sam_file, S3N_sam_file] and then shift all lines underneath '''
        all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,chromosome_number,left,right,names[1],counter)
        ''' in place of GS_bam_30xsubsample_file in the above line do %s with '%(k)' substitution from above list, or maybe seperate function for this '''
        #all_reads_in_window_file_name_updatedforchr17only = NEW_FCN
        number_of_reads_in_window_total = number_of_reads_in_window_compute(all_reads_in_window_file_name) #_updatedforchr17only
        #matrix_str, matrix_int_left_clip, matrix_int_right_clip, matrix_int = make_matrix(number_of_reads_in_window,full_window_length)
        matrix_str, matrix_int_left_clip, matrix_int_right_clip = matrix_read_updater_for_str_int(all_reads_in_window_file_name,coord,window_to_each_side,number_of_reads_in_window_total,full_window_length)

        #####with ( open(all_reads_in_window_file_name,'r') ) as j:
            #####rd_counter = 0
            #####for line in j:
                #rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line)
                #read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
                #read_cleaned_length = len(read_cleaned)
                #read_cleaned_length_computed_end = rd_pos_start + read_cleaned_length

                #read_cleaned_length_arange = arange(rd_pos_start,read_cleaned_length_computed_end)
                #read_content_in_window = intersect1d(window_arange,read_cleaned_length_arange)
                #####rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line)

                #####read_content_in_window = read_content_in_window_fcn(line,i,window_to_each_side)

                #####if len(read_content_in_window) == 0: #''' If read present inside the window! ''' All reads should be present by definition of fetching the reads!
                    #####raise ValueError('primary read outside window')

                #begin_of_intersection = read_content_in_window[0]
                #end_of_intersection = read_content_in_window[len(read_content_in_window)-1]

                #window_begin_intersect_index = where(window_arange == begin_of_intersection)[0][0]
                #window_end_intersect_index = where(window_arange == end_of_intersection)[0][0]
                #rd_begin_intersect_index = where(rd_sequence_ref_mapped_arange == begin_of_intersection)[0][0]
                #rd_end_intersect_index = where(rd_sequence_ref_mapped_arange == end_of_intersection)[0][0]

                #####window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index = clean_read_mapper(line,i,window_to_each_side,rd_sequence_raw,rd_cigar)

                #####matrix_str[counter][window_begin_intersect_index:window_end_intersect_index+1] = list(read_cleaned)[rd_begin_intersect_index:rd_end_intersect_index+1]

                #####break

                #####''' use matrix_int here'''
                #####letter_order = re.findall(r"[SH]",rd_cigar) #to find whether Soft or Hard clipped at beginning or end of the read
                #####if rd_pos_start >= begin_of_intersection and (letter_order[0] == 'S' or letter_order[0] == 'H'): #read already present in window
                #####    ind_rd_start = where(rd_pos_start == window_arange)[0][0]
                #####    matrix_GS_left_clip[counter][ind_rd_start] = 1
                #####if end_of_intersection >= rd_pos_computed_end and (letter_order[len(letter_order)-1] == 'S' or letter_order[len(letter_order)-1] == 'H'): #read already present in window
                #####    ind_rd_end = where(rd_pos_computed_end == window_arange)[0][0]
                #####    matrix_GS_right_clip[counter][ind_rd_end] = 1


                #####rd_counter += 1

        """
        first stack 12 x 2 times and then vstack 'GC_channel = GC_content_dict(entry)' to this 24 block!
        vstack((12channels,12channels,GC_channel)) #check if need dtype=int
        Germline_channels = channels_25_vstacker(,...,)
        Somatic_channels = channels_25_vstacker(,...,)
        NoSV_channels = channels_25_vstacker(,...,)

        channel_matrix_list_Germline.append(Germline_channels)
        channel_matrix_list_Somatic.append(Somatic_channels)
        channel_matrix_list_NoSV.append(NoSV_channels)
        """

        if counter == 1: #then do 10 #then do all #so 10 main windows only w/o embedded windows for the moment
            break

        counter += 1
    print 'counter == 9630:', counter == 9630

''' label object for 3 label cubes '''
''' save as .npy both 3 channel objects and label objects: np.save('outfile1',data_cube), np.save('outfile2', data_array) '''


        #for j in range(number_of_reads_in_window):
        #    read_info_extractor(all_reads_in_window_file_name)



        #all_reads_in_window_file_data = open(all_reads_in_window_file_name,'r')
        #all_reads_in_window_file.close()



        #t0 = time()
        #print '##DEBUG ', 'running sambamba'
        #all_reads_in_window()

        #string1="sambamba view %s %d:%d-%d | cut -f1 >> idFile_GS_%s.txt" %(GS_bam_30xsubsample_file,chromosome_number,left,right,counter) #should be 17 #chr_list[counter] for more dynamic calling
        #call(string1, shell=True)
        #print 'done in ',time()-t0
        #with(open(,'r')) as i:
        #    for line in i:
        #        line = line.rstrip()
        #        string2="sambamba view -F \"read_name == %s\" %s > reads_with_same_name.sam" %(rd_name,GS_bam_30xsubsample_file)

        '''
        t1 = time()
        print '##DEBUG ', 'running grep'
        string2="LC_ALL=C grep -w -F -f idFile_GS_%s.txt < %s >> subset_GS_%s.sam" %(counter,GS_sam_30xsubsample_file,counter) #can convert sam to bam and view with IGV
        call(string2, shell=True)
        print 'done in ',time()-t1
        break
        #string2="grep -f idFile_GS.txt %s >> output_GS.sam" %(GS_sam_30xsubsample_file)
        t2 = time()
        print '##DEBUG ', 'running cat'
        string3="cat subset_GS_%s.sam | grep -v '#'| sort -n -k4,1 > sorted_subset_GS_%s.sam" %(counter,counter)
        call(string3, shell=True)
        print 'done in ',time()-t2
        #os.system(cmd+"\n"+string1+"\n"+string2+"\n"+string2p)
        '''
        #os.system(cmd+"\n"+string1+"\n"+string2)
        #break




''' List to append the vstacked channels '''
''' also include labels '''
''' save as .npy both objects '''
''' array(vstack((histogram_channel_GS,coverage_channel_GS,clipped_rds_left_GS,clipped_rds_right_GS,...histogram_channel_G1,coverage_channel_G1,clipped_rds_left_G1,clipped_rds_right_G1,GC_content)),dtype=int) '''



if __name__ == '__main__':
    main()
