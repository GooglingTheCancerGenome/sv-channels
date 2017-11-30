import re #regex
from numpy import *
import numpy as np
import os
import sys
from time import time
from subprocess import call

#sambamba = 'path to sambamba'
#cmd0 = 'module load sambamba'
#cmd = f'{sambamba} view {input} {output}'
#run(cmd0 +'&&' + cmd, shell=True)

def get_cigar_coord(cigar):
    re1='(\\d+)'	# Integer Number 1
    re2='([a-z])'	# Any Single Word Character (Not Whitespace) 1

    rg = re.compile(re1+re2,re.IGNORECASE|re.DOTALL)
    return [(int(i), v) for i,v in re.findall(rg,cigar)]


def add_dels(seq, cigar):
    '''
    Replacing all instances of deletions
    in Cigar string reported as D by Ds
    '''
    i = 0
    for n,v in get_cigar_coord(cigar):
        if v == 'D':
            seq = seq[:i]+'D'*n+seq[i:]
            #print seq
            i += n
        elif v == 'H':
            pass
        else:
            i += n
    return seq


def rem_ins(seq_with_del, cigar):
    '''
    Removes insertions from deletion compensated string
    '''
    seq = seq_with_del
    i = 0
    for n,v in get_cigar_coord(cigar):
        if v == 'I':
            seq = seq[:i] + seq[i+n:]
            #print seq
            i += n
        elif v == 'H': #or v == 'S'
            pass
        else:
            i += n
    return seq


def rem_softclips(seq_with_rem_ins,cigar):
    '''
    Remove soft clipped sequence, can occur at beginning  and/or end of the read
    '''
    seq = seq_with_rem_ins #with removed insertions
    i = 0
    for n,v in get_cigar_coord(cigar):
        if v == 'S':
            seq = list(seq)
            del seq[i:i+n]
            #print seq
            #i += n # i think sequence is now shorter so don't boost by n!
        elif v == 'H' or v == 'I': #ignore values for insertions #can have a hard clip occuring when have a soft clip in seq too!
            pass
        else: #if D or M
            i += n
    return ''.join(seq)


def read_reconstruct(seq, cigar):
    '''
    Remove insertions after add deletions
    '''
    if seq == '*' or cigar == '*':
        return '*'

    else:
        return rem_softclips(rem_ins(add_dels(seq, cigar), cigar), cigar)


def one_liner(data):
    '''
    Currentline is now chr17
    In the .fasta file the sequence is split over lines,
    this combines all lines into one line
    '''
    currentline = ""
    for line in data.xreadlines():
        if line.startswith('>'):
            line = line.rstrip('\n')
            #print line
            currentline = ""
        else:
            line = line.rstrip('\n')
            currentline += line
    #data.close()

    #print 'currentline[0:5]:', currentline[0:5]
    return currentline


def GC_content_dict(entry):

    GC_content_dict = {'G':1,'C':1,'A':0,'T':0}
    GC_contents = [GC_content_dict[k] for k in [i for i in entry]]

    return GC_contents



def locations_DEL_INS(truth_file):
    start_SV_DEL = []
    end_SV_DEL = []
    SV_type_DEL = []
    chr_list_DEL = []

    start_SV_INS = []
    #end_SV_INS = []
    SV_type_INS = []
    chr_list_INS = []
    with(open(truth_file,'r')) as i:
        for line in i:
            line = line.rstrip()
            columns = line.split("\t")
            if columns[4] == "DEL" :
                start_SV_DEL.append(int(columns[1]))
                end_SV_DEL.append(int(columns[3]))
                #SV_type_DEL.append("DEL")
                chr_list_DEL.append(int(columns[0])) #re.findall(r"(\d+)[AB]",columns[0])
                #print 'chromosome:', chromosome
                #chr_list.append(chromosome)
                #start_end_SV_list.append([int(columns[1]),int(columns[4])])
                #break
            elif str(columns[4]) == "INS" :
                start_SV_INS.append(int(columns[1]))
                #end_SV_INS.append(int(columns[4]))
                #SV_type_INS.append("INS")
                chr_list_INS.append(int(columns[0]))

    start_end_SV_DEL = start_SV_DEL + end_SV_DEL
    start_end_SV_DEL_INS = start_end_SV_DEL + start_SV_INS

    return start_end_SV_DEL, start_end_SV_DEL_INS


def make_window(breakpoint,window_to_each_side):
    left = breakpoint - window_to_each_side
    right = breakpoint + window_to_each_side
    window_arange = arange(left,right)

    return window_arange, left, right


def current_reference_fcn(current_line_on_reference,left,right):
    current_ref = current_line_on_reference[left:right]

    #print 'current_ref:', current_ref
    return current_ref


def make_matrix(number_of_reads_in_window_total,full_window_length):
    matrix_str = array(zeros((number_of_reads_in_window_total,full_window_length)),dtype='|S1')
    matrix_int_left_clip = zeros((number_of_reads_in_window_total,full_window_length))
    matrix_int_right_clip = zeros((number_of_reads_in_window_total,full_window_length))
    #matrix_int = zeros((number_of_reads_in_window_total,full_window_length))

    return matrix_str, matrix_int_left_clip, matrix_int_right_clip #,matrix_int


'''
>>>Don't need sorted_read_and_mates_by_start_position() function as this done by default by sambamba view!<<<
def sorted_read_and_mates_by_start_position(all_reads_in_window_file,name,counter): #sorts all reads located in the window #useful for IGV visualization
    all_reads_in_window_file_actual = all_reads_in_window(bam_file,chr_number,left,right,name,counter) #is this correct?
    cmd_sorted_reads_in_window = "cat %s | grep -v '#'| sort -n -k4,1 > sorted_all_reads_in_window_%s_%s.sam" %(all_reads_in_window_file_actual,name,counter) #all_reads_in_window_file_actual
    call(cmd_sorted_reads_in_window, shell=True)
    sorted_all_reads_in_window_file = "sorted_all_reads_in_window_%s_%s.sam" %(name,counter)

    return sorted_all_reads_in_window_file
'''


'''
Remove the _%s.sam (%s = counter) part of the generated file names and just overwrite a single file in actual implementation.
Currenty, > all_reads_in_window_%s_%s.sam is just useful for testing.
'''

'''
sambamba view -F \"ref_name == '17' and mate_ref_name == '17'\" saves the day as 1) don't want any read artificats where the read and its mate don't map to chr17
and 2) was only looking at start and end position and not on chr number, and
3) that this info is otherwise lost in the sam conversion and then there's is no way to know which chrms dealing with.
So now even if the all_reads_in_window_%s_%s_new.sam file has stars these are coming from a bam filtered on ref and mate with chr17 only! so all is fine!
'''
def all_reads_in_window(bam_file,chr_number,left,right,name,counter): #sambamba sorts read by start position automatically
    #''' FOR HPC USE THE FOLLOWING LINE because on laptop have '*' in the RNAME and RNEXT fields which correspond to sambamba's ref_name and mate_ref_name '''
    ''' cigar != '*' and sequence != '*' doesn't work in sambamba but its ok as have taken these exemption into account and will map '*' in my matrix '''
    ''' and cigar != '*' and sequence != '*' >>> removed from the code as output '*' to the matrix and then will avoid it! So all ok! '''
    cmd_all_reads_in_window = "sambamba view -F \"ref_name == '17' and mate_ref_name == '17'\" %s %d:%d-%d > all_reads_in_window_%s_%s_new_new.sam" %(bam_file,chr_number,left,right,name,counter) #name = GS,G1,S3,S1... #_new added 29.11.17
    #''' The line below here have been using on local files on my laptop '''
    #cmd_all_reads_in_window = "sambamba view %s %d:%d-%d > all_reads_in_window_%s_%s.sam" %(bam_file,chr_number,left,right,name,counter) #name = GS,G1,S3,S1...
    call(cmd_all_reads_in_window, shell=True)
    all_reads_in_window_file = "all_reads_in_window_%s_%s_new_new.sam" %(name,counter) #_new added 29.11.17

    return all_reads_in_window_file

def number_of_reads_in_window_compute(sam_file):
    counter_lines = 0
    sample_file_data = open(sam_file,'r')
    for line in sample_file_data.xreadlines():
        counter_lines += 1
    sample_file_data.close()

    return counter_lines

'''
Subracting: -1 from rd_pos_start read from sam file because it is 1-based and python is 0 based.
'''
def read_info_extractor(element): #(sam_file):
    #reads_data = open(sam_file,'r')
    #with( open(sam_file,'r')) as i:
    #    for line in i:
    columns = element.rstrip().split("\t")
    rd_name = columns[0]
    sam_flag = int(columns[1]) #sam flag
    rd_pos_start = int(columns[3])-1 #array(pos).astype(int)
    rd_cigar = columns[5]
    rd_sequence_raw = columns[9]
    #''' SINCE SIMULATIONS ARE ONLY ON CHROMOSOME 17!'''
    #rd_RNAME = columns[2] #needs to be 17
    #rd_RNEXT = columns[6] #needs to be 17 #want properly paired reads without mapping artifacts! #for Translocation will have a different approach!


    #reads_data.close()

    return rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw #, rd_RNAME, rd_RNEXT


def read_content_in_window_fcn(line,coord,window_to_each_side):
    window_arange, left, right = make_window(coord,window_to_each_side)
    rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line) #, rd_RNAME, rd_RNEXT
    read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
    read_cleaned_length = len(read_cleaned)
    read_cleaned_length_computed_end = rd_pos_start + read_cleaned_length

    read_cleaned_length_arange = arange(rd_pos_start,read_cleaned_length_computed_end)
    read_content_in_window = intersect1d(window_arange,read_cleaned_length_arange)

    return read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange


def clean_read_mapper(line,coord,window_to_each_side,rd_sequence_raw,rd_cigar): #,matrix_str #all_reads_in_window_file_name,full_window_length

    #number_of_reads_in_window = number_of_reads_in_window_compute(all_reads_in_window_file_name)
    #matrix_str, matrix_int = make_matrix(number_of_reads_in_window,full_window_length)
    window_arange = make_window(coord,window_to_each_side)[0]
    read_cleaned_length_arange = read_content_in_window_fcn(line,coord,window_to_each_side)[2]

    read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
    read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange = read_content_in_window_fcn(line,coord,window_to_each_side)
    begin_of_intersection = read_content_in_window[0]
    end_of_intersection = read_content_in_window[len(read_content_in_window)-1]

    window_begin_intersect_index = where(window_arange == begin_of_intersection)[0][0]
    window_end_intersect_index = where(window_arange == end_of_intersection)[0][0]
    rd_begin_intersect_index = where(read_cleaned_length_arange == begin_of_intersection)[0][0]
    rd_end_intersect_index = where(read_cleaned_length_arange == end_of_intersection)[0][0]

    #matrix_str[counter][window_begin_intersect_index:window_end_intersect_index+1] = list(read_cleaned)[rd_begin_intersect_index:rd_end_intersect_index+1]

    return window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index #matrix_str


def matrix_read_updater_for_str_int(all_reads_in_window_file_name,coord,window_to_each_side,number_of_reads_in_window_total,full_window_length): #coord is 'i' in first/outermost for-loop
    matrix_str, matrix_int_left_clip, matrix_int_right_clip = make_matrix(number_of_reads_in_window_total,full_window_length) #matrices initialized #, matrix_int
    #print 'shape(matrix_int_left_clip):', shape(matrix_int_left_clip)
    window_arange = make_window(coord,window_to_each_side)[0]
    with ( open(all_reads_in_window_file_name,'r') ) as j:
        rd_counter = 0
        for line in j:
            rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line) #, rd_RNAME, rd_RNEXT
            #''' Need to have header in the sam,bam files that work with, otherwise seems that sambamba filter will yield a '*' for rd_RNAME, rd_RNEXT which is not good '''
            #if rd_RNAME; don't need this anymore since only getting reads that map to chromosome 17 exclusively!
            read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
            read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange = read_content_in_window_fcn(line,coord,window_to_each_side)
            begin_of_intersection = read_content_in_window[0]
            end_of_intersection = read_content_in_window[len(read_content_in_window)-1]
            if len(read_content_in_window) == 0: #''' If read present inside the window! ''' All reads should be present by definition of fetching the reads!
                raise ValueError('primary read outside window')
            window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index = clean_read_mapper(line,coord,window_to_each_side,rd_sequence_raw,rd_cigar)
            matrix_str[rd_counter][window_begin_intersect_index:window_end_intersect_index+1] = list(read_cleaned)[rd_begin_intersect_index:rd_end_intersect_index+1]
            ''' cigar can also be equal to '*' if unavailable; that's why include in regex expression a '*' '''
            letter_order = re.findall(r"[*MSH]",rd_cigar) #have excluded the '*' with sambamba filter! #to find whether Soft or Hard clipped at beginning or end of the read
            if (letter_order[0] == 'S' or letter_order[0] == 'H'): #read already present in window
                #rd_pos_start >= begin_of_intersection and ###WILL NOT BE MAPPING THE MATES IN ANY CASE, JUST GIVING DISTANCES INSTEAD OF SAY ENTRY OF '1'
                #ind_rd_start = where(rd_pos_start == window_arange)[0][0]
                matrix_int_left_clip[rd_counter][window_begin_intersect_index] = 1
                #matrix_int_left_clip[rd_counter][ind_rd_start] = 1
            if (letter_order[len(letter_order)-1] == 'S' or letter_order[len(letter_order)-1] == 'H'): #read already present in window
                #end_of_intersection >= read_cleaned_length_computed_end and
                #ind_rd_end = where(read_cleaned_length_computed_end == window_arange)[0][0]
                matrix_int_right_clip[rd_counter][window_end_intersect_index] = 1
                #matrix_int_right_clip[rd_counter][ind_rd_end] = 1
            ''' make to comment out self.assertEqual(all(matrix_int_left_clip == zeros((79,10))),True) when if rd_counter == 10: break is commented out here below '''
            #if rd_counter == 10:
            #    break


            rd_counter += 1

            #number_of_read_artifacts_found_deviating_from_chr17 = number_of_reads_in_window_total - rd_counter
    #print 'matrix_str:', matrix_str
    #print 'matrix_int_left_clip:', matrix_int_left_clip
    #print 'matrix_int_right_clip:', matrix_int_right_clip
    return matrix_str, matrix_int_left_clip, matrix_int_right_clip #, number_of_read_artifacts_found_deviating_from_chr17 #,matrix_int


def get_exact_matches(position, matrix_str_updated,current_ref): #current_line_on_reference,left,right
    #current_ref = current_reference_fcn(current_line_on_reference,left,right)
    matrix_str_updated_transpose = matrix_str_updated.transpose()
    tally1 = len(where(matrix_str_updated_transpose[position] == current_ref[position])[0]) #right side is a single letter
    #print 'tally1:', tally1
    tally2 = len(where(matrix_str_updated_transpose[position] == '=')[0]) #ACCORDING TO THE SAM FILE, THE SEQ CAN HAVE AN '=' SIGN AT A BASE POSITION WHICH MEANS THAT ITS EQUAL TO THE REFERENCE THERE!
    #print 'tally2:', tally2
    tally = tally1 + tally2

    #print 'tally:', tally
    return tally


def exact_matches_channel_fcn(matrix_str_updated,current_ref):
    #coverage_channel = []
    #for i in range(len(current_reference)):
    #    tally = len(where(matrix_str_updated_transpose[i] == current_reference[i])[0])
    #    coverage_channel.append(tally)

    exact_matches_channel = [get_exact_matches(position,matrix_str_updated,current_ref) for position in range(len(current_ref))]

    return exact_matches_channel


def get_coverage(position, matrix_str_updated):
    matrix_str_updated_transpose = matrix_str_updated.transpose()
    coverage = len(where(matrix_str_updated_transpose[position] != '0')[0])

    return coverage


def coverage_channel_fcn(matrix_str_updated,current_ref):
    coverage_channel = [get_coverage(position, matrix_str_updated) for position in range(len(current_ref))]

    return coverage_channel


def get_clips_right_left(position, matrix_int_updated): #i think should do the updating inside the function
    matrix_int_updated_transpose = matrix_int_updated.transpose()
    clips_right_left = len(where(matrix_int_updated_transpose[position] == 1)[0])

    return clips_right_left


def clipped_rds_left_fcn(matrix_int_left_updated,current_ref):
    clipped_rds_left_channel = [get_clips_right_left(position, matrix_int_left_updated) for position in range(len(current_ref))]

    return clipped_rds_left_channel


def clipped_rds_right_fcn(matrix_int_right_updated,current_ref):
    clipped_rds_right_channel = [get_clips_right_left(position, matrix_int_right_updated) for position in range(len(current_ref))]

    return clipped_rds_right_channel


''' Need to test the functions below this line '''

#Need to code this one after the final 8 channels
#THINK ABOUT HOW TO GET 12*2 + 1 CHANNELS HERE! NEED TO HAVE PAIR CHANNELS HERE AS WELL!
def channels_12_vstacker(matrix_str_updated,matrix_int_left_updated,matrix_int_right_updated,current_ref):
    exact_matches_channel = exact_matches_channel_fcn(matrix_str_updated,current_ref)
    coverage_channel = coverage_channel_fcn(matrix_str_updated,current_ref)
    clipped_rds_left_channel = clipped_rds_left_fcn(matrix_int_left_updated,current_ref)
    clipped_rds_right_channel = clipped_rds_right_fcn(matrix_int_right_updated,current_ref)
    ''' Need the 8 remaining channels here!
    .
    .
    .
    '''

    vstack_12_channels = array(vstack((exact_matches_channel,coverage_channel,clipped_rds_left_channel,clipped_rds_right_channel)),dtype=int)

    return vstack_12_channels
    ''' remember also to include GC content after calling channels_12_vstacker twice on the properly paired pair of files! '''




''' >>> Need these function here below for the last 8 channels: IGNORE '*' FOR THESE CHANNELS! <<< '''
def find_mates(sam_file, id_, counter): #%(GS_bam_30xsubsample_file,rd_name) #counter? or column?
    cmd_mate = "sambamba view -F \"read_name == %s\" %s > reads_with_same_name_%s.sam" %(id_,sam_file,counter)
    call(cmd_mate, shell=True)

def all_unique_read_ids_in_window(bam_file,chr_number,left,right,name,counter): #replace counter with column?
    cmd_id = "sambamba view %s %d:%d-%d | cut -f1 | uniq > idFile_%s_%s.txt" %(bam_file,chr_number,left,right,name,counter) #name = GS,G1,S3,S1...
    call(cmd_id, shell=True)

    with(open('idFile_%s_%s.txt','r') %(name,counter)) as i:
        for line in i:
            id_ = line.rstrip()
            find_mates(sam_file, id_)




''' >>> FILE LOCATIONS <<<'''

''' I think that the clipped reads used for the NO SV category are from the S2, S3N files; have a script where outputed these results '''

''' Replace these files with local machine locations and group into 3 pairs of 2 windows each '''

data_chr17_fasta_file="chr17_human.fasta"
Truth_set_file = "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed"
''' Make sure this bed file excludes the 24 clips that are very close or maybe keep them in? In reality can have overlapping deletions or those that are close by'''

G1_bam_30x_file = "G1_dedup.bam"
#"/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.bam"
GS_bam_30xsubsample_file = "GS_dedup.subsampledto30x.bam"
#"/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.bam"
G1_sam_30x_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.sam"
GS_sam_30xsubsample_file = "GS_dedup.subsampledto30x.sam"
#"/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.sam"

S2_bam_file = "S2_dedup.bam"
#"/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/S2_dedup.bam"
S3N_bam_file  = "S3N_dedup.bam"
#"/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/S3N_dedup.bam"
S2_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S2_dedup_clippedrds.sam"
S3N_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S3N_dedup_clippedrds.sam"

''' MAKE THE NOSV CATEGORY FILE LOCATIONS!!! '''
