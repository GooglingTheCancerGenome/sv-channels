import logging
import re  # regex
import sys
import os
import pysam
import numpy as np
import matplotlib.pyplot as plt
from subprocess import call, Popen

# DEBUG printing
logging.basicConfig(stream=sys.stderr, level=logging.ERROR)

# sambamba = 'path to sambamba'
# cmd0 = 'module load sambamba'
# cmd = f'{sambamba} view {input} {output}'
# run(cmd0 +'&&' + cmd, shell=True)

HPC_MODE = False

'''on local machine'''
# sambamba 0.6.6 must be in PATH
sambambadir = "/Users/lsantuari/Applications/SAM/sambamba/"


# ''' #uncomment for use on hpc
# sambamba v0.5.8 on the HPC
# cmd = "module load sambamcram/sambamba"
# cmd_hpc = "module load sambamcram/sambamba" #on the HPC cluster #sambamcram/samtools"
# call(cmd_hpc, shell=True)
# '''


def bash_command(cmd):
    Popen(['/bin/bash', '-c', cmd])


def get_cigar_coord(cigar):
    re1 = '(\\d+)'  # Integer Number 1
    re2 = '([a-z])'  # Any Single Word Character (Not Whitespace) 1

    rg = re.compile(re1 + re2, re.IGNORECASE | re.DOTALL)
    return [(int(i), v) for i, v in re.findall(rg, cigar)]


def add_dels(seq, cigar):
    '''
    Replacing all instances of deletions
    in Cigar string reported as D by Ds
    '''
    i = 0
    for n, v in get_cigar_coord(cigar):
        if v == 'D':
            seq = seq[:i] + 'D' * n + seq[i:]
            ##print seq
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
    for n, v in get_cigar_coord(cigar):
        if v == 'I':
            seq = seq[:i] + seq[i + n:]
            ##print seq
            i += n
        elif v == 'H':  # or v == 'S'
            pass
        else:
            i += n
    return seq


def rem_softclips(seq_with_rem_ins, cigar):
    '''
    Remove soft clipped sequence, can occur at beginning  and/or end of the read
    '''
    seq = seq_with_rem_ins  # with removed insertions
    i = 0
    for n, v in get_cigar_coord(cigar):
        if v == 'S':
            seq = list(seq)
            del seq[i:i + n]
            ##print seq
            # i += n # i think sequence is now shorter so don't boost by n!
        elif v == 'H' or v == 'I':  # ignore values for insertions #can have a hard clip occuring when have a soft clip in seq too!
            pass
        else:  # if D or M
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
            ##print line
            currentline = ""
        else:
            line = line.rstrip('\n')
            currentline += line
    # data.close()

    ##print 'currentline[0:5]:', currentline[0:5]
    return currentline


def GC_content_dict_fcn(entry):
    GC_content_dict = {'G': 1, 'C': 1, 'A': 0, 'T': 0}
    GC_contents = [GC_content_dict[k] for k in [i for i in entry]]

    return GC_contents


def locations_DEL_INS(truth_file):
    start_SV_DEL = []
    end_SV_DEL = []
    SV_type_DEL = []
    chr_list_DEL = []

    start_SV_INS = []
    # end_SV_INS = []
    SV_type_INS = []
    chr_list_INS = []
    with(open(truth_file, 'r')) as i:
        for line in i:
            line = line.rstrip()
            columns = line.split("\t")
            if columns[4] == "DEL":
                start_SV_DEL.append(int(columns[1]))
                end_SV_DEL.append(int(columns[3]))
                # SV_type_DEL.append("DEL")
                chr_list_DEL.append(int(columns[0]))  # re.findall(r"(\d+)[AB]",columns[0])
                ##print 'chromosome:', chromosome
                # chr_list.append(chromosome)
                # start_end_SV_list.append([int(columns[1]),int(columns[4])])
                # break
            elif str(columns[4]) == "INS":
                start_SV_INS.append(int(columns[1]))
                # end_SV_INS.append(int(columns[4]))
                # SV_type_INS.append("INS")
                chr_list_INS.append(int(columns[0]))

    start_end_SV_DEL = start_SV_DEL + end_SV_DEL
    start_end_SV_DEL_INS = start_end_SV_DEL + start_SV_INS

    return start_SV_DEL, end_SV_DEL, start_SV_INS


def make_window(breakpoint, window_to_each_side):  # breakpoint is the same as coord
    left = breakpoint - window_to_each_side
    right = breakpoint + window_to_each_side
    window_arange = np.arange(left, right)

    return window_arange, left, right


def current_reference_fcn(current_line_on_reference, left, right):
    current_ref = current_line_on_reference[left:right]

    ##print 'current_ref:', current_ref
    return current_ref


def make_matrix(number_of_reads_in_window_total, full_window_length):
    matrix_str = np.array(np.zeros((number_of_reads_in_window_total, full_window_length)), dtype='|S1')
    matrix_int_left_clip = np.zeros((number_of_reads_in_window_total, full_window_length))
    matrix_int_right_clip = np.zeros((number_of_reads_in_window_total, full_window_length))
    # matrix_int = zeros((number_of_reads_in_window_total,full_window_length))

    return matrix_str, matrix_int_left_clip, matrix_int_right_clip  # ,matrix_int


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


def all_reads_in_window(bam_file, chr_number, left, right, name,
                        counter):  # sambamba sorts read by start position automatically
    # ''' FOR HPC USE THE FOLLOWING LINE because on laptop have '*' in the RNAME and RNEXT fields which correspond to sambamba's ref_name and mate_ref_name '''
    ''' cigar != '*' and sequence != '*' doesn't work in sambamba but its ok as have taken these exemption into account and will map '*' in my matrix '''
    ''' and cigar != '*' and sequence != '*' >>> removed from the code as output '*' to the matrix and then will avoid it! So all ok! '''
    '''
    Need to add +1 here too since sam file is 1-based and get situation earlier when fetched reads from 977 - 987 in sambamba command and got a read
    that started at 827 with 150M so 827 + 150 = 977 which is in the range [977,987] but its 1-based so actually its 826 + 150 which is outside the range!
    so need to add 1 to left side and subtract 1 from right side!
    '''
    left1 = left + 1
    right1 = right - 1
    cmd_all_reads_in_window = "sambamba view -F \"ref_name == '17' and mate_ref_name == '17'\" %s %d:%d-%d > all_reads_in_window_%s_newVer_10k_100window.sam" % (
        bam_file, chr_number, left1, right1, name)  # counter removed #name = GS,G1,S3,S1... #_new added 29.11.17

    if (not HPC_MODE):
        cmd_all_reads_in_window = sambambadir + cmd_all_reads_in_window

    # ''' The line below here have been using on local files on my laptop '''
    # cmd_all_reads_in_window = "sambamba view %s %d:%d-%d > all_reads_in_window_%s_%s.sam" %(bam_file,chr_number,left,right,name,counter) #name = GS,G1,S3,S1...
    call(cmd_all_reads_in_window, shell=True)
    all_reads_in_window_file = "all_reads_in_window_%s_newVer_10k_100window.sam" % (
        name)  # counter removed #_new added 29.11.17
    ##print all_reads_in_window_file

    return all_reads_in_window_file


def number_of_reads_in_window_compute(sam_file):
    counter_lines = 0
    sample_file_data = open(sam_file, 'r')
    for line in sample_file_data.xreadlines():
        counter_lines += 1
    sample_file_data.close()

    return counter_lines


'''
Subracting: -1 from rd_pos_start read from sam file because it is 1-based and python is 0 based.
'''


def read_info_extractor(element):  # (sam_file):
    # reads_data = open(sam_file,'r')
    # with( open(sam_file,'r')) as i:
    #    for line in i:
    columns = element.rstrip().split("\t")
    rd_name = columns[0]
    sam_flag = int(columns[1])  # sam flag
    rd_pos_start = int(columns[3]) - 1  # array(pos).astype(int)
    rd_cigar = columns[5]
    rd_sequence_raw = columns[9]
    # ''' SINCE SIMULATIONS ARE ONLY ON CHROMOSOME 17!'''
    # rd_RNAME = columns[2] #needs to be 17
    # rd_RNEXT = columns[6] #needs to be 17 #want properly paired reads without mapping artifacts! #for Translocation will have a different approach!

    # reads_data.close()

    return rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw  # , rd_RNAME, rd_RNEXT


def read_content_in_window_fcn(line, coord, window_to_each_side):
    window_arange, left, right = make_window(coord, window_to_each_side)
    # print 'window_arange:', window_arange
    rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line)  # , rd_RNAME, rd_RNEXT
    read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
    read_cleaned_length = len(read_cleaned)
    # print 'read_cleaned_length:', read_cleaned_length
    read_cleaned_length_computed_end = rd_pos_start + read_cleaned_length

    read_cleaned_length_arange = np.arange(rd_pos_start, read_cleaned_length_computed_end)
    # print 'read_cleaned_length_arange:', read_cleaned_length_arange
    read_content_in_window = np.intersect1d(window_arange, read_cleaned_length_arange)

    return read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange


def clean_read_mapper(line, coord, window_to_each_side, rd_sequence_raw,
                      rd_cigar):  # ,matrix_str #all_reads_in_window_file_name,full_window_length

    # number_of_reads_in_window = number_of_reads_in_window_compute(all_reads_in_window_file_name)
    # matrix_str, matrix_int = make_matrix(number_of_reads_in_window,full_window_length)
    window_arange = make_window(coord, window_to_each_side)[0]
    read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
    read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange = read_content_in_window_fcn(
        line, coord, window_to_each_side)
    begin_of_intersection = read_content_in_window[0]
    end_of_intersection = read_content_in_window[len(read_content_in_window) - 1]

    window_begin_intersect_index = np.where(window_arange == begin_of_intersection)[0][0]
    window_end_intersect_index = np.where(window_arange == end_of_intersection)[0][0]
    rd_begin_intersect_index = np.where(read_cleaned_length_arange == begin_of_intersection)[0][0]
    rd_end_intersect_index = np.where(read_cleaned_length_arange == end_of_intersection)[0][0]

    # matrix_str[counter][window_begin_intersect_index:window_end_intersect_index+1] = list(read_cleaned)[rd_begin_intersect_index:rd_end_intersect_index+1]

    return window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index  # matrix_str


def matrix_read_updater_for_str_int(all_reads_in_window_file_name, coord, window_to_each_side,
                                    number_of_reads_in_window_total,
                                    full_window_length):  # coord is 'i' in first/outermost for-loop

    # print(all_reads_in_window_file_name)
    matrix_str, matrix_int_left_clip, matrix_int_right_clip = make_matrix(number_of_reads_in_window_total,
                                                                          full_window_length)  # matrices initialized #, matrix_int
    ##print 'shape(matrix_int_left_clip):', shape(matrix_int_left_clip)
    window_arange = make_window(coord, window_to_each_side)[0]
    with (open(all_reads_in_window_file_name, 'r')) as j:
        rd_counter = 0
        for line in j:
            # print 'line:', line #line.rstrip().split("\t")[3]
            rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(
                line)  # , rd_RNAME, rd_RNEXT
            # ''' Need to have header in the sam,bam files that work with, otherwise seems that sambamba filter will yield a '*' for rd_RNAME, rd_RNEXT which is not good '''
            # if rd_RNAME; don't need this anymore since only getting reads that map to chromosome 17 exclusively!
            read_cleaned = read_reconstruct(rd_sequence_raw, rd_cigar)
            read_content_in_window, read_cleaned_length_computed_end, read_cleaned_length_arange = read_content_in_window_fcn(
                line, coord, window_to_each_side)

            if len(
                    read_content_in_window) == 0:  # ''' If read present inside the window! ''' All reads should be present by definition of fetching the reads!
                raise ValueError('primary read outside window')
            begin_of_intersection = read_content_in_window[0]
            end_of_intersection = read_content_in_window[len(read_content_in_window) - 1]
            window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index = clean_read_mapper(
                line, coord, window_to_each_side, rd_sequence_raw, rd_cigar)
            matrix_str[rd_counter][window_begin_intersect_index:window_end_intersect_index + 1] = list(read_cleaned)[
                                                                                                  rd_begin_intersect_index:rd_end_intersect_index + 1]
            ''' cigar can also be equal to '*' if unavailable; that's why include in regex expression a '*' '''
            letter_order = re.findall(r"[*MSH]",
                                      rd_cigar)  # have excluded the '*' with sambamba filter! #to find whether Soft or Hard clipped at beginning or end of the read
            if (letter_order[0] == 'S' or letter_order[0] == 'H'):  # read already present in window
                # rd_pos_start >= begin_of_intersection and ###WILL NOT BE MAPPING THE MATES IN ANY CASE, JUST GIVING DISTANCES INSTEAD OF SAY ENTRY OF '1'
                # ind_rd_start = where(rd_pos_start == window_arange)[0][0]
                # print('Left clipped')
                matrix_int_left_clip[rd_counter][window_begin_intersect_index] = 1
                # matrix_int_left_clip[rd_counter][ind_rd_start] = 1
            if (letter_order[len(letter_order) - 1] == 'S' or letter_order[
                len(letter_order) - 1] == 'H'):  # read already present in window
                # end_of_intersection >= read_cleaned_length_computed_end and
                # ind_rd_end = where(read_cleaned_length_computed_end == window_arange)[0][0]
                # print('Right clipped')
                matrix_int_right_clip[rd_counter][window_end_intersect_index] = 1
                # matrix_int_right_clip[rd_counter][ind_rd_end] = 1
            ''' make to comment out self.assertEqual(all(matrix_int_left_clip == zeros((79,10))),True) when if rd_counter == 10: break is commented out here below '''
            # if rd_counter == 10:
            #    break

            rd_counter += 1

            # number_of_read_artifacts_found_deviating_from_chr17 = number_of_reads_in_window_total - rd_counter
    # print 'matrix_str:', matrix_str
    # print(np.shape(matrix_int_left_clip))
    # for i in range(matrix_int_left_clip.shape[1]):
    #    if(sum(matrix_int_left_clip[:][i]) > 0):
    #        print 'matrix_int_left_clip:', matrix_int_left_clip[:][i]
    # print 'matrix_int_right_clip:', matrix_int_right_clip
    return matrix_str, matrix_int_left_clip, matrix_int_right_clip  # , number_of_read_artifacts_found_deviating_from_chr17 #,matrix_int


def get_exact_matches(position, matrix_str_updated, current_ref):  # current_line_on_reference,left,right
    # current_ref = current_reference_fcn(current_line_on_reference,left,right)
    matrix_str_updated_transpose = matrix_str_updated.transpose()
    tally1 = len(
        np.where(matrix_str_updated_transpose[position] == current_ref[position])[0])  # right side is a single letter
    ##print 'tally1:', tally1
    tally2 = len(np.where(matrix_str_updated_transpose[position] == '=')[
                     0])  # ACCORDING TO THE SAM FILE, THE SEQ CAN HAVE AN '=' SIGN AT A BASE POSITION WHICH MEANS THAT ITS EQUAL TO THE REFERENCE THERE!
    ##print 'tally2:', tally2
    tally = tally1 + tally2

    ##print 'tally:', tally
    return tally


def exact_matches_channel_fcn(matrix_str_updated, current_ref):
    # coverage_channel = []
    # for i in range(len(current_reference)):
    #    tally = len(where(matrix_str_updated_transpose[i] == current_reference[i])[0])
    #    coverage_channel.append(tally)

    exact_matches_channel = [get_exact_matches(position, matrix_str_updated, current_ref) for position in
                             range(len(current_ref))]

    return exact_matches_channel


def get_coverage(position, matrix_str_updated):
    matrix_str_updated_transpose = matrix_str_updated.transpose()
    coverage = len(np.where(matrix_str_updated_transpose[position] != '0')[0])

    return coverage


def coverage_channel_fcn(matrix_str_updated, current_ref):
    coverage_channel = [get_coverage(position, matrix_str_updated) for position in range(len(current_ref))]

    return coverage_channel


def get_clips_right_left(position, matrix_int_updated):  # i think should do the updating inside the function
    matrix_int_updated_transpose = matrix_int_updated.transpose()
    clips_right_left = len(np.where(matrix_int_updated_transpose[position] == 1)[0])

    return clips_right_left


def clipped_rds_left_fcn(matrix_int_left_updated, current_ref):
    clipped_rds_left_channel = [get_clips_right_left(position, matrix_int_left_updated) for position in
                                range(len(current_ref))]

    return clipped_rds_left_channel


def clipped_rds_right_fcn(matrix_int_right_updated, current_ref):
    clipped_rds_right_channel = [get_clips_right_left(position, matrix_int_right_updated) for position in
                                 range(len(current_ref))]

    return clipped_rds_right_channel


# Need to code this one after the final 8 channels
# THINK ABOUT HOW TO GET 12*2 + 1 CHANNELS HERE! NEED TO HAVE PAIR CHANNELS HERE AS WELL!
def channels_12_vstacker(matrix_str_updated, matrix_int_left_updated, matrix_int_right_updated, current_ref):
    ''' Need the 8 remaining channels here!
    .
    .
    .
    '''

    exact_matches_channel = exact_matches_channel_fcn(matrix_str_updated, current_ref)
    coverage_channel = coverage_channel_fcn(matrix_str_updated, current_ref)
    clipped_rds_left_channel = clipped_rds_left_fcn(matrix_int_left_updated, current_ref)
    clipped_rds_right_channel = clipped_rds_right_fcn(matrix_int_right_updated, current_ref)

    # print(clipped_rds_left_channel)
    # print(clipped_rds_right_channel)

    vstack_12_channels = np.array(
        np.vstack((exact_matches_channel, coverage_channel, clipped_rds_left_channel, clipped_rds_right_channel)),
        dtype=int)
    # print(vstack_12_channels)

    return vstack_12_channels
    # ''' remember also to include GC content after calling channels_12_vstacker twice on the properly paired pair of files! '''


def vstack_12_channel_pairer_plus_GC_chanel_fcn(vstack_12_channels1, vstack_12_channels2, GC_contents):
    vstack_12_channel_pairer_plus_GC_chanel = np.vstack((vstack_12_channels1, vstack_12_channels2, GC_contents))

    return vstack_12_channel_pairer_plus_GC_chanel


# def proper_category_pairer_fcn():


''' Need to test the functions below this line '''

''' >>> Need these function here below for the last 8 channels: IGNORE '*' FOR THESE CHANNELS! <<< '''

'''
def find_mates(sam_file, id_, counter):  # %(GS_bam_30xsubsample_file,rd_name) #counter? or column?
    cmd_mate = "sambamba view -F \"read_name == %s\" %s > reads_with_same_name_%s.sam" % (
        id_, sam_file, counter)
    bash_command(cmd_mate)


def all_unique_read_ids_in_window(bam_file, chr_number, left, right, name, counter):  # replace counter with column?
    cmd_id = "sambamba view %s %d:%d-%d | cut -f1 | uniq > idFile_%s_%s.txt" % (
        bam_file, chr_number, left, right, name, counter)  # name = GS,G1,S3,S1...
    bash_command(cmd_id)

    with(open('idFile_%s_%s.txt', 'r') % (name, counter)) as i:
        for line in i:
            id_ = line.rstrip()
            find_mates(sam_file, id_)
'''

'''Extract all position of soft/hard clipped reads in a BAM alignment window [chr:start-end]'''


def get_clipped_positions(bamfile, chr, start, end):
    '''

    :param bamfile: filepath of BAM alignment
    :param chr: chromosome
    :param start: start position of window
    :param end: end position of window
    :return: list of unique (Hard or Soft)-clipped positions
    '''

    assert os.path.isfile(bamfile)
    # print('Reading BAM:%s' % bamfile)
    samfile = pysam.AlignmentFile(bamfile, "r")
    # read_count = samfile.count(chr, start, end)
    iter = samfile.fetch(chr, start, end)
    clipped_pos = set()
    for read in iter:
        # print(str(read))
        if not read.is_unmapped:
            if read.cigartuples is not None:
                if read.cigartuples[0][0] in [4, 5]:
                    # print(str(read))
                    # print('Clipped at the start: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                    # print('Pos:%d, clipped_pos:%d' % (read.reference_start, read.get_reference_positions()[0]))
                    clipped_pos.add(read.get_reference_positions()[0] + 1)
                if read.cigartuples[-1][0] in [4, 5]:
                    # print('Clipped at the end: %s -> %s' % (str(read.cigarstring), str(read.cigartuples)))
                    # print('Pos:%d, clipped_pos:%d' %(read.reference_end, read.get_reference_positions()[-1]))
                    clipped_pos.add(read.get_reference_positions()[-1] + 1)
    samfile.close()
    return list(clipped_pos)


def get_clipped_positions_from_CR_BAM(bamfile):
    '''

    :param bamfile: filepath of BAM alignment with Hard/Soft clipped reads only
    :return: dictionary of list of unique genomic positions per chromosome
    '''

    assert os.path.isfile(bamfile)

    # print('Reading BAM:%s' % bamfile)
    samfile = pysam.AlignmentFile(bamfile, "r")
    print(str(samfile.header))

    # chr name from header
    clipped_pos = dict()
    for el in samfile.header['SQ']:
        clipped_pos[el['SN']] = set()
    # read_count = samfile.count(chr, start, end)
    for read in samfile.fetch():
        if (not read.is_unmapped) and (not read.mate_is_unmapped):
            # assert read.cigartuples[0][0] in [4, 5] or read.cigartuples[-1][0] in [4, 5]
            if read.cigartuples[0][0] in [4, 5]:
                cpos = read.get_reference_positions()[0] + 1
                if cpos not in clipped_pos[read.reference_name]:
                    clipped_pos[read.reference_name].add(cpos)
            if read.cigartuples[-1][0] in [4, 5]:
                cpos = read.get_reference_positions()[-1] + 1
                if cpos not in clipped_pos[read.reference_name]:
                    clipped_pos[read.reference_name].add(cpos)
    samfile.close()
    return clipped_pos


'''In a BAM alignment window [chr:start-end], get array with count of reads per position
for reads with absolute insert size larger than 3 standard deviations as estimate by picard'''


def get_del_reads_per_pos(bamfile, chr, start, end):
    '''

    :param bamfile: filepath of BAM alignment
    :param chr: chromosome
    :param start: start position of window
    :param end: end position of window
    :return: array with counts of reads from properly aligned pairs per position
    '''

    median_insert_size, median_standard_deviation = (352, 66)
    print('Reading BAM:%s\n%s:%d-%d' % (bamfile, chr, start, end))
    samfile = pysam.AlignmentFile(bamfile, "r")
    iter = samfile.fetch(chr, start, end)
    cnt_reads = samfile.count(chr, start, end)
    print('Counted: %d reads' % cnt_reads)
    insert_size = []
    assert (end - start) >= 0
    counter = 0
    cnt_array = np.zeros((end - start,), dtype=int)
    print('cnt_array length: %d' % len(cnt_array))
    for read in iter:
        # print(str(read))
        # print('Paired:%s\tUnmapped:%s\tMate_unmapped:%s' % (read.is_paired, read.is_unmapped, read.mate_is_unmapped))
        # print(str(read.reference_start - read.next_reference_start))
        if (not read.is_unmapped) and read.is_paired and (not read.mate_is_unmapped) \
                and read.next_reference_name == read.reference_name:
            if abs(read.next_reference_start - read.reference_start) > (
                    median_insert_size + 3 * median_standard_deviation):
                # print(str(read))
                # print('Mate %s:%d' % (read.next_reference_name, read.next_reference_start))
                cnt_array[read.reference_start - start] += 1
                counter += 1
    print(str(cnt_array))
    print('Processed: %d reads, cnt_array sum:%d' % (counter, sum(cnt_array)))
    return cnt_array


def plot_channels(start_window, n_windows, X_train, y_train):
    '''
    Function to plot channel from Sonya
    :param start_window:
    :param n_windows:
    :param X_train:
    :param y_train:
    :return:
    '''

    number_channels = X_train.shape[1]
    for i in range(start_window, start_window + n_windows):
        print(y_train[i], 'id:', i)
        for j in range(0, number_channels):
            shift = 0
            start = 0
            if j in [0, 1, 4, 5]:
                shift = -60
            if j in [4, 5, 6, 7]:
                start = 70
            Z = [start + shift + x + 5 * j * 4 for x in X_train[i][j]]
            plt.ylim([-65, 250])
            plt.plot(Z)
        # plt.savefig(y_train[i] + '_' + str(i) + '.png')
        plt.show()


def plot_channels_mtx(ch_mtx, title):

    #print(ch_mtx.shape)
    number_channels = ch_mtx.shape[0] - 1
    for j in range(0, number_channels):
        shift = 0
        start = 0
        if j in [0, 1, 4, 5]:
            shift = -60
        if j in [4, 5, 6, 7]:
            start = 70
        Z = [start + shift + x + 5 * j * 4 for x in ch_mtx[j]]
        plt.ylim([-65, 250])
        plt.plot(Z)
    # plt.savefig(y_train[i] + '_' + str(i) + '.png')
    plt.title(title)
    plt.show()


def load_channels():
    '''
    Load saved channel data and generate plots
    :return:
    '''
    wd = '/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/ChannelMaker/DEL_clipped_win/' + \
         'ChannelMaker_DEL_clipped_win/'
    germline_cube = np.load(wd + 'germline_cube_data_file.npy')
    somatic_cube = np.load(wd + 'somatic_cube_data_file.npy')
    germline_label = np.load(wd + 'germline_label_array_file.npy')
    somatic_label = np.load(wd + 'somatic_label_array_file.npy')

    germline_cube_8channels = germline_cube[:, 0:8, :]
    # germline_cube_8channels.shape
    somatic_cube_8channels = somatic_cube[:, 0:8, :]
    # somatic_cube_8channels.shape

    start_window = 15000
    n_windows = 30

    plot_channels(start_window, n_windows, germline_cube_8channels, germline_label)


''' >>> FILE LOCATIONS <<<'''

''' I think that the clipped reads used for the NO SV category are from the S2, S3N files; have a script where outputed these results '''

''' Replace these files with local machine locations and group into 3 pairs of 2 windows each '''

# G1_bam_30x_file = wd + "SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/" + "G1_dedup.bam"
# "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.bam"
# GS_bam_30xsubsample_file = wd + "SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"
# "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.bam"
# G1_sam_30x_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.sam"
# GS_sam_30xsubsample_file = "GS_dedup.subsampledto30x.sam"
# "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.sam"

# S2_bam_file = wd + "SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/" + "S2_dedup.bam"
# "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/S2_dedup.bam"
# S3N_bam_file = wd + "SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/" + "S3N_dedup.bam"
# "reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/S3N_dedup.bam"
# S2_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S2_dedup_clippedrds.sam"
# S3N_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S3N_dedup_clippedrds.sam"


if HPC_MODE:
    wd = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/"

    data_chr17_fasta_file = wd + "chr17_human.fasta"
    Truth_set_file = wd + "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed"

    G1_bam_30x_file = wd + "reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/" + "G1_dedup.bam"
    GS_bam_30xsubsample_file = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"
    S2_bam_file = wd + "reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/" + "S2_dedup.bam"
    S3N_bam_file = wd + "reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/" + "S3N_dedup.bam"

else:
    wd = "/Users/lsantuari/Documents/Data/HPC/DeepSV/Artificial_data/SURVIVOR-master/Debug/"
    data_chr17_fasta_file = wd + "chr17_human.fasta"
    Truth_set_file = wd + "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed"
    ''' Make sure this bed file excludes the 24 clips that are very close or maybe keep them in? In reality can have overlapping deletions or those that are close by'''

    G1_bam_30x_file = wd + "reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/" + "G1_dedup.bam"
    GS_bam_30xsubsample_file = wd + "reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/" + "GS_dedup.subsampledto30x.bam"
    S2_bam_file = wd + "reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/" + "S2_dedup.bam"
    S3N_bam_file = wd + "reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/" + "S3N_dedup.bam"

''' MAKE THE NOSV CATEGORY FILE LOCATIONS!!! '''
