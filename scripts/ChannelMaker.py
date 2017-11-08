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
            print seq
            i += n
        elif v == 'H':
            pass
        else:
            i += n
    #print 'seq after add_dels function:', seq
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
            print seq
            i += n
        elif v == 'H': #or v == 'S'
            pass
        else:
            i += n
    #print 'seq after rem_ins function:', seq
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
            del seq[i:i+n] #del seq[0:3] and then del seq[sum{D,M}-#s : end]
            print seq
            #i += n # i think sequence is now shorter so don't boost by n!
        elif v == 'H' or v == 'I': #ignore values for insertions #can have a hard clip occuring when have a soft clip in seq too!
            pass
        else: #if D or M
            i += n
    #print 'seq after rem_softclips function:', ''.join(seq)
    return ''.join(seq)


def Ale_reconstruct(seq, cigar):
    '''
    Remove insertions after add deletions
    '''
    #print 'answer:', rem_ins(add_dels(seq, cigar), cigar)
    #return rem_ins(add_dels(seq, cigar), cigar) #rem_softclips(rem_ins(add_dels(seq, cigar), cigar), cigar)
    #print 'answer:', rem_softclips(rem_ins(add_dels(seq, cigar), cigar), cigar) --> doesn't work, last function not initialized!
    #answer = rem_softclips(rem_ins(add_dels(seq, cigar), cigar), cigar) --> this works!
    #print 'answer:', answer
    return rem_softclips(rem_ins(add_dels(seq, cigar), cigar), cigar)

'''
#>>> Ale_reconstruct Testing Area <<<
read4 = 'ATCCCGCTACGGG'
read3 = 'ATCCCGCTACGGGAAAA'
read5 = 'GGGATCCCGCTACGGG'
print 'read5:', read5
read2 = 'GGGATCCCGCTACGGGAAAA' #read string
print 'read2:', read2
read1 = 'ATCCCGCTACGGG'
read6 = 'ATCCCGCTACGGG'

cigar4 = '3H2M1D3M2I2M5D2M4D2M5H'
cigar3 = '3H2M1D3M2I2M5D2M4D2M4S'
cigar5 = '3S2M1D3M2I2M5D2M4D2M4H'
print 'cigar5:', cigar5
cigar2 = '3S2M1D3M2I2M5D2M4D2M4S' #cigar
print 'cigar2:', cigar2
cigar1 = '2M1D3M2I2M5D2M4D2M'
cigar6 = '2M1D3M2I3D2M5D2M4D2M' #haven't seen this kind of cigar string with deletion following insertion; code treats this correctly! checked!

r = 'ATDCCCDDDTADDDDDCGDDDDGG'
#'ATDCCCTADDDDDCGDDDDGG' #expected result
print 'expected result:', r

print get_cigar_coord(cigar6)

r1 = Ale_reconstruct(read6, cigar6)
print r == r1

#sys.exit(0)
'''


def one_liner(data):
    '''
    Currentline is now chr17; in the .fasta file the sequence
    is split over lines, this combines all lines into one line
    '''
    currentline = ""
    for line in data.xreadlines():
        #print line[0:10]
        if line.startswith('>'):
            line = line.rstrip('\n')
            print line
            currentline = ""
        else:
            line = line.rstrip('\n')
            currentline += line
    data.close()

    return currentline

''' >>> Some parameters <<<'''
window = 100 #500
print 'window:', window
full_window = window*2 #+1
print 'full_window:', full_window
max_coverage_depth_from_30x = 160 #hardcoded number based on 6*30; expect at most 100 reads at one location; matrix below is just for testing really.
print 'max_coverage_depth_from_30x:', max_coverage_depth_from_30x
illumina_rd_length_fixed = 150
print 'illumina_rd_length_fixed:', illumina_rd_length_fixed
chromosome_number = 17
print chromosome_number
''' >>> Some parameters <<<'''

data=open("chr17_human.fasta",'r')
current_line = one_liner(data)
print 'len(current_line):', len(current_line)
print current_line[0:10] #current_line[4123760-1:4123759+150] #current_line[0:10] #should agree with reference sequence at chr17


''' INS aren't used to train the network but their occurence is important as leads to clipped reads; so for NO SV category need to take INS into account too '''

Truth_set_file = "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed" #contains DELs and INs 4.815 DELs and 5.185 INS, INS we don't consider for the moment
#lines in this file look like:  17	982	17	1910	INS

start_SV = []
end_SV = []
SV_type = []
chr_list = []

start_SV_INS = []
#end_SV_INS = []
SV_type_INS = []
chr_list_INS = []

with(open(Truth_set_file,'r')) as i:
    for line in i:
        line = line.rstrip()
        columns = line.split("\t")
        if str(columns[3]) == "DEL" :
            start_SV.append(int(columns[1]))
            end_SV.append(int(columns[4]))
            SV_type.append(str(columns[3]))
            chr_list.append(int(columns[0])) #re.findall(r"(\d+)[AB]",columns[0])
            #print 'chromosome:', chromosome
            #chr_list.append(chromosome)
            #start_end_SV_list.append([int(columns[1]),int(columns[4])])
            #break
        elif str(columns[3]) == "INS" :
            start_SV_INS.append(int(columns[1]))
            #end_SV_INS.append(int(columns[4]))
            SV_type_INS.append(str(columns[3]))
            chr_list_INS.append(int(columns[0]))

start_end_SV = start_SV + end_SV #''' list of all start positions and all end positions of all simulated deletions '''
print 'len(start_SV):', len(start_SV)
print 'len(end_SV):', len(end_SV)
print 'len(SV_type) == 4815:', len(SV_type) == 4815 #number of simulated DELs
print 'len(start_end_SV):', len(start_end_SV)
print 'len(start_end_SV) == 9630:', len(start_end_SV) == 9630

start_end_SV_DEL_INS = start_end_SV + start_SV_INS #both DELs and INS
print 'len(start_end_SV_DEL_INS):', len(start_end_SV_DEL_INS)

#start_SV = array(start_SV)
#end_SV = array(end_SV)

#''' Need to have complimentary .sam files to the .bam files! '''
#''' Apparently first need to convert bam to sam to extract forward and reverse strand reads '''

''' I think that the clipped reads used for the NO SV category are from the S2, S3N files; have a script where outputed these results '''


G1_bam_30x_file="/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.bam"
GS_bam_30xsubsample_file="/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.bam"
G1_sam_30x_file="/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline1_mapped/G1/mapping/G1_dedup.sam"
GS_sam_30xsubsample_file="/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Germline2_Somatic1_mapped/GS/mapping/GS_dedup.subsampledto30x.sam"

S2_bam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic2_mapped/S2/mapping/S2_dedup.bam"
S3N_bam_file  = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/reads_chr17_SURV10kDEL_INS_Somatic3_new_mapped/S3N/mapping/S3N_dedup.bam"
S2_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S2_dedup_clippedrds.sam"
S3N_sam_file = "/hpc/cog_bioinf/ridder/users/cshneider/SURVIVOR-master/Debug/chr17_SURV10kDEL_INS_truth_clipped_vcfs_4callers/S3N_dedup_clippedrds.sam"

#''' Need to check that get matrix of height coverage depth and length of 200 for the window size, don't need to sort the sam file'''

''' Since extract paired reads below, need to check for positions so that don't include paired rds in window'''

###
#''' embedd in 2d matrix '''
###

cmd="module load sambamcram/sambamba"
#"module load sambamcram/samtools" #"echo "+string

GC_content_dict = {'G':1,'C':1,'A':0,'T':0} #''' dictionary for GC content '''

rd_pos_starts_GS_list = []
rd_pos_starts_G1_list = []
counter = 0
''' >>> Outer for-loop is the windowing part: fetch all reads that are mapped to
region plus their cousins, recalibrate all read sequences, fit properly into the window, handle cousins <<< '''
for i in start_end_SV[0]: #''' so just testing on one element; need to comment this out! '''
    left = i - window #coordinate 0 of window #if i was 100; [0:201]
    right = i + window #coordinate 200 of window (length of window is 201, -100, 0, +100) where '0' replaces i
    window_arange = arange(left,right)
    current_reference = str(current_line[left:right]) #just to be sure that not modifying the original
    GC_content = [GC_content_dict[k] for k in [i for i in current_reference]] #[GC_content_dict[k] for k in [i for i in 'GCTT']]
    print 'len(GC_content) == 200:', len(GC_content) == 200

    histogram_channel_GS = []
    histogram_channel_G1 = []

    string1="sambamba view %s %d:%d-%d | cut -f1 >> idFile_GS_%s.txt" %(GS_bam_30xsubsample_file,chr_list[counter],left,right,counter) #should be 17
    string2="LC_ALL=C grep -w -F -f idFile_GS_%s.txt < %s >> subset_GS_%s.sam" %(counter,GS_sam_30xsubsample_file,counter) #can convert sam to bam and view with IGV
    #string2="grep -f idFile_GS.txt %s >> output_GS.sam" %(GS_sam_30xsubsample_file)
    string2p="cat subset_GS_%s.sam | grep -v '#'| sort -n -k4,1 > sorted_subset_GS_%s.sam" %(counter) #sorting for easier comparison with IGV visualization
    #file should be headerless anyway since sambamba view is w/o the -h flag i guess as for samtools
    ''' want to sort on start positions so that can compare easier with IGV plot to see if everything is correct '''
    string3="sambamba view %s %d:%d-%d | cut -f1 >> idFile_G1_%s.txt" %(G1_bam_30x_file,chr_list[counter],left,right,counter) #%d is for ints #"python is %s" %("cool")
    string4="LC_ALL=C grep -w -F -f idFile_G1_%s.txt < %s >> subset_G1_%s.sam" %(counter,G1_sam_30x_file,counter)
    string4p="cat subset_G1_%s.sam | grep -v '#'| sort -n -k4,1 > sorted_subset_G1_%s.sam" %(counter)
    #string4="grep -f idFile_G1.txt %s >> output_G1.sam" %(G1_sam_30x_file)
    os.system(cmd+"\n"+string1+"\n"+string2+"\n"+string2p+"\n"+string3+"\n"+string4+"\n"+string4p) #one shell! #a new line character will tell program to wait for next step!

    ###with( open("idFile_GS_%s.txt" %(counter),'r')) as i: #this is the number of reads w/o their pairs; will need to cut away unnessary zeros later
        ###height_GS = len(i.readlines())
        ###print 'height_GS:', height_GS
    #matrix_GS = array(zeros((height_GS,full_window)),dtype='|S1')
    matrix_GS = array(zeros((max_coverage_depth_from_30x,full_window)),dtype='|S1') #''' for purpose of histogram '''
    matrix_G1 = array(zeros((max_coverage_depth_from_30x,full_window)),dtype='|S1') #''' for purpose of histogram '''
    #array entries hold strings of length 1 #this will allow direct string substitution via: M = array(ones((2,3)),dtype='S1'), M = array(ones((2,3)),dtype='S1'), M[1][0:2] = 'HHH'
    #need to see how clipped reads are being displayed!

    #''' Need to also work with column 10 of sam file: actual string! '''
    #''' Need to check if have S or H because for H these bases are not in the actual string whereas for S they are! '''

    #''' Keep to consistent strategy for all reads'''

    #rd_sequence_len_list = [] #for extra checking purposes
    with( open("sorted_subset_GS_%s.sam" %(counter),'r')) as i: #''' can make this into a function to read each file at a time '''
        for line in i:
            #if line.startswith('#'): #taken care of in string2p
            #    pass
            #else:
            line = line.rstrip()
            columns = line.split("\t")
            rd_pos_start = int(columns[3]) #array(pos).astype(int)
            print 'rd_pos_start:', rd_pos_start
            rd_cigar = str(columns[5])
            print 'rd_cigar:', rd_cigar
            rd_sequence_raw = str(columns[9])
            rd_sequence_ref_mapped = Ale_reconstruct(rd_sequence_raw, rd_cigar) #''' simulates how mapped read would look like in IGV '''
            len_rd_sequence_ref_mapped = len(rd_sequence_ref_mapped) #''' now know the actual length of the mapped read'''
            rd_sequence_ref_mapped_arange = arange(rd_pos_start,rd_pos_start+len_rd_sequence_ref_mapped)

            ''' need to sort/filter reads that are pairs in the sam file which are actually outside the window,
            idea is to intersect mapped read coordinates with window coordinates, i think for split read, 1st part end with a clip,
            the 2nd half begins with a clip and the paired rd corrsponding to this split read is also mapped to a location'''

            read_content_in_window = intersect1d(window_arange,rd_sequence_ref_mapped_arange)
            if len(read_content_in_window) != 0: #''' If read present in the window '''
                begin_of_intersection = read_content_in_window[0]
                end_of_intersection = read_content_in_window[len(read_content_in_window)-1]
                window_begin_intersect_index = where(window_arange == begin_of_intersection)[0][0]
                window_end_intersect_index = where(window_arange == end_of_intersection)[0][0]
                rd_begin_intersect_index = where(rd_sequence_ref_mapped_arange == begin_of_intersection)[0][0]
                rd_end_intersect_index = where(rd_sequence_ref_mapped_arange == end_of_intersection)[0][0]

                matrix_GS[counter][window_begin_intersect_index:window_end_intersect_index+1] = list(rd_sequence_ref_mapped)[rd_begin_intersect_index:rd_end_intersect_index+1]

                counter +=1 #''' this counter is for the rows of the matrix which contains the mapped reads in the given window at hand '''

                ''' look for cousin '''

            #''' Also need to vstack all reads constituting a single window! No, because embedded in matrix! '''
            #''' Need correct matrix embedding so that can do proper histogram by aligning to the reference! '''

            #''' Need to account for which part of the read makes it to be considered inside given window!
            #rd can start earlier then window So need to cut it correctly! '''


            else: #these go outside the histogram at the given window!
                ''' for other channels need this statement #then have paired end read outside the window '''

                ''' look for cousin '''
                pass #''' change this when introduce the other channels '''

    matrix_GS_transpose = matrix_GS.transpose()
    matrix_G1_transpose = matrix_G1.transpose()
    for i in range(len(current_reference)):
        tally_GS = len(where(matrix_GS_transpose[i] == current_reference[i])[0])
        tally_G1 = len(where(matrix_G1_transpose[i] == current_reference[i])[0])
        histogram_channel_GS.append(tally_GS)
        histogram_channel_G1.append(tally_G1)




    '''final_25channel_matrix = vstack((histogram_channel_GS,...,histogram_channel_G1,GC_content)) #GC content will be 25th channel '''
'''
First work with GS and G1 and then make functions to deal generically with all the files
Same kind of approach will work with somatic and germline category but need the clipped read approach that have used previously (including INS)
to obtain areas where can form the NO SV category.
'''

'''
ACTUALLY NEED TO HAVE FUNCTIONS THAT CREATE THESE INDIVIDUAL CHANNELS, ESPECIALLY THAT HAVE 6 OR SO FILES THAT ARE USING

Refactor after tested that works!

Actually, don't have indels at the moment; have reference + SVs. Can introduce indels later; but then need to recalibrate the truth coordinates so that know where
simulated SVs occur!


uber_list_channels
uber_list_labels
*saved as .npy*
'''
#matrix_GS[counter_rdsin][internal_start:internal_end+1] = read_string_final

#counter_rdsin += 1
''' Good idea to keep track of reads that made it in perhaps so then can compare to total reads that had! '''
#''' Need to make sure that this is in the right place! '''

#matrix_G1 = array(zeros((max_coverage_depth_from_30x,full_window)),dtype='|S1')



#''' ***** hmm, expressed like so, left can occur before the actual start of the string! so need to convert to window coordinates ***** '''
#read_string_final_final = [left:right+1] #so length is 201
