import unittest
from ChannelMakerNew_lib import *

''' >>> Testing read cleaning / recalibration <<< '''

read7 = 'ATCCCGCTACGGG'
cigar7 = '2M1D3M3D2I2M5D2M4D2M' #unlikely string combination with deletion following insertion but good to test it

read6 = 'ATCCCGCTACGGG'
cigar6 = '2M1D3M2I3D2M5D2M4D2M' #unlikely string combination with deletion following insertion but good to test it

read5 = 'GGGATCCCGCTACGGG'
cigar5 = '3S2M1D3M2I2M5D2M4D2M4H'

read4 = 'ATCCCGCTACGGG'
cigar4 = '3H2M1D3M2I2M5D2M4D2M5H'

read3 = 'ATCCCGCTACGGGAAAA'
cigar3 = '3H2M1D3M2I2M5D2M4D2M4S'

read2 = 'GGGATCCCGCTACGGGAAAA'
cigar2 = '3S2M1D3M2I2M5D2M4D2M4S'

read1 = 'ATCCCGCTACGGG'
cigar1 = '2M1D3M2I2M5D2M4D2M'

expected_read_output2 = 'ATDCCCDDDTADDDDDCGDDDDGG'
expected_read_output1 = 'ATDCCCTADDDDDCGDDDDGG'

''' >>> End of testing read cleaning / recalibration <<< '''

class TestReadRecal(unittest.TestCase):

    def test_cigar(self):
        self.assertEqual(read_reconstruct(read7, cigar7), expected_read_output2)
        self.assertEqual(read_reconstruct(read6, cigar6), expected_read_output2)
        self.assertEqual(read_reconstruct(read5, cigar5), expected_read_output1)
        self.assertEqual(read_reconstruct(read4, cigar4), expected_read_output1)
        self.assertEqual(read_reconstruct(read3, cigar3), expected_read_output1)
        self.assertEqual(read_reconstruct(read2, cigar2), expected_read_output1)
        self.assertEqual(read_reconstruct(read1, cigar1), expected_read_output1)

        self.assertEqual(read_reconstruct('*', cigar1), '*')
        self.assertEqual(read_reconstruct(read1, '*'), '*')
        self.assertEqual(read_reconstruct('*', '*'), '*')

"""
data_chr17_fasta=open("chr17_human.fasta",'r') #need proper path for .fasta file
current_line_on_reference = one_liner(data_chr17_fasta)
data_chr17_fasta.close()
len_data_chr17_fasta = len(current_line_on_reference)
#print len_data_chr17_fasta
"""

class TestOneLinerRefChr17(unittest.TestCase):

    def test_one_liner(self):
        data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
        current_line_on_reference = one_liner(data_chr17_fasta_file_data)
        data_chr17_fasta_file_data.close()
        len_data_chr17_fasta = len(current_line_on_reference)

        self.assertEqual(current_line_on_reference[0:10],'AAGCTTCTCA')
        self.assertEqual(current_line_on_reference[len_data_chr17_fasta-10:len_data_chr17_fasta],'TGGGTGTGGT')

    def test_current_reference_fcn(self):
        data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
        current_line_on_reference = one_liner(data_chr17_fasta_file_data)
        data_chr17_fasta_file_data.close()
        len_data_chr17_fasta = len(current_line_on_reference)

        self.assertEqual(current_reference_fcn(current_line_on_reference,0,10),'AAGCTTCTCA')
        self.assertEqual(current_reference_fcn(current_line_on_reference,len_data_chr17_fasta-10,len_data_chr17_fasta),'TGGGTGTGGT')

        self.assertEqual(current_reference_fcn(current_line_on_reference,977,987),'GCTTGAGCCC')


class TestGCcontentDict(unittest.TestCase):

    def test_GC_content_dict(self):
        GC_content_test1 = 'GCCATTCCG'
        expected_GC_output = [1,1,1,0,0,0,1,1,1]
        self.assertEqual(GC_content_dict(GC_content_test1),expected_GC_output)


    def test_GC_content(self):
        GC_content_test1 = 'GCCATTCCG'
        expected_GC_output = [1,1,1,0,0,0,1,1,1]
        self.assertEqual(len(GC_content_dict(GC_content_test1)),9)


class TestBreakPointPresence(unittest.TestCase):

    def test_break_point_presence(self):
        Truth_set_file = "chr17_somaticallele_10k_INS_DEL.copy.sorted.bed"
        self.assertEqual(len(locations_DEL_INS(Truth_set_file)[0]),9630) #4815; number of simulated DELs * 2 = 9630
        self.assertEqual(len(locations_DEL_INS(Truth_set_file)[1]),14815) #9630 + 5185 INS bps = 14815 total bps for DEL and INS
        #print locations_DEL_INS(Truth_set_file)[0][0:10]


class TestWindow(unittest.TestCase):

        def test_window(self):
            window_to_each_side_example = 100
            self.assertEqual(len(make_window(0,window_to_each_side_example)[0]),200)

class TestMakeMatrix(unittest.TestCase):

        def test_make_matrix(self):
            self.assertEqual(shape(make_matrix(2,3)[0]),(2,3))
            self.assertEqual(shape(make_matrix(2,3)[1]),(2,3))
            self.assertEqual(make_matrix(2,3)[0][0][0],'0')
            self.assertEqual(make_matrix(2,3)[1][0][0],0)
            self.assertEqual(make_matrix(2,3)[2][0][0],0)
            #self.assertEqual(make_matrix(2,3)[3][0][0],0)

''' NOTE: IF SAMBAMBA VIEW -F FILTER EVER WORKS TO EXCLUDE CIGAR AND SEQUENCE WITH '*' THEN WILL HAVE 78 AND NOT 79 IN THE TESTS BELOW! '''
class TestAllReadsinWindow(unittest.TestCase):
        def test_all_reads_in_window(self):
            counter_ex = 0
            all_reads_in_window_example = open(all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0),'r')
            for line in all_reads_in_window_example.xreadlines():
                counter_ex += 1
                if counter_ex == 79:
                    last_line = line.rstrip().split("\t")
            all_reads_in_window_example.close()
            self.assertEqual(counter_ex,79)
            self.assertEqual(last_line[3],'984')

        def test_number_of_reads_in_window(self):
            number_of_reads = number_of_reads_in_window_compute(all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0))
            self.assertEqual(number_of_reads,79)

        def test_read_info_extractor(self): #just for the first line
            all_reads_in_window_example_data = open(all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0),'r')
            line1 = all_reads_in_window_example_data.readline()
            all_reads_in_window_example_data.close()
            rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line1) #, rd_RNAME, rd_RNEXT
            self.assertEqual(rd_name,'chr17B-3249738')
            self.assertEqual(sam_flag,99)
            self.assertEqual(rd_pos_start,828)
            self.assertEqual(rd_cigar,'150M')
            self.assertEqual(len(rd_sequence_raw),150)

        def test_read_content_in_window(self):
            all_reads_in_window_example_data = open(all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0),'r')
            line1 = all_reads_in_window_example_data.readline() #this will keep reading line by line
            all_reads_in_window_example_data.close()
            line7 = 'chr17A-9860714\t73\t17\t843\t60\t140M10S\t=\t843\t0\tCTGGTGGTGGAAACAAGACTGTCCCAGCCTGGGTGATACAGCGAGACCCCATCTCTACCAAAAAATTAAAAATTAGCTGGGCATGGTGGTGCATGCCTGTAGTCCCAGCTATTCACAGTGCTGAGGTGGGAAGATGCTTGCTCATCCCTC\tFAFFFKKKKKKKKKKKKKKKKKK,KKKKKKKKKKKKKKKKFKKKKKKKFKKKKKKKFFKKKKKKKKKKKKKKKKKKKKKKKKKAKKKKKKK7KKKKKKKK<K,KKKKKKKKKKKKKKKAKF7KKKKK<AFKFKKKKKK<7KKFKKK,,,<\tNM:i:0\tMD:Z:140\tAS:i:140\tXS:i:52\tRG:Z:GS_SURVchr17_S01_L001_001\n'
            read_content_in_window_output1, read_cleaned_length_computed_end1, read_cleaned_length_arange1 = read_content_in_window_fcn(line1,982,5)
            read_content_in_window_output2, read_cleaned_length_computed_end2, read_cleaned_length_arange2 = read_content_in_window_fcn(line7,982,5)
            self.assertEqual(all(read_content_in_window_output1 == array([977])),True) #828 + 150 = 978, so 978 - 1 = 977; from the above test function case used.
            self.assertEqual(read_cleaned_length_computed_end1,978)
            self.assertEqual(all(read_content_in_window_output2 == array([977,978,979,980,981])),True) #843-1+140-1=981
            self.assertEqual(read_cleaned_length_computed_end2,982)
            self.assertEqual(all(read_cleaned_length_arange1 == arange(828,978)),True)
            self.assertEqual(all(read_cleaned_length_arange2 == arange(842,982)),True)


class TestCleanReadMapper(unittest.TestCase):
        def test_clean_read_mapper(self):
            line7 = 'chr17A-9860714\t73\t17\t843\t60\t140M10S\t=\t843\t0\tCTGGTGGTGGAAACAAGACTGTCCCAGCCTGGGTGATACAGCGAGACCCCATCTCTACCAAAAAATTAAAAATTAGCTGGGCATGGTGGTGCATGCCTGTAGTCCCAGCTATTCACAGTGCTGAGGTGGGAAGATGCTTGCTCATCCCTC\tFAFFFKKKKKKKKKKKKKKKKKK,KKKKKKKKKKKKKKKKFKKKKKKKFKKKKKKKFFKKKKKKKKKKKKKKKKKKKKKKKKKAKKKKKKK7KKKKKKKK<K,KKKKKKKKKKKKKKKAKF7KKKKK<AFKFKKKKKK<7KKFKKK,,,<\tNM:i:0\tMD:Z:140\tAS:i:140\tXS:i:52\tRG:Z:GS_SURVchr17_S01_L001_001\n'
            rd_name, sam_flag, rd_pos_start, rd_cigar, rd_sequence_raw = read_info_extractor(line7) #, rd_RNAME, rd_RNEXT
            window_begin_intersect_index, window_end_intersect_index, rd_begin_intersect_index, rd_end_intersect_index = clean_read_mapper(line7,982,5,rd_sequence_raw,rd_cigar)
            self.assertEqual(window_begin_intersect_index,0)
            self.assertEqual(window_end_intersect_index,4)
            self.assertEqual(rd_begin_intersect_index,135)
            self.assertEqual(rd_end_intersect_index,139)


class TestMatrixReadUpdater(unittest.TestCase):
        def test_matrix_read_updater_for_str_int(self):
            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0) #all reads here happen to be correctly on chr17 so get 79 rds!
            matrix_str, matrix_int_left_clip, matrix_int_right_clip = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10) #, number_of_read_artifacts_found_deviating_from_chr17
            self.assertEqual(all(matrix_str[0] == array(['G','0','0','0','0','0','0','0','0','0'])),True)
            self.assertEqual(all(matrix_str[6] == array(['G','C','T','T','G','0','0','0','0','0'])),True)
            self.assertEqual(shape(matrix_int_left_clip), (79,10))
            self.assertEqual(shape(matrix_int_right_clip), (79,10))
            #self.assertEqual(all(matrix_int_left_clip == zeros((79,10))),True) #true when the 'rd_counter = 10 then break' condition uncommented
            self.assertEqual(all(matrix_int_right_clip[6] == array([0,0,0,0,1,0,0,0,0,0])),True)
            self.assertEqual(all(matrix_int_left_clip[75] == array([0,0,0,0,1,0,0,0,0,0])),True)
            #print 'matrix_int_left_clip:', matrix_int_left_clip



            #'TG' should be for line1
            #'TGCTTG' should be for line7
            #both lines starting at the very beginning of the window
            #len('TCTCTGTGTTGATTCTGGTGGTGGAAACAAGACTGTCCCAGCCTGGGTGATACAGCGAGACCCCATCTCTACCAAAAAATTAAAAATTAGCTGGGCATGGTGGTGCATGCCTGTAGTCCCAGCTATTCACAGTGCTGAGGTGGGAAGATG')
            #then test the 2nd part of this function for the location of the clipped reads, readline() reveals that will only have a single 1 in the matrix of right clipped only, all zeros for left matrix!

class TestChannels(unittest.TestCase):

        def test_get_exact_matches(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_str_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[0]
            ###print 'matrix_str_updated:', matrix_str_updated
            #print 'matrix_str_updated.transpose()[0]:', matrix_str_updated.transpose()[0]
            current_ref = current_reference_fcn(current_line_on_reference,977,987) #current_line_on_reference = one_liner(data_chr17_fasta) as defined atop of this file
            print 'current_ref:', current_ref
            current_ref_trial = current_reference_fcn(current_line_on_reference,976,986)
            #print 'current_ref_trial:', current_ref_trial #to check that indeed mapping is off by one and that is because start reported in sam file is 1-based and not 0-based!
            tally = get_exact_matches(0, matrix_str_updated,current_ref)
            #print 'tally:', tally

            self.assertEqual(tally,65)

        def test_exact_matches_channel_fcn(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_str_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[0]
            current_ref = current_reference_fcn(current_line_on_reference,977,987)
            exact_matches_channel = exact_matches_channel_fcn(matrix_str_updated,current_ref)
            #print exact_matches_channel
            self.assertEqual(exact_matches_channel[0],65)
            self.assertEqual(len(exact_matches_channel),10)
            self.assertEqual(all(array(exact_matches_channel) == array([65, 64, 64, 64, 71, 61, 63, 62, 61, 61])),True)


        def test_get_coverage(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_str_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[0]

            coverage = get_coverage(4, matrix_str_updated)
            self.assertEqual(coverage,74)


        def test_coverage_channel_fcn(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_str_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[0]
            current_ref = current_reference_fcn(current_line_on_reference,977,987)

            coverage_channel = coverage_channel_fcn(matrix_str_updated,current_ref)
            #print 'coverage_channel:', coverage_channel
            self.assertEqual(coverage_channel[4],74)
            self.assertEqual(len(coverage_channel),10)
            self.assertEqual(all(array(coverage_channel) == array([65, 64, 64, 64, 74, 62, 63, 62, 61, 61])),True)

            """
            first coordinate that based all the tests here on actually corresponds to an insertion! but that's fine for these tests.
            for testing a DEL event because 982 with corresponding 977 and 987 is an insertion event; also has clipped rds; INS does not show drop in coverage as expect.
            all_reads_in_window_file_name1 = all_reads_in_window(GS_bam_30xsubsample_file,17,8426,8436,'GS','DEL_0')
            matrix_str_updated1 = matrix_read_updater_for_str_int(all_reads_in_window_file_name,8431,5,64,10)[0]
            print 'matrix_str_updated1:', matrix_str_updated1
            current_ref1 = current_reference_fcn(current_line_on_reference,8426,8436)

            coverage_channel1 = coverage_channel_fcn(matrix_str_updated,current_ref)
            print 'coverage_channel1:', coverage_channel1
            #self.assertEqual(coverage_channel1[4],74)
            self.assertEqual(len(coverage_channel1),10)
            #self.assertEqual(all(array(coverage_channel1) == array([65, 64, 64, 64, 74, 62, 63, 62, 61, 61])),True)
            """

        def test_get_clips_right_left(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_int_left_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[1]
            ###print 'matrix_int_left_updated:', matrix_int_left_updated
            matrix_int_right_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[2]
            ###print 'matrix_int_right_updated:', matrix_int_right_updated

            clips_right_left_L = get_clips_right_left(4, matrix_int_left_updated)
            print 'clips_right_left_L:', clips_right_left_L
            clips_right_left_R = get_clips_right_left(4, matrix_int_right_updated)
            print 'clips_right_left_R:', clips_right_left_R
            self.assertEqual(clips_right_left_L,8)
            self.assertEqual(clips_right_left_R,10)

        def test_clipped_rds_left_fcn(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_int_left_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[1]
            current_ref = current_reference_fcn(current_line_on_reference,977,987)

            clipped_rds_left_channel = clipped_rds_left_fcn(matrix_int_left_updated,current_ref)

            self.assertEqual(clipped_rds_left_channel[4],8)
            self.assertEqual(len(clipped_rds_left_channel),10)
            self.assertEqual(all(array(clipped_rds_left_channel) == array([0, 0, 0, 0, 8, 0, 0, 0, 0, 0])),True)


        def test_clipped_rds_right_fcn(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()
            #len_data_chr17_fasta = len(current_line_on_reference)

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_int_right_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)[2]
            current_ref = current_reference_fcn(current_line_on_reference,977,987)

            clipped_rds_right_channel = clipped_rds_right_fcn(matrix_int_right_updated,current_ref)

            self.assertEqual(clipped_rds_right_channel[4],10)
            self.assertEqual(len(clipped_rds_right_channel),10)
            self.assertEqual(all(array(clipped_rds_right_channel) == array([0, 0, 0, 0, 10, 0, 0, 0, 0, 0])),True)

class TestVstackChannels(unittest.TestCase):
        def test_channels_12_vstacker(self):
            data_chr17_fasta_file_data = open(data_chr17_fasta_file,'r')
            current_line_on_reference = one_liner(data_chr17_fasta_file_data)
            data_chr17_fasta_file_data.close()

            all_reads_in_window_file_name = all_reads_in_window(GS_bam_30xsubsample_file,17,977,987,'GS',0)
            matrix_str_updated, matrix_int_left_updated, matrix_int_right_updated = matrix_read_updater_for_str_int(all_reads_in_window_file_name,982,5,79,10)
            current_ref = current_reference_fcn(current_line_on_reference,977,987)

            vstack_12_channels = channels_12_vstacker(matrix_str_updated,matrix_int_left_updated,matrix_int_right_updated,current_ref)
            print vstack_12_channels
            '''
            ORDER FROM TOP TO BOTTOM IS: exact_matches_channel,coverage_channel,lipped_rds_left_channel,clipped_rds_right_channel
            '''
            self.assertEqual(all(vstack_12_channels[0] == array(exact_matches_channel_fcn(matrix_str_updated,current_ref),dtype=int)),True)
            self.assertEqual(all(vstack_12_channels[1] == array(coverage_channel_fcn(matrix_str_updated,current_ref),dtype=int)),True)
            self.assertEqual(all(vstack_12_channels[2] == array(clipped_rds_left_fcn(matrix_int_left_updated,current_ref))),True)
            self.assertEqual(all(vstack_12_channels[3] == array(clipped_rds_right_fcn(matrix_int_right_updated,current_ref))),True)


if __name__ == '__main__':
    unittest.main()
