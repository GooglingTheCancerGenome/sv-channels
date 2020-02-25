# check_read
#
#  if read.reference_name == read.next_reference_name \
#    and read.mapping_quality >= minMAPQ \
#    and read.is_reverse != read.mate_is_reverse
#      
# n = 400569

SELECT *
FROM paired_primary_alignments_1
WHERE l_rname = r_rname -- l_rnext = '='
AND l_mapq >= 10 AND r_mapq >= 10
AND l_flag <> r_flag


# check_read_is_proper_paired_forward  # not proper pair -> forward discordant 
#
# if not read.is_unmapped \
#   and not read.mate_is_unmapped \
#   and read.mapping_quality >= minMAPQ \
#   and not read.is_proper_pair \
#   and not read.is_reverse -- ???
#
# n = 396

# check_read_is_proper_paired_reverse
#
# if not read.is_unmapped \
#   and not read.mate_is_unmapped \
#   and read.mapping_quality >= minMAPQ \
#   and not read.is_proper_pair and read.is_reverse -- ???
#
# n = 380

# read quality
# if not read.is_unmapped \
#   and read.mapping_quality >= minMAPQ
#
# n = 399958

================================
SELECT COUNT(*) FROM paired_primary_alignments_1;       # n = 199938
SELECT COUNT(*) FROM paired_secondary_alignments_1;     # n = 0
SELECT COUNT(*) FROM unpaired_alignments_1;             # n = 261

SELECT COUNT(*) FROM unpaired_primary_alignments_1;     # n = 399876
SELECT COUNT(*) FROM unpaired_secondary_alignments_1;   # n = 0
SELECT COUNT(*) FROM unpaired_all_alignments_1;         # n = 400137

SELECT COUNT(*) FROM unpaired_all_alignments_1 WHERE bam_flag(flag, 'segm_unma') = false; # n = 399972; see `is_unmapped`
SELECT COUNT(*) FROM unpaired_all_alignments_1 WHERE rnext = '=';  # n = 400137; see `[next_]referece_name`

\d unpaired_primary_alignments_1
create view bam.unpaired_primary_alignments_1 as 
 select l_virtual_offset as virtual_offset, qname, l_flag as flag, l_rname as rname, l_pos as pos, l_mapq as mapq, 
 l_cigar as cigar, l_rnext as rnext, l_pnext as pnext, l_tlen as tlen, l_seq as seq, l_qual as qual 
 from bam.paired_primary_alignments_1 
 union all 
 select r_virtual_offset as virtual_offset, qname, r_flag as flag, r_rname as rname, r_pos as pos, r_mapq as mapq, 
 r_cigar as cigar, r_rnext as rnext, r_pnext as pnext, r_tlen as tlen, r_seq as seq, r_qual as qual 
 from bam.paired_primary_alignments_1;

\d unpaired_secondary_alignments_1
create view bam.unpaired_secondary_alignments_1 as 
 select l_virtual_offset as virtual_offset, qname, l_flag as flag, l_rname as rname, l_pos as pos, l_mapq as mapq, 
 l_cigar as cigar, l_rnext as rnext, l_pnext as pnext, l_tlen as tlen, l_seq as seq, l_qual as qual 
 from bam.paired_secondary_alignments_1 
 union all 
 select r_virtual_offset as virtual_offset, qname, r_flag as flag, r_rname as rname, r_pos as pos, r_mapq as mapq, 
 r_cigar as cigar, r_rnext as rnext, r_pnext as pnext, r_tlen as tlen, r_seq as seq, r_qual as qual 
 from bam.paired_secondary_alignments_1;
 
 \d unpaired_all_alignments_1
create view bam.unpaired_all_alignments_1 as 
 select * 
 from bam.unpaired_primary_alignments_1 
 union all 
 select * 
 from bam.unpaired_secondary_alignments_1 
 union all 
 select * 
 from bam.unpaired_alignments_1;