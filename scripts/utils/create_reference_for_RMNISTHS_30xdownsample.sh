#!/bin/bash

# This bash script generates the reference needed for the BAM file:
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/NA12878/NIST_NA12878_HG001_HiSeq_300x/RMNISTHS_30xdownsample.bam

# Download reference hs37d5.fa.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.gz
wget ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/human_g1k_v37_decoy.fasta.fai.gz

# Remove the sequences NC_007605 and hs37d5 from the human_g1k_v37_decoy genome
wget https://raw.githubusercontent.com/GooglingTheCancerGenome/sv-callers/sandbox/scripts/remove_sequences_by_id.py
printf "NC_007605\nhs37d5\n" > seqs_to_remove.txt
python remove_sequences_by_id.py human_g1k_v37_decoy.fasta.gz seqs_to_remove.txt > human_g1k_v37_decoy_filtered.fasta

# Keep only the chromosome names in the FASTA headers
awk '{print $1}' human_g1k_v37_decoy_filtered.fasta > human_g1k_v37_decoy_filtered_short_header.fasta

# Download the NIST ERCC sequences
wget https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_Sequence_v1.FASTA

# extract the sequence with header 080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev from the file ercc_and_human_rRNA_and_tagdust.fa
seqtk subseq -l 60 ercc_and_human_rRNA_and_tagdust.fa seq.list > 080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev.fasta

# concatenate files together
cat human_g1k_v37_decoy_filtered_short_header.fasta SRM2374_Sequence_v1.FASTA 080418_Consensus_Vector_Sequence_NIST_SEQUENCING_ASSEMBLY_noRestrict_rev.fasta > b37_human_decoy_reference.dos.fasta

# convert DOS file into UNIX file
dos2unix -n b37_human_decoy_reference.dos.fasta b37_human_decoy_reference.fasta

