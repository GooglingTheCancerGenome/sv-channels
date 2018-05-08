#!/bin/bash -x
#$ -S /bin/bash
#$ -N ExtCR
#$ -l h_rt=04:00:00
#$ -l h_vmem=6G
#$ -q all.q
#$ -e ExtCR.err
#$ -o ExtCR.out
#$ -wd /hpc/cog_bioinf/ridder/users/lsantuari/Datasets/HMF/COLO829/BAM/HMF_COLO829_BAM/Reference
#$ -M mail
#$ -m beas

# Basename for the file
SAMPLENAME="COLO829R_dedup.realigned"
# Input BAM
INBAM="$SAMPLENAME"".bam"
# Output BAM
OUTBAM="$SAMPLENAME"".cr.bam"

# Load sambamba module
module load sambamcram/sambamba/0.6.5

# Extract left and right soft/hard clipped reads
sambamba view -F "cigar =~ /^\d+[S|H]/ or cigar =~ /[S|H]$/" -f bam "$INBAM" > "$OUTBAM"
# Generate bai index file
sambamba index "$OUTBAM" "$OUTBAM"".bai"
# Output flagstat informations
sambamba flagstat "$INBAM" > "$INBAM"".flagstat"
sambamba flagstat "$OUTBAM" > "$OUTBAM"".flagstat"
