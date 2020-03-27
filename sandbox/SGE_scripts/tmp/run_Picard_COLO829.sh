#!/bin/bash -x
#$ -S /bin/bash
#$ -N Picard
#$ -l h_rt=04:00:00
#$ -l h_vmem=16G
#$ -q all.q
#$ -e picard.err
#$ -o picard.out
#$ -cwd
#$ -M mail
#$ -m beas

# COLO829 Tumor 'COLO829T' and Normal 'COLO829R' samples
INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/HMF/COLO829/BAM/HMF_COLO829_BAM/'
INARRAY=(${INPATH}'Tumor/COLO829T_dedup.realigned.bam' ${INPATH}'Reference/COLO829R_dedup.realigned.bam')

for BAMFILE in "${INARRAY[@]}"
do

# Get path and filename
FILEPATH="$(dirname $BAMFILE)"
FILENAME=$(basename "$BAMFILE")

# Get basename
#echo $FILENAME
arrFILE=(${FILENAME//_/ })
BASENAME="${arrFILE[0]}"
#echo $BASENAME

# Create picard output folder if it does not exist
PICARDOUTDIR="$FILEPATH/picard/"
if [ ! -d "$PICARDOUTDIR" ]; then
 mkdir $PICARDOUTDIR      
fi

# Load required modules
module load picard/2.1.0
module load R/3.4.1

# Run Picard tools with CollectInsertSizeMetrics
java -Xmx12g -jar $CLASSPATH/picard.jar CollectInsertSizeMetrics \
	I= $BAMFILE\
	O="$PICARDOUTDIR$BASENAME""_insert_size_metrics.txt" \
	H="$PICARDOUTDIR$BASENAME""_insert_size_histogram.pdf" \
	M=0.5
done
