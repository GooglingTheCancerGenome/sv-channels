#!/bin/bash -x

# This script generate the channel data for the GiaB Tumor/Normal mixture dataset

INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/GiaB/Synthetic_tumor/BAM/'
# Array with BAM input files for Tumor and Normal
INARRAY=(${INPATH}'Tumor/CPCT11111111T_dedup.realigned.bam' ${INPATH}'Reference/CPCT11111111R_dedup.realigned.bam')
# 'cr' files contains only clipped reads to speed up computation, as extracted with Picard (see 'run_Picard' script)
INARRAYCR=(${INPATH}'Tumor/CPCT11111111T_dedup.realigned.cr.bam' ${INPATH}'Reference/CPCT11111111R_dedup.realigned.cr.bam')
OUTARRAY=('Tumor' 'Normal')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then
# Tumor (0) and Normal (1) samples
for i in 0 1; do
	for CHROMOSOME in `seq 1 22` 'X'; do
		for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do
			BAM=${INARRAY[$i]}
			# For scripts clipped_reads.py and clipped_read_pos.py use clipped reads BAM files
			if [ $PRG == 'clipped_reads' ] || [ $PRG == 'clipped_read_pos' ]; then
				BAM=${INARRAYCR[$i]}	
			fi
			qsub -v CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTARRAY[$i]} channel_script_by_name.sge
		done
	done
done

elif [ $RUNALL == 1 ]; then

# Output should be in the Tumor folder
i=0
# This BAM is only used to extract header information
BAM=${INPATH}'Tumor/CPCT11111111T_dedup.realigned.cr.bam'

# ChannelMaker script to generate channel data for real data
PRG='channel_maker_real'
for CHROMOSOME in `seq 1 22` 'X'; do

	qsub -v CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTARRAY[$i]} channel_script_by_name.sge

done

fi
