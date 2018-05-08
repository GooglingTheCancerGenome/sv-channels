#!/bin/bash -x

# This script generate the channel data for the HMF COLO829 dataset

INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/HMF/COLO829/BAM/HMF_COLO829_BAM/'
# Array with BAM input files for Tumor (T) and Normal (R)
INARRAY=(${INPATH}'Tumor/COLO829T_dedup.realigned.bam' ${INPATH}'Reference/COLO829R_dedup.realigned.bam')
# 'cr' files contains only clipped reads to speed up computation, as extracted with Picard (see 'run_Picard' script)
INARRAYCR=(${INPATH}'Tumor/COLO829T_dedup.realigned.cr.bam' ${INPATH}'Reference/COLO829R_dedup.realigned.cr.bam')

SAMPLE="COLO829"
OUTARRAY=("$SAMPLE""/Tumor" "$SAMPLE""/Normal")

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then
# Tumor (0) and Normal (1) samples
for i in 0 1; do
	for CHROMOSOME in `seq 1 22` 'X'; do
		for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do
			BAM=${INARRAY[$i]}
			# With the scripts clipped_reads.py and clipped_read_pos.py, consider the BAM file with clipped reads (CR)
			if [ $PRG == 'clipped_reads' ] || [ $PRG == 'clipped_read_pos' ]; then
				BAM=${INARRAYCR[$i]}	
			fi
			qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTARRAY[$i]} channel_script_by_name.sge
		done
	done
done

elif [ $RUNALL == 1 ]; then

# Output should be in the Tumor folder
i=0
# This BAM is only used to extract header information
BAM=${INPATH}'Tumor/COLO829T_dedup.realigned.cr.bam'

# ChannelMaker script to generate channel data for real data
PRG='channel_maker_real'

for CHROMOSOME in `seq 1 22` 'X'; do

	qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTARRAY[$i]} channel_script_by_name.sge

done

fi