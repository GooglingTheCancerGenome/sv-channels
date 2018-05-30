#!/bin/bash

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real.py script

SVMODE='INV'

INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_'$SVMODE'/samples/'

# The Tumor/Normal sample pair for the somatic category
SOMATIC1='T0'
SOMATIC2='N0'
INARRAY_SOMATIC=(${INPATH}"$SOMATIC1""/BAM/"$SOMATIC1"/mapping/"$SOMATIC1"_dedup.bam" ${INPATH}"$SOMATIC2""/BAM/"$SOMATIC2"/mapping/"$SOMATIC2"_dedup.bam")

# The Tumor/Normal sample pair for the germline category
GERMLINE1='G1'
GERMLINE2='G2'
INARRAY_GERMLINE=(${INPATH}"$GERMLINE1""/BAM/"$GERMLINE1"/mapping/"$GERMLINE1"_dedup.bam" ${INPATH}"$GERMLINE2""/BAM/"$GERMLINE2"/mapping/"$GERMLINE2"_dedup.bam")

# The Tumor/Normal sample pair for the NoSV category
NOSV1="N1"
NOSV2="N2"
INARRAY_NOSV=(${INPATH}"$NOSV1""/BAM/"$NOSV1"/mapping/"$NOSV1"_dedup.bam" ${INPATH}"$NOSV2""/BAM/"$NOSV2"/mapping/"$NOSV2"_dedup.bam")

OUTARRAY=('Tumor' 'Normal')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for SAMPLE in somatic germline; do
	if [ $SAMPLE == 'somatic' ]; then
		INARRAY=(${INARRAY_SOMATIC[@]})
	elif [ $SAMPLE == 'germline' ]; then
		INARRAY=(${INARRAY_GERMLINE[@]})
	elif [ $SAMPLE == 'noSV' ]; then
                INARRAY=(${INARRAY_NOSV[@]})
	fi
	for i in 0 1; do
		for CHROMOSOME in 17; do
			for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do

				BAM=${INARRAY[$i]}
				OUTDIR="Training/"$SAMPLE"/"${OUTARRAY[$i]}
				qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR channel_script_by_name.sge
			done
		done
	done
done

elif [ $RUNALL == 1 ]; then

# Output should be in the Tumor folder
i=0

INARRAY=(${INARRAY_SOMATIC[@]})
# This BAM is only used to extract header information
BAM=${INARRAY[0]}

# ChannelMaker script to generate channel data for Training data
PRG='channel_maker_train'

for SAMPLE in somatic germline; do
	for CHROMOSOME in 17; do
		OUTDIR="Training/"$SAMPLE"/"${OUTARRAY[$i]}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},SVMODEARG=${SVMODE},OUTARG=${OUTDIR} channel_script_by_name.sge
	done
done

fi