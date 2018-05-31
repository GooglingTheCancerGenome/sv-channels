#!/bin/bash

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

SVMODE='INDEL'

INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/WG/run_'$SVMODE'_500K/samples/'

# The germline category
GERMLINE1='G1'
GERMLINE2='G2'
# The NoSV category
NOSV1="N1"
NOSV2="N2"

INARRAY=(${INPATH}"$GERMLINE1""/BAM/"$GERMLINE1"/mapping/"$GERMLINE1"_dedup.bam" \
    ${INPATH}"$GERMLINE2""/BAM/"$GERMLINE2"/mapping/"$GERMLINE2"_dedup.bam" \
    ${INPATH}"$NOSV1""/BAM/"$NOSV1"/mapping/"$NOSV1"_dedup.bam" \
    ${INPATH}"$NOSV2""/BAM/"$NOSV2"/mapping/"$NOSV2"_dedup.bam")

SAMPLE_ARRAY=(${GERMLINE1} ${GERMLINE2} ${NOSV1} ${NOSV2})
OUTARRAY=('Sample')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for SAMPLE in ${SAMPLE_ARRAY[@]; do
	for i in 0; do
		for CHROMOSOME in 17; do
			for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do

				BAM=${INARRAY[$i]}
				OUTDIR="Training/"$SAMPLE"/"${OUTARRAY[$i]}
				JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}
				echo qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
				-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
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
PRG='channel_maker_train_germline'

for CHROMOSOME in `seq 1 22` 'X' 'Y' 'MT'; do
		OUTDIR="Training/"$SAMPLE"/"${OUTARRAY[$i]}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},SVMODEARG=${SVMODE},OUTARG=${OUTDIR} channel_script_by_name.sge
	done
done

fi