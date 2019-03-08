#!/bin/bash

# This script generate the channel data for the germline SVs of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

SVMODE='INDEL'

# INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/WG/run_'$SVMODE'_500K/samples/'
INPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_'$SVMODE'/samples/'

# The germline category
GERMLINE1='G1'
GERMLINE2='G2'
# The NoSV category
NOSV1="N1"
NOSV2="N2"

SAMPLE_ARRAY=(${GERMLINE1} ${NOSV1})

CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for SAMPLE in ${SAMPLE_ARRAY[@]}; do
    for CHROMOSOME in ${CHRARRAY[@]}; do
        for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do

            BAM=${INPATH}"$SAMPLE""/BAM/"$SAMPLE"/mapping/"$SAMPLE"_dedup.bam"
            OUTDIR="Training_"$SVMODE"/"$SAMPLE
            JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}"_"$SVMODE
            echo qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
        done
    done
done


elif [ $RUNALL == 1 ]; then

# Output should be in the Tumor folder
i=0

# This BAM is only used to extract header information
# BAM=${INPATH}"$SAMPLE""/BAM/"$SAMPLE"/mapping/"$SAMPLE"_dedup.bam"

# ChannelMaker script to generate channel data for Training data
PRG='channel_maker_train_germline'
# for SAMPLE in ${SAMPLE_ARRAY[@]}; do

for SAMPLE in G1 N1; do

    if [ $SAMPLE == 'N1' ] || [ $SAMPLE == 'N2' ]; then
        PRG='channel_maker_real_germline'
    fi

    for CHROMOSOME in ${CHRARRAY[@]}; do
        #for CHROMOSOME in 1; do
            BAM=${INPATH}"$SAMPLE""/BAM/"$SAMPLE"/mapping/"$SAMPLE"_dedup.bam"
            OUTDIR="Training_"$SVMODE"/"$SAMPLE
            JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${SVMODE}"_"${PRG}

		    qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},SVMODEARG=${SVMODE},OUTARG=${OUTDIR} \
		    -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge

	done
done

fi