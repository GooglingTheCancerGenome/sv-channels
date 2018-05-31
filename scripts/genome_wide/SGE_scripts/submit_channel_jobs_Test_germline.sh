#!/bin/bash

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

#Add path to BAM files

BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
SAMPLE_ARRAY=('NA12878' 'NA12878' 'GIAB12878' 'PATIENT1' 'PATIENT2')

CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for i in ${#SAMPLE_ARRAY[@]}; do

    SAMPLE=${SAMPLE_ARRAY[$i]}
    BAM=${BAM_ARRAY[$i]}
    OUTDIR=$SAMPLE

    for CHROMOSOME in ${CHRARRAY[@]}; do
        for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do

            JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}

            echo qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
        done
    done
done


elif [ $RUNALL == 1 ]; then

# Output should be in the Tumor folder
i=0


# This BAM is only used to extract header information
BAM=${INPATH}"$SAMPLE""/BAM/"$SAMPLE"/mapping/"$SAMPLE"_dedup.bam"

# ChannelMaker script to generate channel data for Training data
PRG='channel_maker_real_germline'
for SAMPLE in ${SAMPLE_ARRAY[@]}; do
	for CHROMOSOME in ${CHRARRAY[@]}; do
		    OUTDIR=$SAMPLE
		    echo qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR} make_channel.sge
    done
done

fi