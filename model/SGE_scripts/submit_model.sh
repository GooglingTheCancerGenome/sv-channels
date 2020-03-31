#!/bin/bash -x

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

OUTPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/'

#Add path to BAM files
# NA12878_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/GIAB12878/mapping/GIAB12878_dedup.bam'
NA12878_BAM="/hpc/cog_bioinf/ridder/users/akuzniar/NA12878/bam/RMNISTHS_30xdownsample.bam"

NA24385_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/GiaB/HG002_NA24385_son/NIST_Illumina_2x250bps/bam/NA24385/mapping/NA24385_dedup.bam"

CHM1_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM1/bam_GRCh38/CHM1/mapping/CHM1_dedup.bam"
CHM13_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM13/bam_GRCh38/CHM13/mapping/CHM13_dedup.bam"
CHM1_CHM13_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM1_CHM13/bam_GRCh38/CHM1_CHM13_dedup.bam"

#BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA12892' 'NA12891' 'PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA24385_BAM})
#SAMPLE_ARRAY=('NA24385')

BAM_ARRAY=(${NA12878_BAM} ${NA24385_BAM} ${CHM1_CHM13_BAM})
SAMPLE_ARRAY=('NA12878' 'NA24385' 'CHM1_CHM13')

#BAM_ARRAY=(${NA12878_BAM} ${NA24385_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA24385')

#BAM_ARRAY=(${CHM1_CHM13_BAM})
#SAMPLE_ARRAY=('CHM1_CHM13')

#BAM_ARRAY=(${NA12878_BAM} ${NA24385_BAM} ${CHM1_CHM13_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA24385' 'CHM1_CHM13')

#BAM_ARRAY=(${CHM1_BAM} ${CHM13_BAM})
#SAMPLE_ARRAY=('CHM1' 'CHM13')

#BAM_ARRAY=(${NA12878_BAM} ${NA24385_BAM} ${CHM1_BAM} ${CHM13_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA24385' 'CHM1' 'CHM13')

# BAM_ARRAY=(${PATIENT1_BAM} ${PATIENT2_BAM})
# SAMPLE_ARRAY=('PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA24385_BAM} ${CHM1_CHM13_BAM})
#SAMPLE_ARRAY=('NA24385' 'CHM1_CHM13')

#CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')
#CHRARRAY=(`seq 1 22` 'X')

PRG='train_model_with_fit'

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do
#for i in 0; do
	TRAINING_SAMPLE=${SAMPLE_ARRAY[$i]}
	
	for (( j=0; j<${#SAMPLE_ARRAY[@]}; j++)); do

        TEST_SAMPLE=${SAMPLE_ARRAY[$j]}
	
	if [ $TRAINING_SAMPLE != $TEST_SAMPLE ]; then
    		for MODE in "test"; do

        		for WINDOW in 200; do

		    	OUTDIR=$OUTPATH
		    	JOB_NAME=$TRAINING_SAMPLE"_"$TEST_SAMPLE"_win"$WINDOW"_"$MODE"_"${PRG}
		    	qsub -v TRAINING_SAMPLE_ARG=$TRAINING_SAMPLE,TEST_SAMPLE_ARG=$TEST_SAMPLE,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},MODEARG=${MODE},WINDOWARG=${WINDOW} \
				    -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" train_model.sge

        		done
    		done
	fi
	done
done
