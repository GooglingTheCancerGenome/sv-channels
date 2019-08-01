#!/bin/bash -x

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

OUTPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/'

#Add path to BAM files
# NA12878_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/GIAB12878/mapping/GIAB12878_dedup.bam'
NA12878_BAM="/hpc/cog_bioinf/ridder/users/akuzniar/NA12878/bam/RMNISTHS_30xdownsample.bam"

NA24385_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/GiaB/HG002_NA24385_son/NIST_Illumina_2x250bps/bam/NA24385/mapping/NA24385_dedup.bam"

CHM1_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM1/bam/CHM1/mapping/CHM1_dedup.bam"
CHM13_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM13/bam/CHM13/mapping/CHM13_dedup.bam"
CHM1_CHM13_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CHM/CHM1_CHM13/bam/CHM1_CHM13_dedup.bam"

#BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA12892' 'NA12891' 'PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA12878_BAM})
#SAMPLE_ARRAY=('NA12878')

BAM_ARRAY=(${NA24385_BAM})
SAMPLE_ARRAY=('NA24385')

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

PRG='create_windows'

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do
#for i in 0; do

	SAMPLE=${SAMPLE_ARRAY[$i]}
	BAM=${BAM_ARRAY[$i]}

    for MODE in training set; do
    for WINDOW in 200; do

	    for CHROMOSOME in ${CHRARRAY[@]}; do
	    #for CHROMOSOME in 1; do
		    OUTDIR=$OUTPATH
		    JOB_NAME=$SAMPLE"_win"$WINDOW"_"$CHROMOSOME"_"${PRG}
		    qsub -v SAMPLEARG=$SAMPLE,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},MODEARG=${MODE},WINDOWARG=${WINDOW} \
			    -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_windows.sge
    done
    done

done
done