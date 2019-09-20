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

ART_INDEL_HET="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL_HET/samples/G1/BAM/G1/mapping/G1_dedup.bam"
ART_INDEL_HOM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/DeepSV/artificial_data/run_test_INDEL_HOM/samples/G1/BAM/G1/mapping/G1_dedup.bam"

#BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA12892' 'NA12891' 'PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA12878_BAM})
#SAMPLE_ARRAY=('NA12878')

BAM_ARRAY=(${ART_INDEL_HET}, ${ART_INDEL_HOM})
SAMPLE_ARRAY=('ART_INDEL_HET', 'ART_INDEL_HOM')

#BAM_ARRAY=(${NA24385_BAM})
#SAMPLE_ARRAY=('NA24385')

#BAM_ARRAY=(${CHM1_BAM} ${CHM13_BAM})
#SAMPLE_ARRAY=('CHM1' 'CHM13')

#BAM_ARRAY=(${NA12878_BAM} ${NA24385_BAM} ${CHM1_BAM} ${CHM13_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA24385' 'CHM1' 'CHM13')

# BAM_ARRAY=(${PATIENT1_BAM} ${PATIENT2_BAM})
# SAMPLE_ARRAY=('PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA24385_BAM} ${CHM1_CHM13_BAM})
#SAMPLE_ARRAY=('NA24385' 'CHM1_CHM13')

#CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')
#CHRARRAY=(`seq 1 22` 'X' 'Y')
CHRARRAY=('17')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

    SAMPLE=${SAMPLE_ARRAY[$i]}
    BAM=${BAM_ARRAY[$i]}
    OUTDIR=$OUTPATH$SAMPLE

    for PRG in clipped_read_pos clipped_reads split_reads; do
        JOB_NAME=$SAMPLE"_"${PRG}

        qsub -v SAMPLEARG=$SAMPLE,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

done

elif [ $RUNALL == 1 ]; then

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

    SAMPLE=${SAMPLE_ARRAY[$i]}
    BAM=${BAM_ARRAY[$i]}
    OUTDIR=$OUTPATH$SAMPLE

#    LOGDIR=${SAMPLE}"/log"
#    [ ! -d "$LOGDIR" ] && mkdir -p "$LOGDIR"

    for CHROMOSOME in ${CHRARRAY[@]}; do
    #for CHROMOSOME in '1'; do
        for PRG in coverage snv clipped_read_distance; do
	#for PRG in snv; do
            JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}

            qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            	-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
        done
    done
    #mv ${SAMPLE}"*.err" ${LOGDIR}
    #mv ${SAMPLE}"*.out" ${LOGDIR}
done


elif [ $RUNALL == 2 ]; then

# Output should be in the Tumor folder
i=0

SVMODE='INDEL'

# This BAM is only used to extract header information

SLICE=("${SAMPLE_ARRAY[@]:0:2}")

# ChannelMaker script to generate channel data for Training data
PRG='channel_maker_real_germline'
#for SAMPLE in ${SAMPLE_ARRAY[@]}; do

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do
#for i in 0; do

	SAMPLE=${SAMPLE_ARRAY[$i]}
	BAM=${BAM_ARRAY[$i]}

#	LOGDIR=$SAMPLE"/log"
#	echo "creating directory " $LOGDIR
#	[ ! -d ${LOGDIR} ] && mkdir -p $LOGDIR

for WINDOW in 200; do

	for CHROMOSOME in ${CHRARRAY[@]}; do
	#for CHROMOSOME in 1; do
		OUTDIR=$OUTPATH
		JOB_NAME=$SAMPLE"_win"$WINDOW"_"$CHROMOSOME"_"${PRG}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},SVMODEARG=${SVMODE},WINDOWARG=${WINDOW} \
			-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done
	#mv ${SAMPLE}"*.err" ${LOGDIR}
	#mv ${SAMPLE}"*.out" ${LOGDIR}

done
done

elif [ $RUNALL == 3 ]; then

PRG='chr_array'

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

	SAMPLE=${SAMPLE_ARRAY[$i]}
	BAM=${BAM_ARRAY[$i]}

    OUTPUTDIR=${OUTPATH}"/"$SAMPLE"/"$PRG
    [ ! -d $OUTPUTDIR ] && mkdir -p $OUTPUTDIR

	for CHROMOSOME in ${CHRARRAY[@]}; do
	#for CHROMOSOME in 1; do

		OUTDIR=$OUTPATH
		JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR} \
			-N $JOB_NAME -o $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".out" -e $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".err" make_channel.sge
    done

done

fi