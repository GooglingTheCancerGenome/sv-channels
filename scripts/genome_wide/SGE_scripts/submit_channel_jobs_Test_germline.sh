#!/bin/bash -x

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

#Add path to BAM files
NA12878_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/GIAB12878/mapping/GIAB12878_dedup.bam'
NA12892_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/Set4GIAB12892/mapping/Set4GIAB12892_dedup.bam'
NA12891_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/Set4GIAB12891/mapping/Set4GIAB12891_dedup.bam'
PATIENT1_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient1/Patient1.bam"
PATIENT2_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient2/Patient2.bam"

#BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA12892' 'NA12891' 'PATIENT1' 'PATIENT2')

BAM_ARRAY=(${NA12878_BAM})
SAMPLE_ARRAY=('NA12878')

# BAM_ARRAY=(${PATIENT1_BAM} ${PATIENT2_BAM})
# SAMPLE_ARRAY=('PATIENT1' 'PATIENT2')

CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=0

if [ $RUNALL == 0 ]; then

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

    SAMPLE=${SAMPLE_ARRAY[$i]}
    BAM=${BAM_ARRAY[$i]}
    OUTDIR=$SAMPLE

#    LOGDIR=${SAMPLE}"/log"
#    [ ! -d "$LOGDIR" ] && mkdir -p "$LOGDIR"

    # for CHROMOSOME in ${CHRARRAY[@]}; do
    for CHROMOSOME in '1'; do
        #for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do
	for PRG in snv; do
            JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}

            qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            	-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
        done
    done
    #mv ${SAMPLE}"*.err" ${LOGDIR}
    #mv ${SAMPLE}"*.out" ${LOGDIR}
done


elif [ $RUNALL == 1 ]; then

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

	for CHROMOSOME in ${CHRARRAY[@]}; do
	#for CHROMOSOME in 1; do
		OUTDIR=$SAMPLE
		JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},SVMODEARG=${SVMODE} \
			-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    	done
	#mv ${SAMPLE}"*.err" ${LOGDIR}
	#mv ${SAMPLE}"*.out" ${LOGDIR}
done


fi