#!/bin/bash -x

# This script generate the channel data for the somatic and the germline categories of the Training data.
# The NoSV category is also added here, but for the moment it is generated using the channel_maker_real_somatic.py script

FILES="/hpc/cog_bioinf/ridder/users/akuzniar/HMF_data/bam/HGS_*/*.bam"
#echo $DIR_ARRAY

BAM_ARRAY=()
SAMPLE_ARRAY=()

for FILE in $FILES; do
	BAM_ARRAY+=($FILE)
	SAMPLE=`basename $FILE .bam`
	SAMPLE_ARRAY+=($SAMPLE) 
done

#Create comparisons array
filename="OC_comparisons.txt"
COMPARISONS_ARRAY=()
while IFS=' ' read -r line; do
    arrIN=(${line// /})
    COMP_NAME=${arrIN[0]}"_"${arrIN[1]}
    #echo $COMP_NAME
    COMPARISONS_ARRAY+=($COMP_NAME)
done < "$filename"


#printf '%s\n' "${BAM_ARRAY[@]}"
#printf '%s\n' "${SAMPLE_ARRAY[@]}"

SVMODE="INDEL"

#Add path to BAM files
NA12878_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/GIAB12878/mapping/GIAB12878_dedup.bam'
NA12892_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/Set4GIAB12892/mapping/Set4GIAB12892_dedup.bam'
NA12891_BAM='/hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/Set4GIAB12891/mapping/Set4GIAB12891_dedup.bam'
PATIENT1_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient1/Patient1.bam"
PATIENT2_BAM="/hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient2/Patient2.bam"

#BAM_ARRAY=(${NA12878_BAM} ${NA12892_BAM} ${NA12891_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'NA12892' 'NA12891' 'PATIENT1' 'PATIENT2')

#BAM_ARRAY=(${NA12878_BAM} ${PATIENT1_BAM} ${PATIENT2_BAM})
#SAMPLE_ARRAY=('NA12878' 'PATIENT1' 'PATIENT2')

CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')

# Run single channel scripts (0) or ChannelMaker (1)
RUNALL=1

if [ $RUNALL == 0 ]; then

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

    SAMPLE=${SAMPLE_ARRAY[$i]}
    BAM=${BAM_ARRAY[$i]}
    OUTDIR="OC/"$SAMPLE
    [ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"

    # LOGDIR=${SAMPLE}"/log"
    # [ ! -d "$LOGDIR" ] && mkdir -p "$LOGDIR"

    for CHROMOSOME in ${CHRARRAY[@]}; do
        for PRG in clipped_read_pos coverage clipped_read_distance clipped_reads split_read_distance; do
	#for PRG in clipped_reads; do
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


# This BAM is only used to extract header information

SLICE=("${COMPARISONS_ARRAY[@]:0:2}")

# ChannelMaker script to generate channel data for Training data
PRG='channel_maker_real_germline'
#for SAMPLE in ${SAMPLE_ARRAY[@]}; do

for (( i=0; i<${#COMPARISONS_ARRAY[@]}; i++)); do
#for i in 0; do

	SAMPLE=${COMPARISONS_ARRAY[$i]}
	BAM=${BAM_ARRAY[0]}

	#LOGDIR=$SAMPLE"/log"
	#echo "creating directory " $LOGDIR
	#[ ! -d ${LOGDIR} ] && mkdir -p $LOGDIR

	for CHROMOSOME in ${CHRARRAY[@]}; do
	#for CHROMOSOME in 2; do
		OUTDIR="OC_COMPARISONS/"$SAMPLE
		[ ! -d "$OUTDIR" ] && mkdir -p "$OUTDIR"

		JOB_NAME=$SAMPLE"_"$CHROMOSOME"_"${PRG}
		qsub -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},SVMODEARG=${SVMODE} \
			-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    	done
	#mv ${SAMPLE}"*.err" ${LOGDIR}
	#mv ${SAMPLE}"*.out" ${LOGDIR}
done


fi
