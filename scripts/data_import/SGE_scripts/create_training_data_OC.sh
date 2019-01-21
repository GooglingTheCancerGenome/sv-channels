#!/bin/bash -x

# This script generate the channel data for the OC dataset

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
    COMP_NAME=${arrIN[0]}"--"${arrIN[1]}
    #echo $COMP_NAME
    COMPARISONS_ARRAY+=($COMP_NAME)
done < "$filename"

CHRARRAY=(`seq 1 22` 'X' 'Y' 'MT')

for (( i=0; i<${#COMPARISONS_ARRAY[@]}; i++)); do

    SAMPLE=${COMPARISONS_ARRAY[$i]}
    JOB_NAME=$SAMPLE"_lab"

    qsub -v SAMPLEARG=$SAMPLE \
           -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_label.sge
    
done

