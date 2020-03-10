#!/bin/bash -x

SAMPLE=$1
BAM=$2
BAM_SV=$3
OUTPATH=$4

CHRARRAY=(`seq 1 22` 'X' 'Y')

OUTDIR=$OUTPATH$SAMPLE

for PRG in clipped_read_pos clipped_reads split_reads; do

    JOB_NAME=$SAMPLE"_channels"

    qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,BAMARG=$BAM_SV,PRGARG=${PRG},OUTARG=$OUTDIR \
        -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
done

for CHROMOSOME in ${CHRARRAY[@]}; do

    for PRG in coverage snv; do

        JOB_NAME=$SAMPLE"_channels"

        qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

    for PRG in clipped_read_distance; do

        JOB_NAME=$SAMPLE"_channels"

        qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM_SV,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

done

PRG='chr_array'

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do

	SAMPLE=${SAMPLE_ARRAY[$i]}
	BAM=${BAM_ARRAY[$i]}

    OUTPUTDIR=${OUTPATH}"/"$SAMPLE"/"$PRG
    [ ! -d $OUTPUTDIR ] && mkdir -p $OUTPUTDIR

	for CHROMOSOME in ${CHRARRAY[@]}; do
	#for CHROMOSOME in 1; do

		OUTDIR=$OUTPATH
		JOB_NAME=$SAMPLE"_carray"
		JOB_NAME_HOLD=$SAMPLE"_channels"

		qsub -wd $OUTDIR -hold_jid $JOB_NAME_HOLD -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR} \
			-N $JOB_NAME -o $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".out" -e $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".err" make_channel.sge
    done

done