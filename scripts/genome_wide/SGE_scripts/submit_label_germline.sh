#!/usr/bin/env bash

OUTPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/'
PRG="make_label_germline"

for SAMPLE in "NA12878"; do
    for WINDOW in 200; do

		OUTDIR=$OUTPATH$SAMPLE
		JOB_NAME=$SAMPLE"_"$WINDOW"_"${PRG}
		qsub -v SAMPLEARG=${SAMPLE},WINDOWARG=${WINDOW},OUTARG=${OUTDIR} \
			-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" ${PRG}.sge

	done
done