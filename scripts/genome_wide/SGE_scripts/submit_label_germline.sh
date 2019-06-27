#!/usr/bin/env bash

PRG="make_label_germline"

for SAMPLE in "NA12878"; do
    for WINDOW in 200; do

		OUTDIR=$SAMPLE
		JOB_NAME=$SAMPLE"_"$WINDOW"_"${PRG}
		qsub -v SAMPLEARG=${SAMPLE},WINDOWARG=${WINDOW} \
			-N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" ${PRG}.sge

	done
done