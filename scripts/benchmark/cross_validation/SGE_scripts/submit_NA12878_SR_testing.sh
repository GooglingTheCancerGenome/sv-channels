#!/bin/bash -x

for RUN in {1..9}; do
# for RUN in 0; do
        JOB_NAME="CV_"$RUN
        qsub -v RUNARG=$RUN \
             -N $JOB_NAME \
             -o $JOB_NAME".out" \
             -e $JOB_NAME".err" run_NA12878_SR_testing_by_run.sge

done