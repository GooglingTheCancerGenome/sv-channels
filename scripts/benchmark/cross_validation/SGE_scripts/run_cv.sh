#!/bin/bash -x

for MODE in artificial real mixed; do
        JOB_NAME="CV_"$MODE
        qsub -v MODEARG=$MODE \
                -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_cv.sge

done