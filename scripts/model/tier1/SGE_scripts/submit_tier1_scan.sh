#!/bin/bash -x

MODEL="/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/NA12878/cnn/train_NA12878_test_NA12878_cv1/models/results_model_2.hdf5"
OUTPATH='/hpc/cog_bioinf/ridder/users/lsantuari/Processed/DeepSV/channel_data/'

CHRARRAY=(`seq 1 22` 'X' 'Y')

PRG='scan_genome_with_model'

SAMPLE_ARRAY=('NA24385')

for (( i=0; i<${#SAMPLE_ARRAY[@]}; i++)); do
#for i in 0; do
	SAMPLE=${SAMPLE_ARRAY[$i]}

    # for CHROMOSOME in ${CHRARRAY[@]}; do
    for CHROMOSOME in '1'; do

        for WINDOW in 200; do

		    	OUTDIR=${OUTPATH}
		    	JOB_NAME=${SAMPLE}"_win"${WINDOW}"_"${PRG}
		    	qsub -v SAMPLE_ARG=${SAMPLE},CHRARG=${CHROMOSOME},MODELARG=${MODEL},PRGARG=${PRG},OUTDIRARG=${OUTDIR},CHANNELDIRARG=${OUTDIR},WINARG=${WINDOW} \
				    -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" tier1_chromosome_scan.sge

        done
    done
done
