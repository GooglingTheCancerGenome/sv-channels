#!/usr/bin/env bash

#based on https://github.com/GooglingTheCancerGenome/CNN/blob/iss8/run.sh

set -xe

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [SCHEDULER {gridengine,slurm}] [Truth_set_name] [Truth_set_file]"
  exit 1
fi

# set [env] variables
SCH=$1  # scheduler type
TRUTHSET_NAME=$2
TRUTHSET_TRAINING=$(realpath $3)

USERDIR="/hpc/cog_bioinf/ridder/users/lsantuari"
CHANNELDIR=$USERDIR"/Processed/DeepSV/channel_data"

#Named after the truth set used
OUTPUTDIR=$USERDIR"/Processed/DeepSV/two_tier/"$TRUTHSET_NAME
[ ! -d $OUTPUTDIR ] && mkdir -p $OUTPUTDIR

SAMPLE_TRAINING="NA12878"

SAMPLE_TEST="NA24385"
TRUTHSET_TEST=$USERDIR"/Datasets/GiaB/HG002_NA24385_son/NIST_SVs_Integration_v0.6/processed/HG002_SVs_Tier1_v0.6.PASS.vcf.gz"

CHANNELDIR_TRAINING=$CHANNELDIR"/"$SAMPLE_TRAINING
CHANNELDIR_TEST=$CHANNELDIR"/"$SAMPLE_TEST

WINDOW=200
SHIFT=0

CHRARRAY=(`seq 1 22` 'X' 'Y')
CHRLIST=${CHRARRAY[@]}

WORK_DIR="/hpc/cog_bioinf/ridder/users/lsantuari/Git/DeepSV_refactoring/CNN/scripts/model/tier1"

RTIME=30  # runtime in minutes
STIME=1   # sleep X minutes
STARTTIME=$(date +%s)
LOG=xenon.log
JOBS=()  # store jobIDs
NUMEXPR_MAX_THREADS=128  # required by py-bcolz


submit () {  # submit a job via Xenon CLI
  xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log "$1"
}

monitor () {  # monitor a job via Xenon CLI
  xenon -v scheduler $SCH --location local:// list --identifier $1
}

check_jobs () {
    # check if all jobs are done
    for j in ${JOBS[@]}; do
      while true; do
        [ $(monitor $j | cut -f 5 | grep -i true) ] && break || sleep ${STIME}m
      done
    done
}

# source ~/.profile
cd $WORK_DIR
xenon --version

for TRAINSET in positive negative; do  # per training set

    TRAIN_OUTPUT_FILE=$OUTPUTDIR"/"$TRAINSET".npz"

    p=T0_S1_generate_training_data && JOB="python $p.py $TRAINSET \
    -chrlist $CHRLIST \
    -win $WINDOW \
    -truthset $TRUTHSET_TRAINING \
    -inputdir $CHANNELDIR_TRAINING \
    -output $TRAIN_OUTPUT_FILE"

    JOB_ID=$(submit "$JOB")
    JOBS+=($JOB_ID)

done

check_jobs

MODEL_FILE=$OUTPUTDIR"/model.h5"
POS_SET=$OUTPUTDIR"/positive.npz"
NEG_SET=$OUTPUTDIR"/negative.npz"

p=T0_S2_train && JOB="python $p.py \
    -positive $POS_SET \
    -negative $NEG_SET \
    -output $MODEL_FILE"
JOB_ID=$(submit "$JOB")
JOBS+=($JOB_ID)

check_jobs

# Running T0_S3_scan_chromosome

for CHROMOSOME in ${CHRARRAY[@]}; do
#CHROMOSOME=1

    #echo "Running chromosome "$CHROMOSOME"..."

    PRED_OUTPUT=$OUTPUTDIR"/"$CHROMOSOME"_predictions.npz"

    p=T0_S3_scan_chromosome && JOB="python $p.py \
        -inputdir $CHANNELDIR_TEST \
        -window $WINDOW \
        -chr $CHROMOSOME \
        -shift $SHIFT \
        -model $MODEL_FILE \
        -output $PRED_OUTPUT"
    JOB_ID=$(submit "$JOB")
    JOBS+=($JOB_ID)

done

check_jobs

# Running T0_S4_compare

CSV_OUTPUT=$OUTPUTDIR"/results.csv"
BED_OUTPUT=$OUTPUTDIR"/regions_of_interest.bed"

p=T0_S4_compare && JOB="python $p.py \
    -truthset $TRUTHSET_TEST \
    -chrlist $CHRLIST \
    -win $WINDOW \
    -inputdirlist $OUTPUTDIR \
    -output $CSV_OUTPUT \
    -outputbed $BED_OUTPUT"
JOB_ID=$(submit "$JOB")
JOBS+=($JOB_ID)

# check if the last job is done
watch -g -n 10 "xenon -v scheduler $SCH --location local:// list \
  --identifier ${JOBS[-1]} | cut -f 5 | grep -i true"

# collect job accounting info
for j in ${JOBS[@]}; do
  monitor $j >> $LOG
done

ENDTIME=$(date +%s)
echo "Processing took $((ENDTIME - STARTTIME)) seconds to complete."

# output channel/job logs in std{out,err}-[jobid].log
echo "---------------"
echo -e "Log files:"
for f in $(find -type f -name \*.log); do
  echo "### $f ###"
  cat $f
done

# list (channel) outfiles in *.json.gz and *.npy.gz
echo "---------------"
echo -e "Output files:"
#ls
find -type f -name *.npz | grep "." || exit 1
find -type f -name *.h5 | grep "." || exit 1

# check if there are failed jobs
[ $(grep -v "Exit code" $LOG | cut -f 7 | grep -v ^0) ] && exit 1
