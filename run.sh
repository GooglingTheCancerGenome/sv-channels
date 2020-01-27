#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [SCHEDULER {gridengine,slurm}] [BAM file] [SEQID...]"
  exit 1
fi

# set variables
SCH=$1  # scheduler type
BAM=$(realpath -s $1)
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
SEQ_IDS=${@:3}
TWOBIT=${BASE_DIR}/${SAMPLE}.2bit
BIGWIG=${BASE_DIR}/${SAMPLE}.bw
WORK_DIR=scripts/genome_wide
RTIME=10  # runtime in minutes
STIME=1   # sleep X minutes
STARTTIME=$(date +%s)
LOG=xenon.log  # Xenon log file in JSON format
JOBS=()  # store jobIDs
NUMEXPR_MAX_THREADS=128  # required by py-bcolz


submit () {  # submit a job via Xenon CLI
  xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log "$1"
}

monitor () {  # monitor a job via Xenon CLI
  xenon -v --json scheduler $SCH --location local:// list --identifier $1
}

source ~/.profile
cd $WORK_DIR
xenon --version

# submit jobs to output "channel" files (*.json.gz and *.npy.gz)
for s in ${SEQ_IDS[@]}; do  # per chromosome
  p=clipped_read_distance && JOB="python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=clipped_reads && JOB="python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=clipped_read_pos && JOB="python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=split_reads && JOB="python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=snv && JOB="python $p.py -b $BAM -c $s -t $TWOBIT -o $p.npy -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=coverage && JOB="python $p.py -b $BAM -c $s -o $p.npy -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
done

# wait until the jobs are done
for j in ${JOBS[@]}; do
  while true; do
    [[ $(monitor $j | jq '.statuses | .[] | select(.done==true)') ]] && break || sleep ${STIME}m
  done
done

# generate chromosome arrays from the channels above
for s in ${SEQ_IDS[@]}; do
  p=chr_array && JOB="python $p.py -b $BAM -c $s -t $TWOBIT -m $BIGWIG \
    -o $p.npy -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=label_window_pairs_on_split_read_positions && JOB="python $p.py -b $BAM \
    -c $s -w 200 -gt $BEDPE -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=label_window_pairs_on_svcallset && JOB="python $p.py -b $BAM -c $s -w 200 \
    -gt $BEDPE -sv $BASE_DIR/gridss -o $p.json.gz -p . -l $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
done

# wait until the last job is done
while true; do
  [[ $(monitor ${JOBS[-1]} | jq '.statuses | .[] | select(.done==true)') ]] \
    && break || sleep ${STIME}m
done

# collect job accounting info
for j in ${JOBS[@]}; do
  monitor $j >> $LOG
done

ENDTIME=$(date +%s)
echo "Processing took $((ENDTIME - STARTTIME)) seconds to complete."

# output logs in std{out,err}-[jobid].log
echo "---------------"
echo -e "Log files:"
for f in $(find -type f -name \*.log); do
  echo "### $f ###"
  cat $f
done

# list "channel" files
echo "---------------"
echo -e "Output files:"
#ls
find -type f -name \*.json.gz | grep "." || exit 1
find -type f -name \*.npy.gz | grep "." || exit 1

# exit with non-zero if there are failed jobs
[ $(cat $LOG | jq '.statuses | .[] | select(.done==true and .exitCode!=0)') ] \
  && exit 1 || exit 0

