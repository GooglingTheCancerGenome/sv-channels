#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [SCHEDULER {gridengine,slurm}] [BAM file] [SEQID...]"
  exit 1
fi

# set [env] variables
SCH=$1  # scheduler type
BAM=$(realpath $2)
SEQ_IDS=${@:3}
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
TWOBIT=${BASE_DIR}/${SAMPLE}.2bit
BIGWIG=${BASE_DIR}/${SAMPLE}.bw
WORK_DIR=scripts/genome_wide
RTIME=10  # runtime in minutes
STIME=1   # sleep X minutes
LOG=xenon.log
JOBS=()  # store jobIDs


submit () {  # submit a job via Xenon CLI
  xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log "$1"
}

monitor () {  # monitor a job via Xenon CLI
  xenon -v scheduler $SCH --location local:// list --identifier $1
}

source ~/.profile
cd $WORK_DIR
printenv
xenon --version

# write channels into *.json.gz and *.npy.gz files
for s in ${SEQ_IDS[@]}; do  # per chromosome
  p=clipped_read_distance && JOB="python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=clipped_reads && JOB="python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
    
  p=split_reads && JOB="python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=clipped_read_pos && JOB="python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=snv && JOB="python $p.py --bam $BAM --chr $s --twobit $TWOBIT --out $p.npy \
   --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)

  p=coverage && JOB="python $p.py --bam $BAM --chr $s --out $p.npy \
    --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
done

# check if all jobs are done
for j in ${JOBS[@]}; do
  while true; do
    [ $(monitor $j | cut -f 5 | grep -i true) ] && break || sleep ${STIME}m
  done
done

# generate "chromosome arrays" from the channel files
for s in ${SEQ_IDS[@]}; do
  p=chr_array && JOB="python $p.py --bam $BAM --chr $s --twobit $TWOBIT \
    --map $BIGWIG --out $p.npy --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
done

# check if the last job is done
watch -g -n 10 "xenon -v scheduler $SCH --location local:// list \
  --identifier ${JOBS[-1]}) | cut -f 5 | grep -i true"

# collect job accounting info
for j in ${JOBS[@]}; do
  monitor $j >> $LOG
done

# write std{out,err} logs into terminal
echo "---------------"
echo -e "Log files:"
for f in *.log; do
  echo "### $f ###"
  cat $f
done

# check if there are failed jobs
[ $(grep -v "Exit code" $LOG | cut -f 7 | grep -v ^0) ] && exit 1

# list channel outfiles (*.json.gz)
echo "---------------"
echo -e "Output files:"
#ls
find -type f -name *.json.gz | grep "." || exit 1
find -type f -name *.npy.gz | grep "." || exit 1
