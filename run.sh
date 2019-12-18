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
RTIME=5  # runtime in minutes
LOG=xenon.log
JOBS=()  # store jobIDs


submit () {  # submit job via Xenon
  xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr-%j.log --stdout stdout-%j.log "$1"
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
    
  p=chr_array && JOB="python $p.py --bam $BAM --chr $s --twobit $TWOBIT \
   --map $BIGWIG --out $p.npy --outputpath . --logfile $p.log"
  JOB_ID=$(submit "$JOB")
  JOBS+=($JOB_ID)
done

# collect job accounting info
sleep 60
for j in ${JOBS[@]}; do
   xenon -v scheduler $SCH --location local:// list --identifier $j >> $LOG
done

# write stdout/stderr logs into terminal
echo "---------------"
echo -e "Log files:"
for f in *.log; do
  echo "### $f ###"
  cat $f
done

# check if there are failed jobs
cat $LOG
[ $(grep -v "Exit code" $LOG | cut -f 7 | grep -v ^0) ] && exit 1

# list channel outfiles (*.json.gz)
echo "---------------"
echo -e "Output files:"
#ls
find -type f -name *.json.gz | grep "." || exit 1
find -type f -name *.npy.gz | grep "." || exit 1
