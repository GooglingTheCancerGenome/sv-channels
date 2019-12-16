#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# != 2 ]; then
  echo "Usage: $0 [SCHEDULER {local,gridengine,slurm}] [BAM file]"
  exit 1
fi

# set [env] variables
RTIME=5
SCH=$1
BAM=$(realpath $2)
SEQ_IDS=${@:4}
SAMPLE=$(basename $BAM .bam)
TWOBIT=${BASE_DIR}/${SAMPLE}.2bit
WORK_DIR=scripts/genome_wide
JOBS=()  # store jobIDs

source ~/.profile
cd $WORK_DIR
printenv
xenon --version

# write channels into *.json.gz and *.npy.gz files
for p in clipped_read_pos clipped_reads split_reads; do  # calls per BAM
  CMD="python $p.py --bam $BAM --out $p.json.gz --outputpath ."
  JOB_ID=$(xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr.log --stdout stdout.log "$CMD")
  JOBS+=($JOB_ID)
done

for s in ${SEQ_IDS[@]}; do # calls per chromosome given BAM
  p=snv && CMD="python $p.py --bam $BAM --twobit $TWOBIT --chr $s --out $p.npy \
    --outputpath . --logfile $p.log"
  JOB_ID=$(xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr.log --stdout stdout.log "$CMD")
  JOBS+=($JOB_ID)
      
  p=coverage && CMD="python $p.py --bam $BAM --chr $s --out $p.npy \
    --outputpath . --logfile $p.log"
  JOB_ID=$(xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr.log --stdout stdout.log "$CMD")
  JOBS+=($JOB_ID)

  p=clipped_read_distance && python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log
  JOB_ID=$(xenon -v scheduler $SCH --location local:// submit \
    --name $SAMPLE_$p --cores-per-task 1 --inherit-env --max-run-time $RTIME \
    --working-directory . --stderr stderr.log --stdout stdout.log "$CMD")
  JOBS+=($JOB_ID)
done

# fetch job accounting info
sleep 10
for j in ${JOBS[@]}; do
   xenon -v scheduler $SCH --location local:// list --identifier $j
done

# write stdout/stderr logs into terminal
echo "---------------"
echo -e "Log files:"
for f in *.log; do
  echo "### $f ###"
  cat $f
done

# list channel outfiles (*.json.gz)
echo "---------------"
echo -e "Output files:"
#ls
find -type f -name *.json\* | grep "." || exit 1
find -type f -name *.npy\* | grep "." || exit 1
