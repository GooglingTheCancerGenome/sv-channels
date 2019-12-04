#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# != 2 ]; then
  echo "Usage: $0 [SCHEDULER {gridengine,slurm}] [BAM file]"
  exit 1
fi

source ~/.profile
cd $HOME/scripts/genome_wide
printenv

# set [env] variables
SCH=$1
BAM=$2
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
PRGS=(clipped_read_pos clipped_reads split_reads)
#CHROMS=(17)

xenon --version

# create channels per chromosome
for prg in ${PRGS[@]}; do
  xenon -v scheduler $SCH --location local:// submit --name $SAMPLE_$prg \
    --cores-per-task 1 --inherit-env --max-run-time 1 --working-directory . \
    --stderr stderr-%j.log --stdout stdout-%j.log "python $prg.py --bam $BAM --out $prg.json.gz --outputpath ."
done

sleep 10

# write stdout/stderr logs into terminal
echo -e "\nLog files:"
for f in *.log; do
  echo "### $f ###"
  cat $f
done

# list channel outfiles (*.json.gz)
echo -e "\nOutput files:"
#ls
find -type f -name *.json.gz | grep "." || exit 1