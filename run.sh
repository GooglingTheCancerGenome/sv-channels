#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# != 2 ]; then
  echo "Usage: $0 [SCHEDULER {gridengine,slurm}] [BAM file]"
  exit 1
fi

printenv
source ~/.profile
cd genome_wide
conda env update --file environment.yaml

# set [env] variables
SCH=$1
BAM=$2
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)

export SAMPLEARG=$SAMPLE
export BAMARG=$BAM
export OUTARG=$BASE_DIR

xenon --version

# create channels per chromosome
for chn in ${CHANNELS[@]}; do
  export PRGARG=$chn
  for chr in ${CHROMS[@]}; do
    export CHRARG=$chr
    xenon -v scheduler $SCH --location local:// submit --name $SAMPLE_$chn \
      --cores-per-task 1 --inherit-env --max-run-time 1 --working-directory . \
      --stderr stderr-%j.log --stdout stdout-%j.log ./make_channel.sh
  done
done

# list channel outfiles (*.json.gz)
echo -e "\nOutput files:"
cd $BASE_DIR
#ls
find -type f -name \*.gz | grep "." || exit 1

# write stdout/stderr logs into terminal
echo -e "\nLog files:"
cd ..
#ls
for f in $(ls *.log); do
  echo "### $f ###"
  cat $f
done