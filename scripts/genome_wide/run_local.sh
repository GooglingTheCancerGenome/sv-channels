#!/bin/bash -xe

# set [env] variables
BAM=$1
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)

export SAMPLEARG=$SAMPLE
export BAMARG=$BAM
export OUTARG=$BASE_DIR

# create channels per chromosome
for chn in ${CHANNELS[@]}; do
  export PRGARG=$chn
  for chr in ${CHROMS[@]}; do
    export CHRARG=$chr
    ./make_channel.sh
  done
done

cd $BASE_DIR
echo -e "\nOutput files:"
find . -name \*.gz
echo -e "\nLog files:"
find . -name \*.log
