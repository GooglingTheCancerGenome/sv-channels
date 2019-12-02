#!/bin/bash -xe

cd scripts/genome_wide

# set penv] variables
DATA_DIR=../../data/test
BAM=hmz-sv.bam
SAMPLE=$(basename $BAM .bam)
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)

export SAMPLEARG=$SAMPLE
export BAMARG=$DATA_DIR/$BAM
export OUTARG=$DATA_DIR

# create channels per chromosome
for chn in ${CHANNELS[@]}; do
  export PRGARG=$chn
  for chr in ${CHROMS[@]}; do
    export CHRARG=$chr
    ./make_channel.sh
  done
done

cd $DATA_DIR
echo -e "\nOutput files:"
find . -name \*.gz
echo -e "\nLog files:"
find . -name \*.log