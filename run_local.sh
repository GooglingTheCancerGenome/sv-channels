#!/bin/bash -xe

# install deps into new conda env
ENV_NAME=cm
cd scripts/genome_wide
conda env create -n $ENV_NAME -f environment.yaml && conda activate

# set env variables
DATA_DIR=../../data/test
BAM=hmz-sv.bam
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)

# create channels per chromosome
for chn in ${CHANNELS[@]}; do
  SAMPLEARG=$(basename $BAM .bam)
  export SAMPLEARG=$SAMPLE
  export BAMARG=$DATA_DIR/$BAM
  export PRGARG=$chn
  export OUTARG=$DATA_DIR
  for chr in ${CHROMS[@]}; do
    export CHRARG=$chr
    ./make_channel.sh
  done
done

# print logs into STDOUT
echo -e "\nLog files:"
for f in $(find $DATA_DIR -name *.log); do
  echo -e "\n### $f ###\n"
  cat $f
done