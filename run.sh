#!/bin/bash -xe

source ~/.profile
export BRANCH=iss6
git clone -b $BRANCH https://github.com/GooglingTheCancerGenome/CNN.git
cd CNN/scripts/genome_wide
conda install --file environment.yaml

DATA_DIR=../../data/test
BAM=hmz-sv.bam
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)
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

echo -e "\nLog files:"
for f in $(find $DATA_DIR -name *.log); do
  echo -e "\n### $f ###\n"
  cat $f
done