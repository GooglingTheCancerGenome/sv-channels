#!/bin/bash -xe

source ~/.profile
export BRANCH=iss6
git clone -b $BRANCH https://github.com/GooglingTheCancerGenome/CNN.git
cd CNN/scripts/genome_wide
conda env create -n cm -f environment.yaml && conda activate cm

DATA_DIR=data
BAM=hmz-sv.bam

for PRG in (clipped_read_pos clipped_reads split_reads); do
  SAMPLEARG="$(basename $BAM .bam)"
  JOB_NAME="$SAMPLEARG_$PRG"
  export SAMPLEARG=$SAMPLE
  export BAMARG=$BAM
  export PRGARG=$PRG
  export OUTARG=$DATA_DIR
  ./make_channel.sge
done
