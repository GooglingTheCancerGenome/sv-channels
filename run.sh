#!/bin/bash -xe

# clone repo & install deps into conda env
source ~/.profile
export BRANCH=iss6
git clone -b $BRANCH https://github.com/GooglingTheCancerGenome/CNN.git
cd CNN/scripts/genome_wide
conda install --file environment.yaml

# set env variables
SCH=$1
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
    xenon -vvv scheduler $SCH --location local:// submit --name $SAMPLE_$chn \
      --cores-per-task 1 --inherit-env --max-run-time 1 --working-directory . \
      --stderr stderr-%j.log --stdout stdout-%j.log make_channel.sh
  done
done

# write stdout/stderr logs into terminal
echo -e "\nLog files:"
for f in $(find $DATA_DIR -name \*.log); do
  echo -e "\n### $f ###\n"
  cat $f
done