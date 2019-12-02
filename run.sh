#!/bin/bash -xe

# clone repo & install deps into current conda env
BRANCH=iss6
source ~/.profile
git clone -b $BRANCH https://github.com/GooglingTheCancerGenome/CNN.git
cd CNN/scripts/genome_wide
conda env update --file environment.yaml

# set [env] variables
SCH=$1
DATA_DIR=../../data/test
BAM=hmz-sv.bam
SAMPLE=$(basename $BAM .bam)
CHANNELS=(clipped_read_pos clipped_reads split_reads)
CHROMS=(17)

export SAMPLEARG=$SAMPLE
export BAMARG=$DATA_DIR/$BAM
export OUTARG=$DATA_DIR

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
cd $DATA_DIR
#ls
find . -name \*.gz

# write stdout/stderr logs into terminal
echo -e "\nLog files:"
cd ../../scripts/genome_wide
#ls
for f in $(find . -name \*.log); do
  echo "### $f ###"
  cat $f
done