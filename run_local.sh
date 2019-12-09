#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# == 0 ]; then
  echo "Usage: $0 [BAM file]"
  exit 1
fi

# set [env] variables
BAM=$(realpath $1)
SAMPLE=$(basename $BAM .bam)
SEQ_IDS=(17.1 17.2)  # ($(seq 1 22) X Y MT)
WORK_DIR=scripts/genome_wide

source ~/.profile
cd $WORK_DIR
printenv
ls -alh

# output channels in *.json.gz files
for p in clipped_reads split_reads; do  # clipped_read_pos
  python $p.py --bam $BAM --out $p.json.gz --outputpath . \
    --logfile $p.log
done

for p in coverage clipped_read_distance; do
  for s in ${SEQ_IDS[@]}; do
    python $p.py --bam $BAM --chr $s --out $p.json.gz --outputpath . \
      --logfile $p.log
    done
done

ls -alh
echo -e "\nOutput files:"
find . -name \*.gz
echo -e "\nLog files:"
find . -name \*.log
