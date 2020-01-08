#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "2" ]; then
  echo "Usage: $0 [BAM file] [SEQID...]"
  exit 1
fi

# set [env] variables
BAM=$(realpath $1)
SEQ_IDS=${@:2}
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
TWOBIT=${BASE_DIR}/${SAMPLE}.2bit
BIGWIG=${BASE_DIR}/${SAMPLE}.bw
WORK_DIR=scripts/genome_wide
NUMEXPR_MAX_THREADS=128  # required by py-bcolz

source ~/.profile
cd $WORK_DIR
printenv

# write channels into *.json.gz and *.npy.gz files
for s in ${SEQ_IDS[@]}; do  # per chromosome
  p=clipped_read_distance && python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log
  p=clipped_reads && python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log
  p=split_reads && python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log
  p=clipped_read_pos && python $p.py --bam $BAM --chr $s --out $p.json.gz \
    --outputpath . --logfile $p.log
  p=snv && python $p.py --bam $BAM --chr $s --twobit $TWOBIT --out $p.npy \
    --outputpath . --logfile $p.log
  p=coverage && python $p.py --bam $BAM --chr $s --out $p.npy \
    --outputpath . --logfile $p.log
  p=chr_array && python $p.py --bam $BAM --chr $s --twobit $TWOBIT \
    --map $BIGWIG --out $p.npy --outputpath . --logfile $p.log
done

echo -e "\nLog files:"
find . -name \*.log
echo -e "\nOutput files:"
find -type f -name *.json.gz | grep "." || exit 1
find -type f -name *.npy.gz | grep "." || exit 1
