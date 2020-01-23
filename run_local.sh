#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "2" ]; then
  echo "Usage: $0 [BAM file] [SEQID...]"
  exit 1
fi

# set variables
BAM=$(realpath -s $1)
BASE_DIR=$(dirname $BAM)
SAMPLE=$(basename $BAM .bam)
SEQ_IDS=${@:2}
TWOBIT=${BASE_DIR}/${SAMPLE}.2bit
BIGWIG=${BASE_DIR}/${SAMPLE}.bw
TSV=${BASE_DIR}/${SAMPLE}.tsv
BEDPE=${BASE_DIR}/${SAMPLE}.bedpe
WORK_DIR=scripts/genome_wide
NUMEXPR_MAX_THREADS=128  # required by py-bcolz

#source ~/.profile
cd $WORK_DIR
printenv

# convert the truth set in TSV to BEDPE file
# N.B.: consider only INS and DEL types
awk '{OFS="\t"}{if($5 ~ /DEL|INS/){print $1,$2,$2+1,$1,$4,$4+1,$5}}' \
  $TSV > $BEDPE

# write channels into *.json.gz and *.npy.gz files
for s in ${SEQ_IDS[@]}; do  # per chromosome
  p=clipped_read_distance && python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log
  p=clipped_reads && python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log
  p=clipped_read_pos && python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log
  p=split_reads && python $p.py -b $BAM -c $s -o $p.json.gz -p . -l $p.log
  p=snv && python $p.py -b $BAM -c $s -t $TWOBIT -o $p.npy -p . -l $p.log
  p=coverage && python $p.py -b $BAM -c $s -o $p.npy -p . -l $p.log
  p=chr_array && python $p.py -b $BAM -c $s -t $TWOBIT -m $BIGWIG -o $p.npy \
    -p . -l $p.log
  p=label_window_pairs_on_split_read_positions && python $p.py -b $BAM -c $s \
    -w 200 -gt $BEDPE -o $p.json.gz -p . -l $p.log
  p=label_window_pairs_on_svcallset && python $p.py -b $BAM -c $s \
    -w 200 -gt $BEDPE -o $p.json.gz -p . -l $p.log
done

echo -e "\nLog files:"
find -type f -name *.log
echo -e "\nOutput files:"
find -type f -name *.json.gz | grep "." || exit 1
find -type f -name *.npy.gz | grep "." || exit 1
