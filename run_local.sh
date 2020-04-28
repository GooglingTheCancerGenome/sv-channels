#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "2" ]; then
  echo "Usage: $0 [BAM file] [SEQID...]"
  exit 1
fi

# set variables
BAM=$(realpath -s "$1")
BASE_DIR=$(dirname "$BAM")
SAMPLE=$(basename "$BAM" .bam)
SV_TYPES=(DEL)
SV_CALLS=(gridss)  # to speed up exclude callers: manta delly lumpy
KFOLD=2            # k-fold cross validation
WIN_SZ=200  # in bp
SEQ_IDS=(${@:2})
SEQ_IDS_CSV=$(IFS=, ; echo "${SEQ_IDS[*]}")  # stringify
PREFIX="${BASE_DIR}/${SAMPLE}"
TWOBIT="${PREFIX}.2bit"
BIGWIG="${PREFIX}.bw"
BEDPE="${PREFIX}.bedpe"  # truth set
BED="${PREFIX}.bed" # chromosome regions
WORK_DIR=scripts/genome_wide

# convert SV calls (i.e. truth set and sv-callers output) in VCF to BEDPE files
for vcf in $(find data -name "*.vcf" | grep -E "test"); do
  prefix=$(basename $vcf .vcf)
  bedpe="${BASE_DIR}/${prefix}.bedpe"
  scripts/R/vcf2bedpe.R -i "$vcf" -o "$bedpe"
done

cd $WORK_DIR
printenv

# run per BAM file
p=clipped_reads
python $p.py -b "$BAM" -c "${SEQ_IDS_CSV}" -o $p.json.gz -p . -l $p.log

p=clipped_read_pos
python $p.py -b "$BAM" -c "${SEQ_IDS_CSV}" -o $p.json.gz -p . -l $p.log


p=split_reads
python $p.py -b "$BAM" -c "${SEQ_IDS_CSV}" -o $p.json.gz -ob $p.bedpe.gz -p . -l $p.log


# write channels into *.json.gz and *.npy.gz files
for s in "${SEQ_IDS[@]}"; do  # per chromosome

  p=clipped_read_distance
  python $p.py -b "$BAM" -c $s -o $p.json.gz -p . -l $p.log

  p=snv
  python $p.py -b "$BAM" -c $s -t "$TWOBIT" -o $p.npy -p . -l $p.log

  p=coverage
  python $p.py -b "$BAM" -c $s -o $p.npy -p . -l $p.log

  p=chr_array
  python $p.py -b "$BAM" -c $s -t "$TWOBIT" -m "$BIGWIG" -o $p.npy -p . -l $p.log

done

for sv in "${SV_TYPES[@]}"; do

    for c in "${SV_CALLS[@]}"; do
        p=label_windows
        python $p.py -b "${BED}" -w $WIN_SZ -c "${SEQ_IDS_CSV}" -gt "$BEDPE" \
          -s $sv -sv "$BASE_DIR/$c" -o labels.json.gz -p . -l $p.log

        p=create_window_pairs
        out="labels/win$WIN_SZ/$sv/$c"
        lb="$out/labels.json.gz"
        python $p.py -b "$BAM" -c "${SEQ_IDS_CSV}" -lb "$lb" -ca . -w $WIN_SZ \
          -p "$out" -l $p.log

        p=add_win_channels
        pfix="$out/windows/windows"
        iwin="${pfix}.npz"
        owin="${pfix}_en.npz"
        log="${pfix}_en.log"
        python $p.py -b "$BAM" -w "$WIN_SZ" -i ${iwin} -o ${owin} -l $log
        mv ${iwin} ${iwin}.bkup
        mv ${owin} ${iwin}

        p=train_model_with_fit
        python $p.py --training_sample_name ${SAMPLE} --training_sample_folder . \
          --test_sample_name ${SAMPLE} --test_sample_folder . -k $KFOLD -p "$out" \
          -s $sv -l $p.log
    done
done

echo -e "\nLog files:"
find -type f -name "*.log"
echo -e "\nOutput files:"
find -type f -name "*.json.gz" | grep "." || exit 1
find -type f -name "*.npy.gz" | grep "." || exit 1
