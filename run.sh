#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [SCHEDULER {local,gridengine,slurm}] [BAM file] [SEQID1 ...N]"
  exit 1
fi

# set variables
SCH=$1  # scheduler type
BAM="$(realpath -s "$2")"
BASE_DIR="$(dirname "$BAM")"
SAMPLE="$(basename "$BAM" .bam)"
SEQ_IDS=(${@:3})
SEQ_IDS_CSV="$(IFS=, ; echo "${SEQ_IDS[*]}")"  # stringify
SV_TYPES=(DEL)  # INS INV DUP CTX)
SV_CALLS=(split_reads gridss)  # manta delly lumpy)
CV_MODES=(kfold chrom)  # cross validation modes
KMERS=19
MAX_MISMATCH=2
KFOLD=2  # k-fold cross validation
EPOCHS=1  # epochs
WIN_SZ=25  # window size in bp
PREFIX="$BASE_DIR/$SAMPLE"
FASTA="$PREFIX.fasta"
TWOBIT="$PREFIX.2bit"
BIGWIG="$PREFIX.bw"
BEDPE="$PREFIX.bedpe"
BED="$PREFIX.bed"
EXCL_LIST="$BASE_DIR/ENCFF001TDO.bed"
REF_REG="$BASE_DIR/reference_N_regions.bed"
JOBS=()  # array of job IDs
JOBS_LOG=jobs.json  # job accounting log
RTIME=20  # runtime in minutes
SLTIME=1  # sleep X minutes
STIME=$(date +%s)
CONDA_ENV="sv-channels"  # conda env in gtcg/xenon-* docker images

# define functions
submit () {  # submit a job via Xenon CLI
  local xenon="xenon scheduler $SCH "
  local exec=$1
  local jobname=$2

  if [ "$SCH" == 'local' ]; then
    xenon+="exec --cores-per-task 1 "
  else
    xenon+="--location local:// submit --name '$jobname' --cores-per-task 1 \
      --stderr stderr-%j.log --stdout stdout-%j.log "
  fi

  xenon+="--inherit-env --max-run-time $RTIME --working-directory . "
  exec=$(echo $exec | sed 's/ / -- /')  # workaround argparse
  $xenon $exec
}

monitor () {  # monitor a job via Xenon CLI
  if [ "$SCH" == 'local' ]; then
    return
  fi

  xenon --json scheduler $SCH --location local:// list --identifier $1
}

waiting () {  # wait until all jobs are done
  if [ "$SCH" == 'local' ]; then
    return
  fi

  for j in "${JOBS[@]}"; do
    while true; do
      [[ $(monitor $j | grep -v "WARN" | jq '.statuses | .[] | select(.done==true)') ]] && \
        break || sleep ${SLTIME}m
    done
  done
}


# activate conda env
eval "$(conda shell.bash hook)"
conda activate $CONDA_ENV
conda list

# compute genome mappability
p=genmap
cmd="./$p.sh \"$FASTA\" \"$BIGWIG\" $KMERS $MAX_MISMATCH"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)

# convert SV calls (i.e. truth set and sv-callers output) in VCF to BEDPE files
cd scripts/R
p=vcf2bedpe
for vcf in $(find "$BASE_DIR" -mindepth 3 -name "*.vcf"); do
    prefix="$(basename "$vcf" .vcf)"
    bedpe="$BASE_DIR/$prefix.bedpe"
    cmd="./$p.R -i \"$vcf\" -o \"$bedpe\""
    JOB_ID=$(submit "$cmd" "$p-$prefix")
    JOBS+=($JOB_ID)
done

waiting

# submit jobs to output "channel" files (*.json.gz and *.npy.gz)
cd ../genome_wide
p=clipped_reads
cmd="python $p.py \
  -b \"$BAM\" \
  -c \"$SEQ_IDS_CSV\" \
  -o $p.json.gz \
  -p . \
  -l $p.log"
JOB_ID=$(submit "$cmd" $p)
JOBS+=($JOB_ID)

p=clipped_read_pos
cmd="python $p.py \
  -b \"$BAM\" \
  -c \"$SEQ_IDS_CSV\" \
  -o $p.json.gz \
  -p . \
  -l $p.log"
JOB_ID=$(submit "$cmd" $p)
JOBS+=($JOB_ID)

p=split_reads
cmd="python $p.py \
  -b \"$BAM\" \
  -c \"$SEQ_IDS_CSV\" \
  -o $p.json.gz \
  -ob $p.bedpe.gz \
  -p . \
  -l $p.log"
JOB_ID=$(submit "$cmd" $p)
JOBS+=($JOB_ID)

for s in "${SEQ_IDS[@]}"; do  # per chromosome
  p=clipped_read_distance
  cmd="python $p.py \
    -b \"$BAM\" \
    -c $s \
    -o $p.json.gz \
    -p . \
    -l $p.log"
  JOB_ID=$(submit "$cmd" "$p-$s")
  JOBS+=($JOB_ID)

  p=snv
  cmd="python $p.py \
    -b \"$BAM\" \
    -c $s \
    -t \"$TWOBIT\" \
    -o $p.npy \
    -p . \
    -l $p.log"
  JOB_ID=$(submit "$cmd" "$p-$s")
  JOBS+=($JOB_ID)

  p=coverage
  cmd="python $p.py \
    -b \"$BAM\" \
    -c $s \
    -o $p.npy \
    -p . \
    -l $p.log"
  JOB_ID=$(submit "$cmd" "$p-$s")
  JOBS+=($JOB_ID)
done

waiting

# generate chromosome arrays from the channels as well as label window pairs
for s in "${SEQ_IDS[@]}"; do
  p=chr_array
  cmd="python $p.py \
    -b \"$BAM\" \
    -c $s \
    -t \"$TWOBIT\" \
    -m \"$BIGWIG\" \
    -o $p.npy \
    -p . \
    -l $p.log"
  JOB_ID=$(submit "$cmd" "$p-$s")
  JOBS+=($JOB_ID)
done

waiting

# Create labels
for c in "${SV_CALLS[@]}"; do
    for sv in "${SV_TYPES[@]}"; do
        p=label_windows
        cmd="python $p.py \
          -b \"$BED\" \
          -c \"$SEQ_IDS_CSV\" \
          -w $WIN_SZ \
          -gt \"$BEDPE\" \
          -s $sv \
          -sv \"$c\" \
          -o labels.json.gz \
          -p . \
          -l $p.log"
        JOB_ID=$(submit "$cmd" "$p-$c-$sv")
        JOBS+=($JOB_ID)
    done
done

waiting

# Create windows
for c in "${SV_CALLS[@]}"; do
    for sv in "${SV_TYPES[@]}"; do
        p=create_window_pairs
        out="cnn/win$WIN_SZ/$c/windows/$sv"
        lb="$out/labels.json.gz"
        cmd="python $p.py \
          -b \"$BAM\" \
          -c \"$SEQ_IDS_CSV\" \
          -lb \"$lb\" -ca . \
          -w $WIN_SZ \
          -p \"$out\" \
          -l $p.log"
        JOB_ID=$(submit "$cmd" "$p-$c-$sv")
        JOBS+=($JOB_ID)
    done
done

waiting

# Add window channels
for c in "${SV_CALLS[@]}"; do
    for sv in "${SV_TYPES[@]}"; do
        p=add_win_channels
        out="cnn/win$WIN_SZ/$c/windows/$sv"
        prefix="$out/windows"
        infile="$prefix.npz"
        outfile="${prefix}_en.npz"
        log="${prefix}_en.log"
        cmd="python $p.py \
          -b \"$BAM\" \
          -w $WIN_SZ \
          -i \"$infile\" \
          -o \"$outfile\" \
          -l \"$log\""
        JOB_ID=$(submit "$cmd" "$p-$c-$sv")
        JOBS+=($JOB_ID)
    done
done

waiting

exit 0
# Train and test model
for sv in "${SV_TYPES[@]}"; do
    for c in "${SV_CALLS[@]}"; do
        for cv in "${CV_MODES[@]}"; do
            p=train
            out_dir="cnn/win$WIN_SZ/$c"
            train_dir="$out_dir/windows/$sv/windows_en.npz"
            cmd="python $p.py \
              --training_sample_name \"$SAMPLE\" \
              --training_windows \"$train_dir\" \
              --test_sample_name \"$SAMPLE\" \
              --test_windows \"$train_dir\" \
              -k $KFOLD \
              -e $EPOCHS \
              -p \"$out_dir\" \
              -s $sv \
              -cv $cv \
              -l $p.log"
            JOB_ID=$(submit "$cmd" "$p-$sv-$c")
            JOBS+=($JOB_ID)
        done
    done
done

waiting

# merge SV calls
cd ../R
p=merge_sv_calls
for c in "${SV_CALLS[@]}"; do
    for m in "${CV_MODES[@]}"; do
        win_dir="../genome_wide/cnn/win$WIN_SZ"
        sv_dir="$win_dir/$c/$m"
        bedpe="$win_dir/sv-channels.$c.$m.$SAMPLE"
        cmd="$p.R \
          -i \"$sv_dir\" \
          -f \"$EXCL_LIST\" \
          -n \"$REF_REG\" \
          -m $c \
          -o \"$bedpe\""
        JOB_ID=$(submit "$cmd" "$p-$c-$m")
        JOBS+=($JOB_ID)
    done
done

waiting

# convert BEDPE to VCF
cd ../utils
for c in "${SV_CALLS[@]}"; do
    for m in "${CV_MODES[@]}"; do
        p=bedpe_to_vcf
        win_dir="../genome_wide/cnn/win$WIN_SZ"
        prefix="$win_dir/sv-channels.$c.$m.$SAMPLE"
        bedpe="$prefix.bedpe"
        vcf="$prefix.vcf"
        cmd="python $p.py \
          -i \"$bedpe\" \
          -b \"$TWOBIT\" \
          -o \"$vcf\" \
          -s \"$SAMPLE\""
        JOB_ID=$(submit "$cmd" "$p-$c-$m")
        JOBS+=($JOB_ID)
    done
done

waiting

ETIME=$(date +%s)
echo "Processing ${#JOBS[@]} jobs took $((ETIME - STIME)) sec to complete."

# collect job accounting info
for j in "${JOBS[@]}"; do
  monitor $j >> $JOBS_LOG
done
cat $JOBS_LOG

# output logs in std{out,err}-[jobid].log
echo "----------"
echo "Log files:"
for f in $(find -type f -name "*.log"); do
  echo "### $f ###"
  cat $f
done

# list "channel" files
echo "-------------"
echo "Output files:"
#ls
find -type f -name "*.json.gz" | grep "." || exit 1
find -type f -name "*.npy.gz" | grep "." || exit 1

# exit with non-zero if there are failed jobs
[[ $(jq ".statuses | .[] | select(.done==true and .exitCode!=0)" $JOBS_LOG) ]] \
  && exit 1 || exit 0
