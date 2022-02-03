#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -ne "4" ]; then
  echo "Usage: $0 [SCHEDULER] [BAM file] [SEQIDS] [SVTYPES]"
  exit 1
fi

# define variables

# set variables
SCH=$1  # scheduler types: local, gridengine or slurm
BAM="$(realpath -s "$2")"
BASE_DIR="$(dirname "$BAM")"
SAMPLE="$(basename "$BAM" .bam)"
PREFIX="$BASE_DIR/$SAMPLE"
FASTA="$PREFIX.fasta"
EXPAND=250
GAP=10
SIGNALS_DIR="sv-channels"
SVCALLER="manta"

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

# convert SV calls (i.e. truth set and sv-callers output) in VCF to BEDPE files
cd scripts/R
p=vcf2bedpe
for vcf in $(find "$BASE_DIR" -mindepth 1 -name "*.vcf" | grep -v "htz-sv*"); do
    prefix="$(basename "$vcf" .vcf)"
    bedpe="$BASE_DIR/$prefix.bedpe"
    cmd="./$p.R -i \"$vcf\" -o \"$bedpe\""
    JOB_ID=$(submit "$cmd" "$p-$prefix")
    JOBS+=($JOB_ID)
done

waiting

# extract SV signals
cd ././../svchannels
p=extract_signals
cmd="./$p.py --out-dir \"$SIGNALS_DIR\" \"$FASTA\" \"$BAM\""
JOB_ID=$(submit "$cmd" "$p-$prefix")
JOBS+=($JOB_ID)

waiting

# extract SV signals
cd ././../svchannels
p=generate_channels
bedpe="$BASE_DIR/$SVCALLER.bedpe"
cmd="./$p.py \"$SIGNALS_DIR\" \"$bedpe\" --expand \"$EXPAND\" --gap \"$GAP\""
JOB_ID=$(submit "$cmd" "$p-$prefix")
JOBS+=($JOB_ID)

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
ls
# find -type f -name "*.json.gz" | grep "." || exit 1


# exit with non-zero if there are failed jobs
[[ $(jq ".statuses | .[] | select(.done==true and .exitCode!=0)" $JOBS_LOG) ]] \
  && exit 1 || exit 0
