#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -lt "3" ]; then
  echo "Usage: $0 [SCHEDULER {local,gridengine,slurm}] [BWA index file] [BAM file] [output directory]"
  exit 1
fi

MY_ENV=sv-channels
SCH=$1
BWA_INDEX=$2
BAM=$3
OUTDIR=$4
RTIME=1440

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


# activate conda env
eval "$(conda shell.bash hook)"
conda activate $MY_ENV
conda list

cmd="sh ./realign_clipped_reads.sh \"$BWA_INDEX\" \"$BAM\" \"$OUTDIR\""
JOB_ID=$(submit "$cmd" realign)
JOBS+=($JOB_ID)