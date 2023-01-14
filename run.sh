#!/usr/bin/env bash

set -e

# check input arg(s)
if [ $# -ne "1" ]; then
  echo "Usage: $0 [SCHEDULER]"
  exit 1
fi

# set variables
SCH=$1  # scheduler types: local, gridengine or slurm
JOBS=()  # array of job IDs
JOBS_LOG=jobs.json  # job accounting log
RTIME=5  # runtime in minutes
SLTIME=1  # sleep X minutes
STIME=$(date +%s)
CONDA_ENV="wf"

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

# Convert VCF files to BEDPE format
p=vcf2bedpe
cmd="Rscript svchannels/utils/R/$p.R -i data/test.vcf -o data/test.bedpe"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)

p=vcf2bedpe
cmd="Rscript svchannels/utils/R/$p.R -i data/vcf/manta_out/manta.vcf -o test/manta.bedpe"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)

# Extract signals
p=svchannels
cmd="$p extract-signals data/test.fasta data/test.bam -o test"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)
waiting

# Generate channels
p=svchannels
cmd="$p generate-channels --reference data/test.fasta test channels test/manta.bedpe"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)
waiting

# Label SVs
p=svchannels
cmd="$p label -f data/test.fasta.fai -o labels channels/sv_positions.bedpe data/test.bedpe"
JOB_ID=$(submit "$cmd" "$p")
JOBS+=($JOB_ID)
waiting

# Train the model
p=svchannels
cmd="$p train channels/channels.zarr.zip labels/labels.json.gz -m model.keras"
JOB_ID=$(submit "$cmd" "$p")
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

# exit with non-zero if there are failed jobs
[[ $(jq ".statuses | .[] | select(.done==true and .exitCode!=0)" $JOBS_LOG) ]] \
  && exit 1 || exit 0
