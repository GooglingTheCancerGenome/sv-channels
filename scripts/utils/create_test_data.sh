#!/usr/bin/env bash

set -xe

REF=$(realpath -s "$1")
REFNAME=$(basename "$REF" .fa)
READ_LEN=$2  # read length of the GEM mappability track
OUTDIR=$3
GENOME=${REF}
BED=../data/seqs.bed
BW=${REFNAME}.${READ_LEN}mer.bw
CHRSIZES=chr_sizes.txt

mkdir -p ${OUTDIR}
cd ${OUTDIR}

# prepare the GEM mappability track in BigWig format
../mappability/run_gem.sh ${REFNAME} ${REF} ${READ_LEN} .

# get chromosome sizes
samtools faidx -o ${REFNAME}.fai ${REF}
cut -f1,2 ${REFNAME}.fai > ${CHRSIZES}

# extract chromosome regions
../create_test_fasta.sh ${GENOME} ${BED} ${BW} ${CHRSIZES} ${OUTDIR}
