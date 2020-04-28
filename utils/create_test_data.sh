#!/usr/bin/env bash

REF=$(realpath -s "$1")
REFNAME=$(basename "$REF" .fa)

# specify read length the GEM mappability track
READ_LENGTH=$2

OUTDIR=$3
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}
cd ${OUTDIR}

# prepare the GEM mappability track in BigWig format:
sh ../mappability/run_gem.sh ${REFNAME} ${REF} ${READ_LENGTH} .

# get chromosome sizes
samtools faidx -o ${REFNAME}.fai ${REF}
cut -f1,2 ${REFNAME}.fai > sizes.genome

GENOME=${REF}
BED=../data/seqs.bed
BW=${REFNAME}.${READ_LENGTH}mer.bw
CHRSIZES=sizes.genome
OUTDIR=.

# extract chromosome regions
sh ../create_test_fasta.sh ${GENOME} ${BED} ${BW} ${CHRSIZES} ${OUTDIR}