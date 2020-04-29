#!/usr/bin/env bash

set -xe

GENOME_NAME=$1
GENOME=$2
READ_LEN=$3
OUTDIR=$4
THREADS=$(nproc)

cd $OUTDIR

gem-indexer -i ${GENOME} -o ${GENOME_NAME}

gem-mappability -I ${GENOME_NAME}.gem -l ${READ_LEN} \
  -o ${GENOME_NAME}.${READ_LEN}mer -T ${THREADS}

gem-2-wig -I ${GENOME_NAME}.gem -i ${GENOME_NAME}.${READ_LEN}mer.mappability \
  -o ${GENOME_NAME}.${READ_LEN}mer

wigToBigWig ${GENOME_NAME}.${READ_LEN}mer.wig \
  ${GENOME_NAME}.${READ_LEN}mer.sizes ${GENOME_NAME}.${READ_LEN}mer.bw

