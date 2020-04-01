#!/bin/bash -x

GENOME_NAME=$1
GENOME=$2
READ_LENGTH=$3
OUTDIR=$4

THREADS=8

cd $OUTDIR

gem-indexer -i ${GENOME} -o ${GENOME_NAME}

gem-mappability -I ${GENOME_NAME}.gem -l ${READ_LENGTH} -o ${GENOME_NAME}.${READ_LENGTH}mer -T ${THREADS}

gem-2-wig -I ${GENOME_NAME}.gem -i ${GENOME_NAME}.${READ_LENGTH}mer.mappability -o ${GENOME_NAME}.${READ_LENGTH}mer

wigToBigWig ${GENOME_NAME}.${READ_LENGTH}mer.wig ${GENOME_NAME}.${READ_LENGTH}mer.sizes ${GENOME_NAME}.${READ_LENGTH}mer.bw