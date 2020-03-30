#!/usr/bin/env bash

REF=$(realpath -s "$1")
REFNAME=$(basename "$REF" .fa.gz)

# extract the FASTA file:
gunzip $REF

OUTDIR=$2
[ ! -d $OUTDIR ] && mkdir -p $OUTDIR
cd $OUTDIR

# specify read length the GEM mappability track
READ_LENGTH=150

# prepare the GEM mappability track in BigWig format:
sh ../mappability/run_gem.sh $REF $REF.fa $READ_LENGTH .

# get chromosome sizes
samtools faidx ref.fa
cut -f1,2 $REF.fa > sizes.genome

GENOME=$REF.fa
BED=../data/seqs.bed
BW=$REF.bw
CHRSIZES=sizes.genome
OUTDIR=.

# extract chromosome regions
../create_test_fasta.sh $GENOME $BED $BW $CHRSIZES $OUTDIR