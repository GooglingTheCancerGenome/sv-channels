#!/usr/bin/env bash

GENOME=$1
BED=$2
BW=$3
CHRSIZES=$4
OUTDIR=$5

# conda activate sv-gen_utils

CMDDIR=$PWD
[ ! -d $OUTDIR ] && mkdir -p $OUTDIR
cd $OUTDIR

bedtools getfasta -fi $GENOME -bed $BED -fo test.tmp.fasta
awk '{if($1~/^>/){split($1,a,":"); print a[1]}else{print}}' test.tmp.fasta > test.fasta

faToTwoBit test.fasta test.2bit

python $CMDDIR/bigwig_from_bed.py --bigwig $BW \
    --bigwigout test.bw \
    --bed $BED \
    --chromsizes $CHRSIZES