#!/usr/bin/env bash

set -xe

# check input arg(s)
if [ $# -ne "4" ]; then
  echo "Usage: $0 [FASTA] [BIGWIG] [k-mers] [# mismatches]"
  exit 1
fi

FASTA=$1
BIGWIG=$2
KMERS=$3
MAX_MISMATCH=$4
BASE_DIR="$(dirname "$FASTA")"
FAI="$FASTA.fai"
INDEX="$BASE_DIR/genmap_index"
MAP="$BASE_DIR/genmap_K${KMERS}_E${MAX_MISMATCH}"

rm -fr "$INDEX" "$MAP*" "$FAI*"
samtools faidx -o "$FAI" "$FASTA"
cut -f 1,2 "$FAI" > "$FAI.sizes"
genmap index -F "$FASTA" -I "$INDEX"
genmap map -K "$KMERS" -E "$MAX_MISMATCH" -I "$INDEX" -O "$MAP" -bg
bedGraphToBigWig "$MAP.bedgraph" "$FAI.sizes" "$BIGWIG"
