#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=32G

BAM=$(realpath -s "$1")
SAMPLE=$(basename "$BAM" .bam)
OUTDIR=$2
TMPDIR=$OUTDIR"/"$SAMPLE"/tmp"

# create tmp folder
[ ! -d $TMPDIR ] && mkdir -p $TMPDIR

echo activating...
conda activate samtools

echo clipped to split reads...
python split_clipped_reads.py \
    -b $BAM \
    -ob no_clipped.bam \
    -of clipped.fq \
    -l split_clipped_reads.log

echo sorting BAM without clipped reads...
samtools sort -o no_clipped_sorted.bam no_clipped.bam

echo indexing BAM without clipped reads...
samtools index no_clipped_sorted.bam

echo mapping soft clipped sequences...
bwa mem -L 0,0 ../../data/test.fasta clipped.fq > clipped.sam

echo converting SAM to BAM...
samtools view -b clipped.sam -o clipped.bam

echo sorting BAM with split reads...
samtools sort -o clipped_sorted.bam clipped.bam

echo indexing BAM with split reads..
samtools index clipped_sorted.bam

OUTBAM=$OUTDIR"/"$SAMPLE"/"$SAMPLE"_split.bam"

echo merging the two BAMs..
samtools merge $OUTBAM clipped_sorted.bam no_clipped_sorted.bam

echo indexing the final BAM output...
samtools index $OUTBAM