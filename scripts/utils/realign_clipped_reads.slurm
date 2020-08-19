#!/bin/bash
#SBATCH --time=48:00:00
#SBATCH --mem=32G

BAM=$(realpath -s "$1")
SAMPLE=$(basename "$BAM" .bam)
OUTDIR=$2
TMPDIR=$OUTDIR"/"$SAMPLE"/tmp"

# create tmp folder
[ ! -d $TMPDIR ] && mkdir -p $TMPDIR

# echo activating...
# conda activate utils

echo clipped to split reads...
python split_clipped_reads.py \
    -b $BAM \
    -ob ${TMPDIR}/no_clipped.bam \
    -of ${TMPDIR}/clipped.fq \
    -l ${TMPDIR}/split_clipped_reads.log

echo sorting BAM without clipped reads...
samtools sort -o ${TMPDIR}/no_clipped_sorted.bam ${TMPDIR}/no_clipped.bam

echo indexing BAM without clipped reads...
samtools index ${TMPDIR}/no_clipped_sorted.bam

echo mapping soft clipped sequences...
bwa mem -L 0,0 ../../data/test.fasta ${TMPDIR}/clipped.fq > ${TMPDIR}/clipped.sam

echo converting SAM to BAM...
samtools view -b ${TMPDIR}/clipped.sam -o ${TMPDIR}/clipped.bam

echo sorting BAM with split reads...
samtools sort -o ${TMPDIR}/clipped_sorted.bam ${TMPDIR}/clipped.bam

echo indexing BAM with split reads..
samtools index ${TMPDIR}/clipped_sorted.bam

OUTBAM=$OUTDIR"/"$SAMPLE"/"$SAMPLE"_split.bam"

echo merging the two BAMs..
samtools merge $OUTBAM ${TMPDIR}/clipped_sorted.bam ${TMPDIR}/no_clipped_sorted.bam

echo indexing the final BAM output...
samtools index $OUTBAM