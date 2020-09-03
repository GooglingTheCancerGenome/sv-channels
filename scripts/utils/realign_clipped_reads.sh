#!/bin/bash

BWA_INDEX=$(realpath -s "$1")
BAM=$(realpath -s "$2")
SAMPLE=$(basename "$BAM" .bam)
OUTDIR=$3
OUTBAM_FILENAME=$4
OUTBAM=${OUTDIR}/${OUTBAM_FILENAME}
TMPDIR=${OUTDIR}"/tmp"

# create tmp folder
[ ! -d $TMPDIR ] && mkdir -p $TMPDIR

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
bwa mem -L 0,0 ${BWA_INDEX} ${TMPDIR}/clipped.fq > ${TMPDIR}/clipped.sam

echo converting SAM to BAM...
samtools view -b ${TMPDIR}/clipped.sam -o ${TMPDIR}/clipped.bam

echo sorting BAM with split reads...
samtools sort -o ${TMPDIR}/clipped_sorted.bam ${TMPDIR}/clipped.bam

echo indexing BAM with split reads..
samtools index ${TMPDIR}/clipped_sorted.bam

echo merging the two BAMs..
samtools merge $OUTBAM ${TMPDIR}/clipped_sorted.bam ${TMPDIR}/no_clipped_sorted.bam

echo indexing the final BAM output...
samtools index $OUTBAM