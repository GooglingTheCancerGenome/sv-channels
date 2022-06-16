#!/bin/bash

#set -euo pipefail
SAMPLE=HG01053
#for i in `seq 1 22`;
for i in X;
do 
	CHR=chr${i}
	sbatch -J ${SAMPLE}.${CHR}.loco.eval.svchan -p cpu  -e ${SAMPLE}.${CHR}.loco.eval.err -o ${SAMPLE}.${CHR}.loco.eval.out --export=CHR=${CHR},SAMPLE=${SAMPLE} run_evaluation_locso.slurm
done
