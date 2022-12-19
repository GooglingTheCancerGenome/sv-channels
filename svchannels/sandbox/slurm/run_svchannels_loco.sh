#!/bin/bash

#set -euo pipefail
#for i in `seq 2 21`;
for i in 1;
do 
	CHR=chr${i}
	sbatch -J ${CHR}.loco.eval.svchan -p gpu --gpus-per-node=RTX6000:1  -e ${CHR}.loco.eval.err -o ${CHR}.loco.eval.out --export=CHR=${CHR} run_evaluation_loco.slurm
done
