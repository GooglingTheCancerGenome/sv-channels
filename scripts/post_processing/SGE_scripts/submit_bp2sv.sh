#!/bin/bash -x
#$ -S /bin/bash
#$ -N BP2SV
#$ -l h_rt=24:00:00
#$ -l h_vmem=64G
#$ -q all.q
#$ -e BP2SV.err
#$ -o BP2SV.out
#$ -wd /hpc/cog_bioinf/ridder/users/lsantuari/Processed/Benchmark_GPU/Git/CNN/scripts
#$ -M L.Santuari@umcutrecht.nl
#$ -m beas

echo "Load python"
source activate channel-maker
python "post_processing/breakpoint2sv.py"