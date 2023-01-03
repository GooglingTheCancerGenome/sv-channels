#!/bin/bash
#SBATCH -t 96:00:00
#SBATCH --mem=256G


#debug="/usr/bin/time -v"
debug="echo"

PATH=$PATH:/hpc/compgen/users/bpedersen/conda/bin


echo ${CHR}

eval "$(conda shell.bash hook)"
sleep 2
conda activate sv-channels
#mamba activate sv-channels

RESDIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/results"
BAMDIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/BAM"
SVCHANNELSDIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/GitHub/sv-channels"

REFERENCEDIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/reference"
SVDATADIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/sv-channels-data"
TRUTHSETSDIR="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/truth-sets/1KG"

#mkdir -p ${RESDIR}/loco

set -euo pipefail

#echo 'running extract_signals'

#${debug} python ${SVCHANNELSDIR}/svchannels/extract_signals.py ${REFERENCEDIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa ${BAMDIR}/${SAMPLE}.bam -o ${SAMPLE}

#echo 'running generate_channels'

#${debug} python ${SVCHANNELSDIR}/svchannels/generate_channels.py ${SAMPLE}  ${SVDATADIR}/${SAMPLE}/manta.bedpe

#${debug} python ${SVCHANNELSDIR}/svchannels/label_windows.py \
#	--fai ${REFERENCEDIR}/GRCh38_full_analysis_set_plus_decoy_hla.fa.fai \
#	--ground_truth ${TRUTHSETSDIR}/${SAMPLE}.DEL.bedpe

sample_array=(HG03992 HG00420 HG01881 HG01114 HG02924 NA06991 HG02018)
validation=HG01053

#TEST=${SAMPLE}
#HG02924
#HG01114
#HG01881
#HG00420
#HG03992
#NA06991
SVTYPE=DEL

#echo ${TEST}
training_array=( ${sample_array[@]} )

outdir=${RESDIR}/models_loco/validation_${validation}_loco
mkdir -p ${outdir}

program=train_loco

windows=""
labels=""
TRAINED=""
samples=""
for s in ${training_array[@]};
do
        windows+=${RESDIR}/${s}/sv_chan.zarr.zip,
	labels+=${RESDIR}/${s}/labels.json.gz,
        TRAINED+=${s}_
	samples+=${s},
done

models_list=""
chr_list=""
for i in $(seq 1 21); 
do 
	models_list+=/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/results/models_loco/validation_HG01053_loco/chr${i}.manta_model.keras,
	chr_list+=chr${i},
done

# remove last comma
win_len=$(( ${#chr_list} - 1 ))
chr_list=${chr_list:0:$win_len}

win_len=$(( ${#models_list} - 1 ))
models_list=${models_list:0:$win_len}

win_len=$(( ${#windows} - 1 ))
windows=${windows:0:$win_len}

win_len=$(( ${#labels} - 1 ))
labels=${labels:0:$win_len}

win_len=$(( ${#TRAINED} - 1 ))
TRAINED=${TRAINED:0:$win_len}

win_len=$(( ${#samples} - 1 ))
samples=${samples:0:$win_len}

hparams="/hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/results/models/validation_HG01053_test_HG00420/manta_model.hparams.npy"

cmd="${debug} python ${SVCHANNELSDIR}/svchannels/${program}.py \
        --windows ${windows} \
	--labels ${labels} \
	-test ${CHR} \
	-s ${SVTYPE} \
        -m ${outdir}/manta_model.keras \
        -l ${outdir}/manta_model.log \
        -p ${hparams}"

echo 'running create_model'

${cmd}


OUTPUTDIR=${outdir}/results
NREGIONS=${REFERENCEDIR}/reference_N_regions.bed
TWOBIT=${REFERENCEDIR}/GRCh38_full_analysis_set_plus_decoy_hla.2bit
ENCODE=${REFERENCEDIR}/ENCFF001TDO.bed

debug="/usr/bin/time -v"
#debug=echo

p=predict_loco

cmd="${debug} python ${SVCHANNELSDIR}/svchannels/${p}.py \
        -m ${models_list} \
	--labels ${labels} \
        -n ${chr_list} \
        -sm ${samples} \
        --windows ${windows} \
        -s $SVTYPE \
        -o $OUTPUTDIR \
	--sv_channels $SVCHANNELSDIR \
        -fn $NREGIONS \
	--encode_blacklist $ENCODE \
        -tb $TWOBIT"

echo 'running predict_loco'
${cmd}

#$debug Rscript ${SVCHANNELSDIR}/svchannels/plot_evaluation.R \
#	-f ${OUTPUTDIR}/sv-channels.DEL \
#	--manta /hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/Manta/${SAMPLE}/manta.vcf \
#	--gridss /hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/GRIDSS2/${SAMPLE}.vcf \
#	-t /hpc/compgen/projects/googling-the-cancer-genome/sv-channels/analysis/SPIDER/1KG/truth-sets/1KG/${SAMPLE}.DEL.vcf \
#	-n ${SAMPLE} \
#	--svchannels ${OUTPUTDIR}/sv-channels.DEL.vcf

