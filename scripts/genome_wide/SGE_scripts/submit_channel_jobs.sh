#!/bin/bash -x

SAMPLE=$1
BAM=$2
BAM_SV=$3
TWOBIT=$4
MAP=$5
OUTPATH=$6

CHRARRAY=(`seq 1 22` 'X' 'Y')

if [ SAMPLE == 'CHM1_CHM13' ]
then
    for i in 0..${#CHRARRAY[@]}; do
        CHRARRAY[$i]='chr'CHRARRAY[$i]
    done
fi

OUTDIR=$OUTPATH"/"$SAMPLE

[ ! -d $OUTDIR ] && mkdir -p $OUTDIR

for PRG in clipped_read_pos clipped_reads split_reads; do

    JOB_NAME=$SAMPLE"_channels"

    qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,BAMARG=$BAM_SV,PRGARG=${PRG},OUTARG=$OUTDIR,CHRLIST=${CHRARRAY[@]} \
        -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
done

<<COMMENT
for CHROMOSOME in ${CHRARRAY[@]}; do

    for PRG in coverage; do

        JOB_NAME=$SAMPLE"_channels"

        qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

    for PRG in snv; do

        JOB_NAME=$SAMPLE"_channels"

        qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=$OUTDIR,TWOBIT=$TWOBIT \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

    for PRG in clipped_read_distance; do

        JOB_NAME=$SAMPLE"_channels"

        qsub -wd $OUTDIR -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM_SV,PRGARG=${PRG},OUTARG=$OUTDIR \
            -N $JOB_NAME -o $JOB_NAME".out" -e $JOB_NAME".err" make_channel.sge
    done

done
COMMENT


PRG='chr_array'

OUTPUTDIR=${OUTPATH}"/"$SAMPLE"/"$PRG
[ ! -d $OUTPUTDIR ] && mkdir -p $OUTPUTDIR

for CHROMOSOME in ${CHRARRAY[@]}; do
#for CHROMOSOME in 1; do

    OUTDIR=$OUTPATH
    JOB_NAME=$SAMPLE"_carray"
    JOB_NAME_HOLD=$SAMPLE"_channels"

    qsub -wd $OUTDIR -hold_jid $JOB_NAME_HOLD -v SAMPLEARG=$SAMPLE,CHRARG=$CHROMOSOME,BAMARG=$BAM,PRGARG=${PRG},OUTARG=${OUTDIR},TWOBIT=$TWOBIT,MAP=$MAP \
        -N $JOB_NAME -o $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".out" -e $OUTDIR"/"$SAMPLE"/"${PRG}"/"$JOB_NAME".err" make_channel.sge
done