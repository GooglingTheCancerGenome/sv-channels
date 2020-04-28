#!/bin/bah

for s in DEL INS INV DUP TRA;do
for c in split_reads gridss manta lumpy delly; do

INPUTDIR=$1/${s}/${c}/
cat ${INPUTDIR}train_model_with_fit/train_*_test_*_cv*/predictions/cnn_predictions.bedpe | sortBed -i - | awk '{if($1==$4){print $0"\t*\t*\t" $5-$2}else{print $0"\t*\t*\t"0}}' > ${INPUTDIR}cnn_predictions.bedpe

done
done
