#!/usr/bin/env bash

# 1st argument: sv-callers output folder
# 2nd argument test data folder of CNN

cp $1/all.vcf $2/vcf/all.vcf
for s in gridss manta delly lumpy; do
	cp $1/$s_out/$s.vcf $2/vcf/$s_out/$s.vcf
done
