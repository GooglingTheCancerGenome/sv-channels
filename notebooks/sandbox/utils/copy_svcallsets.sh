#!/bin/bash
cp $1/all.vcf $2/all.vcf
for svcaller in gridss manta lumpy delly; do
	cp $1/$svcaller\_out/$svcaller.vcf $2/$svcaller\_out/$svcaller.vcf
done
