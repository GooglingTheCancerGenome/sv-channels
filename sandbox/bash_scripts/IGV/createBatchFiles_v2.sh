#!/bin/bash

# Requires a XML session file saved from IGV
# Takes in input a BED file with positions

slop=50

for i in *.xml;do
    base=$(basename $i .xml)
    mkdir ${base}_snaps;

    a=$base.igv;
    #echo "load $(pwd -P)/$i">$a;
    echo "new" >>$a;
    echo "maxPanelHeight 1250" >>$a
    echo "genome hg19" >>$a;
    #echo "load /hpc/cog_bioinf/diagnostiek/projects/na12878_wgs_trio/GIAB12878/mapping/GIAB12878_dedup.bam" >>$a;
    echo "load /hpc/cog_bioinf/ridder/users/lsantuari/Datasets/CretuStancu2017/Patient1/Patient1.bam" >>$a;

    echo "snapshotDirectory $(pwd -P)/${base}_snaps">>$a;
    #echo "maxPanelHeight 750" >>$a

    awk '{print "goto "$1":"$2-$slop"-"$3+$slop"\nsort position\nviewaspairs\ncollapse\nsnapshot "$1"-"$2"-"$3"_"$4".png"}' $base.bed >>$a;
    echo "exit" >>$a
done
echo "Launch IGV, then use 'Tools > Run Batch Script' to generate images"
