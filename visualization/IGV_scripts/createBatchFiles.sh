#!/bin/bash

slop=100

for i in *.xml;do
    base=$(basename $i .xml)
    mkdir ${base}_snaps;

    a=$base.igv;
    #echo "load $(pwd -P)/$i">$a;
    echo "new" >>$a;
    echo "maxPanelHeight 1250" >>$a;
    echo "genome hg19" >>$a;
    echo "load /Users/lsantuari/Documents/mount_points/hpc_giab/GIAB12878/mapping/GIAB12878_dedup.bam" >>$a;
    echo "load /Users/lsantuari/Documents/mount_points/hpc_ont_giab/NA12878_ONT_sorted.bam" >>$a;
    echo "load /Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/Mills2011_nanosv_full_inclusion.unique.bed" >>$a;
    echo "load /Users/lsantuari/Documents/Processed/NA12878/Overlap_diagrams/lumpy-Mills2011-call-set.bed" >>$a;
    echo "snapshotDirectory $(pwd -P)/${base}_snaps">>$a;
    #echo "maxPanelHeight 750" >>$a

    awk -v var=$slop '{print "goto "$1":"$2-var"-"$5+var"\nsort position\nviewaspairs\ncollapse\nsnapshot "$1"-"$2"-"$5"_"$7".png"}' $base.bedpe >>$a;
    echo "exit" >>$a
done
echo "Launch IGV, then use 'Tools > Run Batch Script' to generate images"
