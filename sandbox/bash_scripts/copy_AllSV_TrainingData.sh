#!/bin/bash

DATE="091118"
SVMODE="INDEL"
#for SVMODE in INDEL INV DUP TRA; do
echo "copying "$SVMODE
for SAMPLE in G1 N1; do
	echo "sample "$SAMPLE
	DIR="../TrainingData_"$DATE"/"$SVMODE"/"$SAMPLE"/ChannelData"
	[ ! -d $DIR ] && mkdir -p $DIR
	cp "../Training_"$SVMODE"/"$SAMPLE"/channel_maker_"*"_germline/"*"npy.gz" $DIR"/"$SAMPLE".npy.gz"
	cp -r "../Training_"$SVMODE"/"$SAMPLE"/channel_maker_"*"_germline/label" $DIR"/../LabelData"
done
done
