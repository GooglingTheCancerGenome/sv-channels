#!/bin/bash

DATE=""
for SAMPLE in G1 N1; do
	DIR="../TrainingData_"$DATE"/"$SAMPLE"/ChannelData"
	[ ! -d $DIR ] && mkdir -p $DIR
	cp "../Training/"$SAMPLE"/channel_maker_*_germline/*npy.gz" "../TrainingData_"$DATE"/"$SAMPLE"/ChannelData/"
	cp "../Training/"$SAMPLE"/channel_maker_*_germline/label" "../TrainingData_"$DATE"/"$SAMPLE"/LabelData"
done