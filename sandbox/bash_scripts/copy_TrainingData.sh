#!/bin/bash

DATE="041018"
SVMODE="TRA"
for SAMPLE in G1 N1; do
        DIR="../TrainingData_"$DATE"/"$SAMPLE"/ChannelData"
        [ ! -d $DIR ] && mkdir -p $DIR
        cp "../Training_"$SVMODE"/"$SAMPLE"/channel_maker_"*"_germline/"*"npy.gz" "../TrainingData_"$DATE"/"$SAMPLE"/ChannelData/"
        cp -r "../Training_"$SVMODE"/"$SAMPLE"/channel_maker_"*"_germline/label" "../TrainingData_"$DATE"/"$SAMPLE"/LabelData"
done