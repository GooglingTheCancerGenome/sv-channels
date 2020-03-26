#!/bin/bash

DATE="270519"
#for SAMPLE in NA12878 NA12891 NA12892 PATIENT1 PATIENT2; do

        DIR="../TestData_"$DATE"/"$SAMPLE"/ChannelData"
        [ ! -d $DIR ] && mkdir -p $DIR
        cp "../"$SAMPLE/channel_maker_real_germline/*npy.gz "../TestData_"$DATE"/"$SAMPLE"/ChannelData/"
        DIR="../TestData_"$DATE"/"$SAMPLE"/MultiLabelData"
        [ ! -d $DIR ] && mkdir -p $DIR
        cp "../"$SAMPLE/label_npy/*.gz "../TestData_"$DATE"/"$SAMPLE"/MultiLabelData/"
done