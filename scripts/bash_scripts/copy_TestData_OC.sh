#!/bin/bash

DATE="070119"
BASEDIR="../OC_COMPARISONS/"

for SAMPLE in `ls $BASEDIR`; do
	echo "Saving $SAMPLE..."
	DIR="../TestData_"$DATE"/"$SAMPLE"/ChannelData"
	[ ! -d $DIR ] && mkdir -p $DIR
	cp "$BASEDIR/"$SAMPLE/channel_maker_real_germline/*npy.gz "../TestData_"$DATE"/"$SAMPLE"/ChannelData/"
	LABELDIR="../TestData_"$DATE"/"$SAMPLE"/MultiLabelData"
        [ ! -d $LABELDIR ] && mkdir -p $LABELDIR
	cp "$BASEDIR/"$SAMPLE/label/*.gz "../TestData_"$DATE"/"$SAMPLE"/MultiLabelData/"
done
