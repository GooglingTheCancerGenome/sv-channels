#!/usr/bin/env bash

set -xe

# If directory for script used to store per chromosome output files does not exist, create it
[ ! -d ${OUTARG}"/"${PRGARG} ] && mkdir -p ${OUTARG}"/"${PRGARG}

# which sample is considered
echo ${SAMPLEARG}
# which script is being executed
echo ${PRGARG}
# which BAM file is considered
echo ${BAMARG}
# which chromosome is considered
echo ${CHRARG}
# which window is considered
echo ${WINDOWARG}

echo "Running ${PRGARG} for ${BAMARG} ${CHRARG}"

# coverage.py will output file with extension '.npy.bz2'
if [ ${PRGARG} == 'coverage' ]; then
        python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
                --out ${PRGARG}".npy" \
                --outputpath ${OUTARG} \
                --logfile ${PRGARG}".log"

# snv.py will output file with extension '.npz.gz'
elif [ ${PRGARG} == 'snv' ]; then
        python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
                --out ${PRGARG}".npy" \
                --outputpath ${OUTARG} \
                --logfile ${PRGARG}".log"

# ChannelMaker for real data. Flag 'train' is only used to output the NoSV category for training, that for now it is
# not implemented in the channel_make_train.py script
elif [ ${PRGARG} == 'channel_maker_real_germline' ]; then
        python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
		        --sample ${SAMPLEARG} \
		        --train True \
		        --svmode ${SVMODEARG} \
                --out ${SAMPLEARG}"_"${CHRARG}".npy" \
                --outputpath ${OUTARG} \
                --window ${WINDOWARG} \
                --logfile ${SAMPLEARG}"_"${CHRARG}".log"

# ChannelMaker for training data
elif [ ${PRGARG} == 'channel_maker_train_germline' ]; then
        python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
                --sample ${SAMPLEARG} \
		        --train True \
		        --svmode ${SVMODEARG} \
                --out ${OUTARG}"/"${PRGARG}"/"$SAMPLEARG".npy.gz" \
                --logfile ${OUTARG}"/"${PRGARG}"/"${CHRARG}"_${PRGARG}.log"

elif [ ${PRGARG} == 'chr_array' ]; then
        python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
                --sample ${SAMPLEARG} \
                --out ${SAMPLEARG}"_"${CHRARG}".hdf5" \
                --outputpath ${OUTARG} \
                --logfile ${SAMPLEARG}"_"${CHRARG}".log"

elif [ ${PRGARG} == 'clipped_read_distance' ]; then
# generic execution
	python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --chr ${CHRARG} \
                --out ${PRGARG}".json.gz" \
                --outputpath ${OUTARG} \
                --logfile ${PRGARG}".log"

else
# generic execution
	python "${PRGARG}.py" \
                --bam ${BAMARG} \
                --out ${PRGARG}".json.gz" \
                --outputpath ${OUTARG} \
                --logfile ${PRGARG}".log"
fi