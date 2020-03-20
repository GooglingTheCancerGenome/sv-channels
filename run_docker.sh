#!/usr/bin/env bash

set -e

# check input arg(s)
if [ $# -lt "1" ]; then
  echo "Usage: $0 [scheduler]"
  exit 1
fi

HOME=/home/xenon
PREFIX=gtcg/xenon
TAG=dev
BAM=data/test/chr12_chr22.44Mb_45Mb.GRCh37.bam
CHROM="12 22"
PORT=10023
SCH=$1  # slurm or gridengine

docker run -d -p $PORT:22 --name $SCH $PREFIX-$SCH:$TAG
sleep 10
docker ps -a
docker cp ./ $SCH:$HOME
docker exec -t $SCH chown -R xenon.xenon $HOME
docker exec -u xenon -t $SCH ./install.sh
docker exec -u xenon -t $SCH ./run.sh $SCH $BAM $CHROM
