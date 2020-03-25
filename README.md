[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/CNN.svg?branch=iss32)](https://travis-ci.org/GooglingTheCancerGenome/CNN)

# ChannelMaker

**Install & execute**

```bash
conda update -y conda  # update Conda
conda env create -n cm -f environment.yaml
conda activate cm

SCH=gridengine  # or slurm
BAM=data/test/chr12_chr22.44Mb_45Mb.GRCh37.bam
CHRS="12 22"
./run_local.sh $BAM $CHRS  # run locally
./run.sh $SCH $BAM         # submit jobs to GridEngine/Slurm cluster
```
