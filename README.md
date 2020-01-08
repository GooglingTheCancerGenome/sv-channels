[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/CNN.svg?branch=dev)](https://travis-ci.org/GooglingTheCancerGenome/CNN)

# ChannelMaker

**Install & execute**

```bash
conda update -y conda  # update Conda
conda env create -n cm -f environment.yaml
conda activate cm

./run_local.sh data/test/chr22.bam chr22       # run locally
./run.sh gridengine data/test/chr22.bam chr22  # submit jobs to GE cluster
./run.sh slurm data/test/chr22.bam chr22       # submit jobs to Slurm cluster
```
