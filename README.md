# sv-channels

[![DOI](https://zenodo.org/badge/DOI/10.000/FIXME.svg)](https://doi.org/10.000/FIXME)
[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-channels.svg?branch=master)](https://travis-ci.org/GooglingTheCancerGenome/sv-channels)

**sv-channels** is a Deep Learning workflow for calling structural variants (SVs) in short read alignment data (BAM format) using 1-dimensional Convolutional Neural Networks. It was tested on a benchmark dataset with [three cell lines](https://zenodo.org/record/4001614#.YCz4ixNKhlM): two samples from the [Genome in a Bottle](https://www.nist.gov/programs-projects/genome-bottle) consortium (NA12878 and NA24385) and one synthetic diploid sample (CHM1_CHM13).

**sv-channels** supports five major types of SVs:
- deletions (DEL);
- insertions (INS);
- inversions (INV);
- tandem duplications (DUP);
- inter-chromosomal translocations (CTX).

# Dependencies

-   [Python 3](https://www.python.org/)
-   [Conda](https://conda.io/) - package/environment management system
-   [Xenon CLI](https://github.com/NLeSC/xenon-cli) - command-line interface to compute and storage resources

The [workflow](doc/sv-channels.svg) includes the following (major) steps:

Model training:
- A
- B

**1. Clone this repo.**

```bash
git clone https://github.com/GooglingTheCancerGenome/sv-channels.git
cd sv-channels
```

**2. Install dependencies.**

```bash
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# create & activate new env with installed deps
conda env create -n sv-channels -f environment.yaml
conda activate sv-channels
```

**3. Execution.**

-   **input files**:
    - blabla
-   **output files**:
    - blabla

```bash
SCH=local  # gridengine or slurm
BAM=data/test.bam
SEQIDS="12,22"
SVTYPES="INV,DEL,INS,INV,DUP,CTX"
./run.sh $SCH $BAM $SEQIDS $SVTYPES # run jobs locally or on a compute cluster
```

# Contributing

If you want to contribute to the development of _sv-channels_,
have a look at the [CONTRIBUTING.md](CONTRIBUTING.md).

# License

Copyright (c) 2020, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
