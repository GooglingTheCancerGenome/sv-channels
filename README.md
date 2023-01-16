# sv-channels

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4584797.svg)](https://doi.org/10.5281/zenodo.4584797)
[![Build Status](https://travis-ci.org/GooglingTheCancerGenome/sv-channels.svg?branch=master)](https://travis-ci.org/GooglingTheCancerGenome/sv-channels)

*sv-channels* is a Deep Learning workflow for filtering structural variants (SVs) in short read alignment data using 
one-dimensional Convolutional Neural Networks (CNN).

*sv-channels* currently supports deletions (DEL) called with [Manta](https://github.com/Illumina/manta).

The workflow includes the following key steps:

**Transform read alignments into channels**

For each pair of SV breakpoints, a 2D Numpy array called *window-pair* is constructed. The shape of a window is
[window_size*2+padding, number_of_channels], where the genomic interval encompassing each window is centered on the 
breakpoint position with a context of [-window_size/2, +window_size/2]. From all the reads overlapping this genomic 
interval and from the relative segment subsequence of the reference sequence *number_of_channels* channels are 
constructed, where each channel encode a signal that can be used for SV calling. The list of channels can be 
found [here](doc/channels_list.md). The two windows are joined as *linked-windows* with a zero padding 2D array of 
shape [10, *number_of_channels*] in between to avoid artifacts related to the CNN kernel at the interface between the 
two windows. The window-pairs are labelled as *DEL* when the breakpoint positions overlap the DEL callset used as 
the ground truth and *noDEL* otherwise.

**Labelling**

Window-pairs are labelled as DEL (a true deletion) or noDEL (a false positive call)
based on the overlap of the DEL breakpoints of the window-pair with the truth set.

**Model training**

The labelled *window-pairs* are used to train a 1D CNN to classify Manta SVs as either DEL (true deletions)
or noDEL (false positives).

**Scoring Manta DELs**

The model is run on the window-pairs of a test sample. The SV qualities for the Manta DELs (QUAL) of the test
sample are substituted with the posterior probabilities obtained by the model.

## Dependencies

-   [Conda](https://conda.io/)
-   listed in [`environment.yaml`](/environment.yaml)

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

-   **input**:
    - read alignment in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format
    - reference genome used to map the reads in [FASTA](https://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml) format
-   **output**:
    - SV callset in [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format

### Run sv-channels on test data
Install sv-channels:
```commandline
python setup.py install
```
Process test set:
1. Extract signals:
```commandline
svchannels extract-signals data/test.fasta data/test.bam -o test
```
2. Convert VCF files (Manta callset and truth set) to BEDPE format:
```commandline
Rscript svchannels/utils/R/vcf2bedpe.R -i data/test.vcf \
                                       -o data/test.bedpe
Rscript svchannels/utils/R/vcf2bedpe.R -i data/vcf/manta_out/manta.vcf \
                                       -o test/manta.bedpe
```
3. Generate channels:
```commandline
svchannels generate-channels --reference data/test.fasta test channels test/manta.bedpe
```
4. Label SVs:
```commandline
svchannels label -f data/test.fasta.fai -o labels channels/sv_positions.bedpe data/test.bedpe
```
5. Train the model:
```commandline
svchannels train channels/channels.zarr.zip labels/labels.json.gz \
    -m model.keras
```
6. Score SVs. Note that _channels_ should be the channels folder generated from the hold-out test sample and _manta.vcf_ should be the Manta callset called on the hold-out test sample. For the Continuous Integration, we are using the same test data as in the training step for testing purpose.
```commandline
svchannels score channels model.keras data/vcf/manta_out/manta.vcf sv-channels.vcf
```

## Contributing

If you want to contribute to the development of _sv-channels_,
have a look at the [CONTRIBUTING.md](CONTRIBUTING.md).

## License

Copyright (c) 2023, Netherlands eScience Center

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
