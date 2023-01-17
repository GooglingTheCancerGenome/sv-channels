# sv-channels

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4584797.svg)](https://doi.org/10.5281/zenodo.4584797)
[![CI](https://github.com/GooglingTheCancerGenome/sv-channels/actions/workflows/ci.yaml/badge.svg)](https://github.com/GooglingTheCancerGenome/sv-channels/actions/workflows/ci.yaml)

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
# install Mamba
conda install -n base -c conda-forge -y mamba
# create a new environment with dependencies & activate it
mamba env create -n sv-channels -f environment.yaml
conda activate sv-channels
# install svchannels CLI
python setup.py install
```

**3. Execution.**

-   **input**:
    - read alignment in [BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) format
    - reference genome used to map the reads in [FASTA](https://www.ncbi.nlm.nih.gov/BLAST/fasta.shtml) format
-   **output**:
    - SV callset in [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) format

_1. Convert VCF files (Manta callset and truth set) to BEDPE format._
```
Rscript svchannels/utils/R/vcf2bedpe.R -i data/test.vcf \
                                       -o data/test.bedpe
Rscript svchannels/utils/R/vcf2bedpe.R -i data/vcf/manta_out/manta.vcf \
                                       -o test/manta.bedpe
```

_2. Extract signals._
```
svchannels extract-signals data/test.fasta data/test.bam -o test
```

_3. Generate channels._
```
svchannels generate-channels --reference data/test.fasta test channels test/manta.bedpe
```

_4. Label SVs._
```
svchannels label -f data/test.fasta.fai -o labels channels/sv_positions.bedpe data/test.bedpe
```

_5. Train the model._
```
svchannels train channels/channels.zarr.zip labels/labels.json.gz -m model.keras
```

_6. Score SVs._
```
svchannels score channels model.keras data/vcf/manta_out/manta.vcf sv-channels.vcf
```

_7. Evaluate results._
```
faToTwoBit data/test.fasta data/test.2bit
grep 'DEL' data/test.bedpe > data/test.DEL.bedpe
python svchannels/utils/python/bedpe_to_vcf.py -i data/test.DEL.bedpe -o data/test.proper.vcf -b data/test.2bit
Rscript svchannels/utils/R/plot_evaluation_test_data.R \
    --svchannels sv-channels.vcf \
    --manta data/vcf/manta_out/manta.vcf \
    --gridss data/vcf/gridss_out/gridss.vcf \
    -t data/test.proper.vcf \
    -f sv-channels.metric
```

| caller      | calls | TP | FP   | precision (%) | recall (%) | F1 score (%) |
|-------------|-------|----|------|---------------|------------|-----------|
| gridss      | 1422 | 1174 | 248 | 82.56         | 59.78      | 69.35     |
| manta       | 1530 | 1494 | 36  | 97.65         | 76.07      | 85.52     |
| sv-channels | 1528 | 1528 | 34  | 97.77         | 76.07      |  85.57    |

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
