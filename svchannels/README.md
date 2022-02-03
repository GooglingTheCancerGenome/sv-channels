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

```bash
SCH=local  # gridengine or slurm
BAM=../data/test.bam
./run.sh $SCH $BAM # run jobs locally or on a compute cluster
```
