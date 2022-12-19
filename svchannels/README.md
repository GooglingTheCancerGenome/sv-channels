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
### Run sv-channels on test data
Install sv-channels:
```commandline
python setup.py install
```
Process test set:
1. Extract signals:
```commandline
svchannels extract-signals data/test.fasta data/test.bam test
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
svchannels label -v channels/sv_positions.bedpe -g data/test.bedpe -f data/test.fasta.fai -p labels
```
5. Train the model:
```commandline
svchannels train channels/channels.zarr.zip labels/labels.json.gz channels/channels.zarr.zip labels/labels.json.gz \
    -m model.keras
```
6. Score SVs:
```commandline
svchannels score channels model.keras data/vcf/manta_out/manta.vcf sv-channels.vcf
```

