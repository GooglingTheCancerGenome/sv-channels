# Test data

**1. Install dependencies.**

```bash
# create & activate new env
conda env create -n utils -f environment.yaml
conda activate utils
```
**2. Download the hs37d5 reference genome.**

```bash
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip -d hs37d5.fa.gz
```

**3. Generate mappability track.**

```bash
./create_test_data.sh hs37d5.fa 100 test_data
```

**4. Generate artificial short-read alignments with SVs.**

https://github.com/GooglingTheCancerGenome/sv-gen
