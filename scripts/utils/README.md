To create the test data:

`cd utils`

download the hs37d5 reference genome from:

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz

and extract it:

`gunzip hs37d5.fa.gz`

 
```
conda env create -n test-env -f environment.yaml
conda activate test-env

sh create_test_data.sh hs37d5.fa test_data
```