The model *model.keras* was trained on the following high-coverage 1KG samples:

- HG01053
- HG01114
- HG01881
- HG02018
- HG02924
- HG03992
- NA06991

with Manta (v1.1.0) DEL calls labelled using the following truth set:

[1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20210124.SV_Illumina_Integration/1KGP_3202.Illumina_ensemble_callset.freeze_V1.vcf.gz)

Training was performed with the following snakemake workflow:

[Snakemake](https://github.com/GooglingTheCancerGenome/sv-channels/blob/master/svchannels/cross-validations/workflow_0/Snakefile) 
