From [1000Genomes_data_indexes](https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/1000_genomes_project/README.1000genomes.GRCh38DH.alignment):

###Reference Genome - GRCh38 with Alts

ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

This was sourced from  
ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_full_plus_hs38d1_analysis_set.fna

The HLA sequence in the file came from Heng Li's bwakit distribution  
http://sourceforge.net/projects/bio-bwa/files/bwakit/bwakit-0.7.12_x64-linux.tar.bz2/download

```commandline
faToTwoBit GRCh38_full_analysis_set_plus_decoy_hla.fa GRCh38_full_analysis_set_plus_decoy_hla.2bit
```

```commandline
python ../github/sv-channels/scripts/utils/Ns_to_bed.py \
    -b GRCh38_full_analysis_set_plus_decoy_hla.Ns.bed \
    -t GRCh38_full_analysis_set_plus_decoy_hla.2bit \
    -c chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY
```

[ENCFF001TDO.bed](https://www.encodeproject.org/files/ENCFF001TDO/) is based on the hg19 assembly

versions:
- ucsc-fatotwobit v377
- bedtools v2.30.0

```commandline
gunzip hg38-blacklist.v2.bed.gz
cut -f1-3 hg38-blacklist.v2.bed > hg38-blacklist.v2.3fields.bed
cat hg38-blacklist.v2.3fields.bed ../../../reference/GRCh38_full_analysis_set_plus_decoy_hla.Ns.bed | \
    sortBed -i - | mergeBed -i - > exclude.N-and-encodebl.GRCh38.bed
```

