From [1000Genomes_data_indexes](https://github.com/igsr/1000Genomes_data_indexes/blob/master/data_collections/1000_genomes_project/README.1000genomes.GRCh38DH.alignment):

### Reference Genome - GRCh38 with Alts

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

### Generation of test data
Extract the following regions from [hs37d5](ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz)
using [seqkit subseq](https://bioinf.shenwei.me/seqkit/). The file _chr12_chr22_40Mb_to_44Mb.bed_ contains:

```text
12   44000000        46000000
22   44000000        46000000
```

```commandline
seqkit subseq \
    --bed chr12_chr22_44Mb_to_46Mb.bed \
    hs37d5.fa.gz \
    > chr12_chr22_44Mb_to_46Mb.fasta
```
Add Ns to the sequences randomly:
```commandline
cp chr12_chr22_44Mb_to_46Mb.fasta test.fasta
for p in $(seq 2000000 | shuf -n 10); do
  seqkit mutate -p "$p:N" test.fasta -o test-$p.fasta
  mv test-$p.fasta test.fasta
done
sed -e 's/12_44000001-46000000:./12/' -e 's/22_44000001-46000000:./22/' test.fasta > test.fasta.renamed
mv test.fasta.renamed test.fasta
```
To save the regions containing Ns in BED format:
```commandline
seqkit locate -i -P -r -p "N+" --bed test.fasta -o test.Ns.bed
```
[sv-gen](https://github.com/GooglingTheCancerGenome/sv-gen) version 1.0.0 was used to generate the following data:
- htz-sv.bam
- htz-sv.bam.bai
- htz-sv.vcf
- hmz-sv.bam
- hmz-sv.bam.bai
- hmz-sv.vcf

with the following file _analysis.yaml_:
```yaml
---
threads: -1  # Samtools & BWA (default: -1 = use available cores)

# I/O files
input:
  fasta: data/test.fasta  # filepath of ref. genome (haploid)
  seqids: [12,22]  # zero or more SeqIDs (e.g. chromosomes)

output:
  basedir: data/out  # relative or absolute path
  genotype:   # diploid genomes
    - hmz     # homozygous
    - hmz-sv  # homozygous with SVs
    - htz-sv  # heterozygous with SVs

# registered I/O file extensions
filext:
  fasta: .fasta
  fasta_idx:
    - .fasta.ann  # BWA v0.6.x index files
    - .fasta.amb  #
    - .fasta.bwt  #
    - .fasta.pac  #
    - .fasta.sa   #
  fastq: .fq
  bam: .bam
  bam_idx: .bam.bai
  bed: .bed
  vcf: .vcf

simulation:
  # SURVIVOR parameters
  config: survivor.cfg
  svtype:
    dup: [1000, 50, 500]    # duplication: [count, min_len, max_len]
    inv: [1000, 50, 500]      # inversion: ""
    tra: [1000, 50, 500]    # translocation: ""
    indel: [2000, 50, 500]    # insertion+deletion: ""
    invdel: [0, 600, 800]   # inversion+deletion: ""
    invdup: [0, 600, 800]   # inversion+duplication: ""
  # ART parameters
  seed: 1000
  profile: HSXt
  coverage: [30]  # [cov1, cov2, ...] MEDIAN_COVERAGE=37 (picard CollectWgsMetrics)
  read:
    length: [150]  # [len1, len2, ...]
  insert:
    stdev: 100      # standard deviation of the fragment length (bp)
    length: [500]  # [len1, len2, ...]
```

