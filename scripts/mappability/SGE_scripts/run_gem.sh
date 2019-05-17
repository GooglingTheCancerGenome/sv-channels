#!/bin/bash -x
#$ -S /bin/bash
#$ -N GEM
#$ -l h_rt=24:00:00
#$ -l h_vmem=24G
#$ -q all.q
#$ -e gem.err
#$ -o gem.out
#$ -wd /hpc/cog_bioinf/ridder/users/lsantuari/Datasets/Mappability
#$ -pe threaded 8

#gem-indexer -i /hpc/cog_bioinf/GENOMES/Homo_sapiens.GRCh37.GATK.illumina/Homo_sapiens.GRCh37.GATK.illumina.fasta -o GRCh37

#READ_LENGTH=151
READ_LENGTH=250

gem-mappability -I GRCh37.gem -l $READ_LENGTH -o GRCh37.$READ_LENGTH\mer -T 8

gem-2-wig -I GRCh37.gem -i GRCh37.$READ_LENGTH\mer.mappability -o GRCh37.$READ_LENGTH\mer

wigToBigWig GRCh37.$READ_LENGTH\mer.wig GRCh37.$READ_LENGTH\mer.sizes GRCh37.$READ_LENGTH\mer.bw

bigWigToBedGraph GRCh37.$READ_LENGTH\mer.bw  GRCh37.$READ_LENGTH\mer.bedGraph