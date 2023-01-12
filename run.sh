#!/usr/bin/env bash

set -e

MY_ENV=wf

# activate conda env
eval "$(conda shell.bash hook)"
conda activate "$MY_ENV"
conda list

# Extract signals
svchannels extract-signals data/test.fasta data/test.bam -o test
# Convert VCF files (Manta callset and truth set) to BEDPE format
Rscript svchannels/utils/R/vcf2bedpe.R -i data/test.vcf \
                                       -o data/test.bedpe
Rscript svchannels/utils/R/vcf2bedpe.R -i data/vcf/manta_out/manta.vcf \
                                       -o test/manta.bedpe
# Generate channels
svchannels generate-channels --reference data/test.fasta test channels test/manta.bedpe
# Label SVs
svchannels label -f data/test.fasta.fai -o labels channels/sv_positions.bedpe data/test.bedpe
# Train the model
svchannels train channels/channels.zarr.zip labels/labels.json.gz -m model.keras
# Score SVs
svchannels score channels model.keras data/vcf/manta_out/manta.vcf sv-channels.vcf
