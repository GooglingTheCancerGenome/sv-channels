rm(list = ls())

suppressPackageStartupMessages(require(argparser))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))

# The palette with black:
cbbPalette <-
  c(
    "#000000",
    "#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7"
  )

#default_sample <- 'HG01053'
default_sample <- 'HG00420'

make_percent <- function(x) {
  signif(x * 100, digits = 4)
}

infer_svtype <- function(gr)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "TRA",
      # Using TRA instead of ITX or BP
      ifelse(
        gr$insLen >= abs(gr$svLen) * 0.7,
        "INS",
        ifelse(
          strand(gr) == strand(partner(gr)),
          "INV",
          ifelse(xor(
            start(gr) < start(partner(gr)), strand(gr) == "-"
          ), "DEL",
          "DUP")
        )
      )
    )
  return(gr)
}

load_vcf <- function(vcf_file, svtype, caller)
{
  
  message(paste('Loading ', svtype, 'calls for', caller))
  
  # Load VCF file
  vcf_gr <-
    VariantAnnotation::readVcf(vcf_file)
  
  # Keep only SVs that passed the filtering (PASS or .)
  vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("PASS", ".")]
  
  #if(caller=='truth_set')
  #{
  #vcf_gr <- vcf_gr[which(unlist(info(vcf_gr)$ALGORITHMS)!='manta')]
  #}
  
  if (caller == 'lumpy')
  {
    # Read evidence support as a proxy for QUAL
    support <- unlist(info(vcf_gr)$SU)
    fixed(vcf_gr)$QUAL <- support
  } else if (caller == 'delly')
  {
    # Split-read support plus Paired-end read support as a proxy for QUAL
    sr_support <- info(vcf_gr)$SR
    sr_support[is.na(vcf_gr)] <- 0
    fixed(vcf_gr)$QUAL <-
      sr_support + info(vcf_gr)$PE
  }
  
  vcf_gr <- breakpointRanges(vcf_gr)
  vcf_gr <- infer_svtype(vcf_gr)
  
  # Select only one SV type
  print(table(vcf_gr$svtype))
  vcf_gr <- vcf_gr[which(vcf_gr$svtype == svtype)]
  message(length(vcf_gr))
  # Select SVs >= 50 bp
  if (!svtype %in% c('TRA', 'INS'))
  {
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]
  }
  # select only chromosomes 1 to 22
  vcf_gr <-
    vcf_gr[seqnames(vcf_gr) %in% paste('chr', c(1:22), sep = '')]
  vcf_gr <- vcf_gr[vcf_gr$partner %in% names(vcf_gr)]
  
  return(vcf_gr)
}