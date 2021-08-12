suppressPackageStartupMessages(require(UpSetR))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))

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