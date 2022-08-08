rm(list = ls())
suppressPackageStartupMessages(require(argparser))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))


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

load_bed <- function(bed_file)
{
  bed_regions <- rtracklayer::import(bed_file)
  # set NCBI seqlevels
  seqlevelsStyle(bed_regions) <- "UCSC"
  return(bed_regions)
}

load_vcf <- function(vcf_file, svtype, caller, filter.regions)
{
  # Load VCF file
  vcf_gr <-
    VariantAnnotation::readVcf(vcf_file)
  
  # Keep only SVs that passed the filtering (PASS or .)
  vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("PASS", ".")]
  
  if (caller == 'lumpy')
  {
    # Read evidence support as a proxy for QUAL
    fixed(vcf_gr)$QUAL <- unlist(info(vcf_gr)$SU)
    
  } else if (caller == 'delly')
  {
    # Split-read support plus Paired-end read support as a proxy for QUAL
    sr_support <- info(vcf_gr)$SR
    sr_support[is.na(vcf_gr)] <- 0
    fixed(vcf_gr)$QUAL <- sr_support + info(vcf_gr)$PE
  }
  
  vcf_gr <- StructuralVariantAnnotation::breakpointRanges(vcf_gr)
  vcf_gr <- infer_svtype(vcf_gr)
  
  # Select only one SV type
  vcf_gr <- vcf_gr[which(vcf_gr$svtype == svtype)]
  
  message(paste('Loading', length(vcf_gr), svtype, 'calls for', caller, sep=" "))
  
  # Select SVs >= 50 bp
  if (! svtype %in% c('CTX', 'INS'))
  {
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]
  }
  
  # Filter regions
  vcf_gr <- filter_regions('ENCODE blacklist', vcf_gr, load_bed(filter.regions), mode = 'remove')
  
  # select only chromosomes 1 to 22
  vcf_gr <-
    vcf_gr[seqnames(vcf_gr) %in% paste('chr', c(1:22), sep = '')]
  vcf_gr <- vcf_gr[vcf_gr$partner %in% names(vcf_gr)]
  
  return(vcf_gr)
}

filter_regions <-
  function(filter.name, regions_to_filter, ref_regions, mode = 'remove')
  {
    print(length(regions_to_filter))
    if (mode == 'keep')
    {
      result <-
        regions_to_filter[overlapsAny(regions_to_filter, ref_regions) &
                            overlapsAny(partner(regions_to_filter), ref_regions),]
    } else if (mode == 'remove') {
      result <- regions_to_filter[! (
        overlapsAny(regions_to_filter, ref_regions) |
          overlapsAny(partner(regions_to_filter), ref_regions)
      ),]
    }
    message(paste(length(result), 'calls after filtering for', filter.name, sep=" "))
    return(result)
  }

svtype <- "DEL"
bed_filter <- '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/reference/exclude.N-and-encodebl.bed'

sv_regions <- list()
for(default_sample in c("HG00736", "HG00737", "HG00738"))
{
  truth_set <- paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/truth_sets_DEL/',default_sample,'.DEL.vcf',sep='')
  gridss <- paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/GRIDSS2/',default_sample,'.vcf.gz',sep='')
  manta <- paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/Manta/',default_sample,'/manta.vcf',sep='')
  svchannels <- paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results_locso/simplearc-hyperopt/trio-tnsv/sv-channels.',default_sample,'.DEL.vcf',sep='')
  
  sv_regions[paste(default_sample,'truth_set',sep='_')] <-
    load_vcf(truth_set, svtype, "truth_set", bed_filter)
  sv_regions[paste(default_sample,'gridss',sep='_')] <-
    load_vcf(gridss, svtype, "gridss", bed_filter)
  sv_regions[paste(default_sample,'manta',sep='_')] <-
    load_vcf(manta, svtype, "manta", bed_filter)
  sv_regions[paste(default_sample,'svchannels',sep='_')] <-
    load_vcf(svchannels, svtype, "svchannels", bed_filter)
  
}

names(sv_regions)
require(UpSetR)

hits <- list()

for(i in 1:(length(sv_regions)))
{
  if(names(sv_regions)[i]!='HG00738_truth_set')
  {
    print(names(sv_regions)[i])
    my.hits <- findBreakpointOverlaps(
      sv_regions[[i]],
      sv_regions[['HG00738_truth_set']],
      # read pair based callers make imprecise calls.
      # A margin around the call position is required when matching with the truth set
      maxgap = 100,
      # Since we added a maxgap, we also need to restrict the mismatch between the
      # size of the events. We don't want to match a 100bp deletion with a
      # 5bp duplication. This will happen if we have a 100bp margin but don't also
      # require an approximate size match as well
      sizemargin = 0.25,
      ignore.strand = FALSE,
      # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
      # by restricting the margin based on the size of the event, we can make sure
      # that simple events actually do overlap
      restrictMarginToSizeMultiple = 0.5
    )
    hits[[names(sv_regions)[i]]] <- subjectHits(my.hits)
  }
}

hits[['HG00738_gridss']]

upset(fromList(hits), nsets = length(hits), order.by = "freq")
