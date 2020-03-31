suppressPackageStartupMessages(require(StructuralVariantAnnotation))

#SV type assignment based on
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
apply_svtype <- function(gr)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "BP",
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
  gr
}

# Load the HG002_SVs_Tier1_v0.6.bed confidence regions
load_confidence_regions <- function(sample)
{
  if (sample == 'NA24385')
  {
    bed.file <-
      "/Users/lsantuari/Documents/Data/germline/NA24385/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed"
  } else{
    bed.file <-
      "/Users/lsantuari/Documents/Data/ENCODE/ENCFF001TDO.bed"
  }
  confidence_regions_gr <- rtracklayer::import(bed.file)
  seqlevelsStyle(confidence_regions_gr) <- "NCBI"
  confidence_regions_gr
}

# keep only SVs with both breakpoints overlapping confidence regions
remove_blacklist <- function(gr, confidence_regions_gr, sample)
{
  if (sample %in% c('NA24385'))
  {
    gr[overlapsAny(gr, confidence_regions_gr) &
         overlapsAny(partner(gr), confidence_regions_gr), ]
  } else{
    gr[!(
      overlapsAny(gr, confidence_regions_gr) |
        overlapsAny(partner(gr), confidence_regions_gr)
    ), ]
  }
  
}


load_sv_caller_vcf <-
  function(vcf_file,
           confidence_regions_gr,
           sample_name,
           sv_caller)
  {
    # vcf_file <- truth_set_file[[sample]]
    sv_callset_vcf <-
      VariantAnnotation::readVcf(vcf_file)
    
    sv_callset_vcf <- sv_callset_vcf[rowRanges(sv_callset_vcf)$FILTER%in%c("PASS",".")]
    
    if (sv_caller == 'lumpy')
    {
      # Read evidence support as a proxy for QUAL
      support <- unlist(info(sv_callset_vcf)$SU)
      fixed(sv_callset_vcf)$QUAL <- support
    } else if (sv_caller == 'delly')
    {
      # Split-read support plus Paired-end read support as a proxy for QUAL
      sr_support <- info(sv_callset_vcf)$SR
      sr_support[is.na(sr_support)] <- 0
      fixed(sv_callset_vcf)$QUAL <-
        sr_support + info(sv_callset_vcf)$PE
    }
    
    bpgr <- breakpointRanges(sv_callset_vcf)
    #begr <- breakendRanges(sv_callset_vcf)
    #gr <- sort(c(bpgr, begr))
    gr <- apply_svtype(bpgr)
    
    # if (sv_caller %in% c('gridss', 'manta')) {
    #   gr <- apply_svtype(gr)
    # }
    
    # Select DEL
    gr <- gr[which(gr$svtype == "DEL")]
    # Remove NAs
    # gr <- gr[!is.na(gr$svLen)]
    # Select DEL >= 50 bp
    gr <- gr[abs(gr$svLen) >= 50]
    gr <- remove_blacklist(gr, confidence_regions_gr, sample_name)
    gr
    
  }

load_truth_set_vcf <-
  function(vcf_file, confidence_regions_gr, sample)
  {
    # vcf_file <- truth_set_file[[sample]]
    sv_callset_vcf <-
      VariantAnnotation::readVcf(vcf_file)
    bpgr <- breakpointRanges(sv_callset_vcf)
    begr <- breakendRanges(sv_callset_vcf)
    gr <- sort(c(bpgr, begr))
    #print(length(gr))
    gr <- gr[which(gr$svtype == "DEL")]
    if (sample %in% c('CHM1', 'CHM13', 'CHM1_CHM13'))
    {
      gr <- gr[gr$svLen >= 50]
      seqlevelsStyle(gr) <- 'NCBI'
      gr <-
        keepSeqlevels(gr, c(1:22, 'X', 'Y'), pruning.mode = "coarse")
    } else if (sample %in% c('NA24385')) {
      gr <- gr[gr$FILTER == 'PASS']
      gr <- gr[gr$svLen <= (-50)]
    } else{
      gr <- gr[gr$svLen <= (-50)]
    }
    
    #gr <- apply_svtype(gr)
    #print(length(gr))
    # gr <- gr[which(gr$svtype == "DEL")]
    #print(length(gr))
    gr <- remove_blacklist(gr, confidence_regions_gr, sample)
    #print(length(remove_blacklist(gr, confidence_regions_gr, sample)))
    gr
  }

load_bedpe <- function(bedpe_file,
                       confidence_regions_gr,
                       sample)
{
  # bedpe_file <- truth_set_file[['NA12878']]
  bpgr <- pairs2breakpointgr(rtracklayer::import(bedpe_file))
  gr <- sort(bpgr)
  gr <- remove_blacklist(gr, confidence_regions_gr, sample)
  gr
}

make_percent <- function(x){
  signif(x*100,digits = 4)
}

load_callsets <- function(sample)
{
  
  gr[[sample]] <- list()
  
  sv_caller_list <- c('gridss', 'manta', 'lumpy', 'delly')
  
  # Load sv_callers results
  for (sv_caller in sv_caller_list)
    #for(sv_caller in c('manta'))
  {
    # sv_caller <- 'gridss'
    print(paste('Loading', sv_caller))
    if (test_sample %in% c('NA12878', 'NA24385'))
    {
      vcf_file <-
        file.path(
          '/Users/lsantuari/Documents/Data/germline/',
          test_sample,
          'SV/Filtered',
          paste(sv_caller, '.vcf', sep = '')
        )
    } else if (test_sample == 'CHM1_CHM13')
    {
      vcf_file <-
        file.path(
          '/Users/lsantuari/Documents/Data/germline/CHM/SV',
          test_sample,
          'Filtered',
          paste(sv_caller, '.vcf', sep = '')
        )
    }
    
    # CHM1_CHM13 mapped with bwa mem
    #vcf_file <- file.path('/Users/lsantuari/Documents/Data/germline/CHM/SV',sample,'Filtered',paste(sv_caller,'.vcf',sep=''))
    print(paste('Loading ', vcf_file, sep = ''))
    gr[[test_sample]][[sv_caller]] <-
      load_sv_caller_vcf(vcf_file, confidence_regions_gr, test_sample, sv_caller)
  }
  
  
  for (c in sv_caller_list)
  {
    print(paste(c, length(gr[[test_sample]][[c]]), 'SVs'))
  }
  
  return(gr)
         
}

bedpe_to_bed <- function(breakpointgr, output.bed)
{
  require(bedr)
  bedpe.file <- StructuralVariantAnnotation::breakpointgr2bedpe(breakpointgr)
  bp1_bed <- bedpe.file[,c(1:3)]
  names(bp1_bed) <- c('chr', 'start', 'end')
  bp2_bed <- bedpe.file[,c(4:6)]
  names(bp2_bed) <- c('chr', 'start', 'end')
  bed_df <- rbind(bp1_bed, bp2_bed)
  bed_df$chr <- as.character(bed_df$chr)
  bed_df$start <- as.integer(bed_df$start)
  bed_df$end <- as.integer(bed_df$end)
  bed_df <- bedr.sort.region(bed_df, check.chr = FALSE)
  write.table(file = output.bed, bed_df, row.names = FALSE, col.names = FALSE, quote = FALSE, sep='\t')
  
}