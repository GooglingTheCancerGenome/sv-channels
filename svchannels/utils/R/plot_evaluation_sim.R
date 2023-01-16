rm(list = ls())
suppressPackageStartupMessages(require(argparser))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))

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

load_bed <- function(bed_file)
{
  bed_regions <- rtracklayer::import(bed_file)
  # set NCBI seqlevels
  #seqlevelsStyle(bed_regions) <- "NCBI"
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
  
  message(paste('Loading', length(vcf_gr), svtype, 'calls for', caller, sep =
                  " "))
  
  # Select SVs >= 50 bp
  if (!svtype %in% c('CTX', 'INS'))
  {
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]
  }
  
  # Filter regions
  vcf_gr <-
    filter_regions('ENCODE blacklist', vcf_gr, load_bed(filter.regions), mode = 'remove')
  
  return(vcf_gr)
}

filter_regions <-
  function(filter.name,
           regions_to_filter,
           ref_regions,
           mode = 'remove')
  {
    print(length(regions_to_filter))
    if (mode == 'keep')
    {
      result <-
        regions_to_filter[overlapsAny(regions_to_filter, ref_regions) &
                            overlapsAny(partner(regions_to_filter), ref_regions), ]
    } else if (mode == 'remove') {
      result <- regions_to_filter[!(
        overlapsAny(regions_to_filter, ref_regions) |
          overlapsAny(partner(regions_to_filter), ref_regions)
      ), ]
    }
    message(paste(length(result), 'calls after filtering for', filter.name, sep =
                    " "))
    return(result)
  }

svtype <- 'DEL'

params_vec <- c("coverage", "insert-size", "read-length")
gt_list <- c("hmz", "htz")
# This is the file "exclude.N-and-encodebl.hg19.bed" with only chromosome 10 and 12
excl_list <-
  "/Users/lsantuari/Documents/GitHub/sv-channels/data/exclude.N-and-encodebl.10_12_hs37d5.bed"

# directory with results
res_dir <-
  "/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/UMCU_hpc/results/workflow_sim/results_sim"
# directory with GRIDSS callsets
gridss_dir <- "/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/tmp/simulated-data-results_v0"

truth_set <- list()
for (t in params_vec) {
  print(t)
  truth_set[[t]] <- list()
  for (gt in gt_list) {
    print(gt)
    vcf_file <-
      file.path(res_dir,
                "truth_set_sim",
                t,
                paste(gt, "-sv.proper.vcf", sep = ""))
    truth_set[[t]][[gt]] <-
      load_vcf(vcf_file, svtype, 'truth_set', excl_list)
  }
}

results <- data.frame()

# coverage
p <- "coverage"
cov_range <- c("5", "10", "15", "30", "45", "60", "75", "90")
isize <- "500"
rlen <- "150"
for (cov in cov_range) {
  for (gt in gt_list) {
    sv_callset <- list()
    for (caller in c("manta", "gridss", "sv-channels")) {
      if(caller != "gridss"){
      vcf_file <- file.path(
        res_dir,
        p,
        paste("c", cov, sep = ""),
        paste("i", isize, sep = ""),
        paste("r", rlen, sep = ""),
        gt,
        gsub('-', '', caller),
        paste(caller, ".vcf", sep = "")
      )
      }else{
        vcf_file <- file.path(
          gridss_dir,
          p,
          paste("cov", cov, sep = ""),
          paste(gt,'-sv',sep=''),
          paste(caller, ".vcf", sep = "")
        )
      }
      
      sv_callset[[caller]] <- load_vcf(vcf_file, svtype, caller, excl_list)
      
      sv_callset[[caller]]$truth_matches <-
        countBreakpointOverlaps(
          sv_callset[[caller]],
          truth_set[[p]][[gt]],
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
          restrictMarginToSizeMultiple = 0.5,
          # Some callers make duplicate calls and will sometimes report a variant multiple
          # times with slightly different bounds. countOnlyBest prevents these being
          # double-counted as multiple true positives.
          countOnlyBest = TRUE
        )
      
    }
    sv_callset <- sv_callset[sapply(sv_callset, length) != 0]
    for (c in names(sv_callset))
    {
      sv_callset[[c]]$caller <- c
    }
    
    caller <- 'sv-channels'
    sv_callset[[caller]] <- sv_callset[[caller]][sv_callset[[caller]]$QUAL > 0.5, ]
    sv_callset[[caller]] <- sv_callset[[caller]][hasPartner(sv_callset[[caller]])]
    
    sv_regions_unlisted <- unlist(GRangesList(sv_callset))
    res.df <- as.data.frame(sv_regions_unlisted) %>%
      dplyr::select(caller, truth_matches, QUAL) %>%
      dplyr::group_by(caller) %>%
      dplyr::arrange(QUAL) %>%
      dplyr::summarise(calls = n(),
                       TP = sum(truth_matches > 0)) %>%
      dplyr::group_by(caller) %>%
      dplyr::mutate(
        FP = calls - TP,
        precision = signif(TP / calls, digits = 4),
        recall = signif(TP / length(truth_set[[p]][[gt]]), digits = 4)
      )
    res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
    res.df$precision <- make_percent(res.df$precision)
    res.df$recall <- make_percent(res.df$recall)
    res.df$F1_score <- make_percent(res.df$F1_score)
    results <- 
      rbind(results,
            cbind(data.frame(parameter=p, genotype=gt, coverage=cov, 
                             insert_size=isize, read_length=rlen), as.data.frame(res.df))
            )
    
  }
}

# insert-size
p <- "insert-size"
isize_list <- c("200", "250", "300", "400", "500", "600")
cov <- "30"
rlen <- "150"
for (isize in isize_list) {
  for (gt in gt_list) {
    sv_callset <- list()
    for (caller in c("manta", "gridss", "sv-channels")) {
      if(caller != "gridss"){
      vcf_file <- file.path(
        res_dir,
        p,
        paste("c", cov, sep = ""),
        paste("i", isize, sep = ""),
        paste("r", rlen, sep = ""),
        gt,
        gsub('-', '', caller),
        paste(caller, ".vcf", sep = "")
      )
      }else{
        vcf_file <- file.path(
          gridss_dir,
          p,
          paste("isize", isize, sep = ""),
          paste(gt,'-sv',sep=''),
          paste(gt,'-sv',sep=''),
          paste(caller, "_out", sep = ""),
          paste(caller, ".vcf", sep = "")
        )
      }
      if(file.exists(vcf_file))
      {
      sv_callset[[caller]] <- load_vcf(vcf_file, svtype, caller, excl_list)
      
      sv_callset[[caller]]$truth_matches <-
        countBreakpointOverlaps(
          sv_callset[[caller]],
          truth_set[[p]][[gt]],
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
          restrictMarginToSizeMultiple = 0.5,
          # Some callers make duplicate calls and will sometimes report a variant multiple
          # times with slightly different bounds. countOnlyBest prevents these being
          # double-counted as multiple true positives.
          countOnlyBest = TRUE
        )
      }
    }
    sv_callset <- sv_callset[sapply(sv_callset, length) != 0]
    for (c in names(sv_callset))
    {
      sv_callset[[c]]$caller <- c
    }
    
    caller <- 'sv-channels'
    sv_callset[[caller]] <- sv_callset[[caller]][sv_callset[[caller]]$QUAL > 0.5, ]
    sv_callset[[caller]] <- sv_callset[[caller]][hasPartner(sv_callset[[caller]])]
    
    sv_regions_unlisted <- unlist(GRangesList(sv_callset))
    res.df <- as.data.frame(sv_regions_unlisted) %>%
      dplyr::select(caller, truth_matches, QUAL) %>%
      dplyr::group_by(caller) %>%
      dplyr::arrange(QUAL) %>%
      dplyr::summarise(calls = n(),
                       TP = sum(truth_matches > 0)) %>%
      dplyr::group_by(caller) %>%
      dplyr::mutate(
        FP = calls - TP,
        precision = signif(TP / calls, digits = 4),
        recall = signif(TP / length(truth_set[[p]][[gt]]), digits = 4)
      )
    res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
    res.df$precision <- make_percent(res.df$precision)
    res.df$recall <- make_percent(res.df$recall)
    res.df$F1_score <- make_percent(res.df$F1_score)
    results <- 
      rbind(results,
            cbind(data.frame(parameter=p, genotype=gt, coverage=cov, 
                             insert_size=isize, read_length=rlen), as.data.frame(res.df))
      )
    
  }
}


# read-length
p <- "read-length"
rlen_list <- c("36", "50", "75", "100", "150", "250")
cov <- "30"
isize <- "500"
for (rlen in rlen_list) {
  for (gt in gt_list) {
    sv_callset <- list()
    for (caller in c("manta", "gridss", "sv-channels")) {
      if(caller != "gridss"){
      vcf_file <- file.path(
        res_dir,
        p,
        paste("c", cov, sep = ""),
        paste("i", isize, sep = ""),
        paste("r", rlen, sep = ""),
        gt,
        gsub('-', '', caller),
        paste(caller, ".vcf", sep = "")
      )
      }else{
        vcf_file <- file.path(
          gridss_dir,
          p,
          paste("rlen", rlen, sep = ""),
          paste(gt,'-sv',sep=''),
          paste(caller, ".vcf", sep = "")
        )
      }
      if(file.exists(vcf_file))
      {
        sv_callset[[caller]] <- load_vcf(vcf_file, svtype, caller, excl_list)
        
        sv_callset[[caller]]$truth_matches <-
          countBreakpointOverlaps(
            sv_callset[[caller]],
            truth_set[[p]][[gt]],
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
            restrictMarginToSizeMultiple = 0.5,
            # Some callers make duplicate calls and will sometimes report a variant multiple
            # times with slightly different bounds. countOnlyBest prevents these being
            # double-counted as multiple true positives.
            countOnlyBest = TRUE
          )
      }
    }
    sv_callset <- sv_callset[sapply(sv_callset, length) != 0]
    for (c in names(sv_callset))
    {
      sv_callset[[c]]$caller <- c
    }
    
    caller <- 'sv-channels'
    if(caller %in% names(sv_callset))
    {
    sv_callset[[caller]] <- sv_callset[[caller]][sv_callset[[caller]]$QUAL > 0.5, ]
    sv_callset[[caller]] <- sv_callset[[caller]][hasPartner(sv_callset[[caller]])]
    }
    if(length(sv_callset)>0)
    {
    sv_regions_unlisted <- unlist(GRangesList(sv_callset))
    res.df <- as.data.frame(sv_regions_unlisted) %>%
      dplyr::select(caller, truth_matches, QUAL) %>%
      dplyr::group_by(caller) %>%
      dplyr::arrange(QUAL) %>%
      dplyr::summarise(calls = n(),
                       TP = sum(truth_matches > 0)) %>%
      dplyr::group_by(caller) %>%
      dplyr::mutate(
        FP = calls - TP,
        precision = signif(TP / calls, digits = 4),
        recall = signif(TP / length(truth_set[[p]][[gt]]), digits = 4)
      )
    res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
    res.df$precision <- make_percent(res.df$precision)
    res.df$recall <- make_percent(res.df$recall)
    res.df$F1_score <- make_percent(res.df$F1_score)
    results <- 
      rbind(results,
            cbind(data.frame(parameter=p, genotype=gt, coverage=cov, 
                             insert_size=isize, read_length=rlen), as.data.frame(res.df))
      )
    }
  }
}

#write.csv(results, file = file.path(res_dir, 'results.csv'))

results <- read.csv(file = file.path(res_dir, 'results.csv'))

df1 <- results[results$parameter=="insert-size",c("parameter","genotype","insert_size", "caller", "precision", "recall", "F1_score")]
df2 <- results[results$parameter=="read-length",c("parameter","genotype","read_length", "caller", "precision", "recall", "F1_score")]
df3 <- results[results$parameter=="coverage",c("parameter","genotype","coverage", "caller", "precision", "recall", "F1_score")]
names(df1)[3] <- "value"
names(df2)[3] <- "value"
names(df3)[3] <- "value"

results_prec <- rbind(df1,df2,df3)

results$coverage <- as.factor(results$coverage)
results$insert_size <- as.factor(results$insert_size)
results$read_length <- as.factor(results$read_length)

results <- results[results$caller=='sv-channels',]
results$genotype[results$genotype=="hmz"] <- "homozygous"
results$genotype[results$genotype=="htz"] <- "heterozygous"

library(ggplot2)
bp1 <- ggplot(results[results$parameter=="coverage",], aes(x=coverage, y=recall, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp1 <- bp1 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.x=element_blank()) + scale_fill_manual(values=c("#619CFF"))

bp2 <- ggplot(results[results$parameter=="insert-size",], aes(x=insert_size, y=recall, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp2 <- bp2 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.y=element_blank(),
                                                                 axis.title.x=element_blank()) + scale_fill_manual(values=c("#619CFF"))

bp3 <- ggplot(results[results$parameter=="read-length",], aes(x=read_length, y=recall, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.title.y = element_blank(),
                                                      axis.title.x=element_blank(),
                                                      legend.position = "none") + scale_fill_manual(values=c("#619CFF"))
bp3 <- bp3 + facet_grid(vars(genotype), vars(parameter)) 

require(ggpubr)
my.bp <- ggarrange(bp1, bp2, bp3, 
          labels = c("A", "B", "C"),
          ncol = 3, nrow = 1)
          #common.legend = TRUE, legend="bottom")
my.bp
ggsave(plot=my.bp,
       width = 178, height = as.integer(178/2), dpi = 300,
       filename = file.path(res_dir, "recall.sv-channels.png"), units = "mm")

bp1 <- ggplot(results[results$parameter=="coverage",], aes(x=coverage, y=precision, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp1 <- bp1 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.x=element_blank()
                                                                 ) + scale_fill_manual(values=c("#619CFF"))

bp2 <- ggplot(results[results$parameter=="insert-size",], aes(x=insert_size, y=precision, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp2 <- bp2 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.x=element_blank(),
                                                                 axis.title.y=element_blank()) + scale_fill_manual(values=c("#619CFF"))

bp3 <- ggplot(results[results$parameter=="read-length",], aes(x=read_length, y=precision, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.title.y = element_blank(),
                                                      axis.title.x=element_blank(),
                                                      legend.position = "none") + scale_fill_manual(values=c("#619CFF"))
bp3 <- bp3 + facet_grid(vars(genotype), vars(parameter)) 

my.bp <- ggarrange(bp1, bp2, bp3, 
                   labels = c("A", "B", "C"),
                   ncol = 3, nrow = 1)
                   #common.legend = TRUE, legend="bottom")
my.bp
ggsave(plot=my.bp,
       width = 178, height = as.integer(178/2), dpi = 300,
       filename = file.path(res_dir, "precision.sv-channels.png"), units = "mm")

bp1 <- ggplot(results[results$parameter=="coverage",], aes(x=coverage, y=F1_score, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp1 <- bp1 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.x=element_blank()
                                                                 ) + scale_fill_manual(values=c("#619CFF"))

bp2 <- ggplot(results[results$parameter=="insert-size",], aes(x=insert_size, y=F1_score, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity")
bp2 <- bp2 + facet_grid(vars(genotype), vars(parameter)) + theme(legend.position = "none",
                                                                 axis.title.x=element_blank(),
                                                                 axis.title.y=element_blank()) + scale_fill_manual(values=c("#619CFF"))

bp3 <- ggplot(results[results$parameter=="read-length",], aes(x=read_length, y=F1_score, group=caller, fill=caller)) + 
  geom_bar(position="dodge", stat="identity") + theme(axis.title.y = element_blank(),
                                                      axis.title.x=element_blank(),
                                                      legend.position = "none") + scale_fill_manual(values=c("#619CFF"))
bp3 <- bp3 + facet_grid(vars(genotype), vars(parameter)) 

my.bp <- ggarrange(bp1, bp2, bp3, 
                   labels = c("A", "B", "C"),
                   ncol = 3, nrow = 1)
                   #common.legend = TRUE, legend="bottom")
my.bp
ggsave(plot=my.bp,
       width = 178, height = as.integer(178/2), dpi = 300,
       filename = file.path(res_dir, "F1_score.sv-channels.png"), units = "mm")
