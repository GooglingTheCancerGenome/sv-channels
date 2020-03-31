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
           sample,
           sv_caller)
  {
    # vcf_file <- truth_set_file[[sample]]
    sv_callset_vcf <-
      VariantAnnotation::readVcf(vcf_file)
    
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
    begr <- breakendRanges(sv_callset_vcf)
    gr <- sort(c(bpgr, begr))
    if (sv_caller %in% c('gridss', 'manta')) {
      gr <- apply_svtype(gr)
    }
    # Select DEL
    gr <- gr[which(gr$svtype == "DEL")]
    # Remove NAs
    # gr <- gr[!is.na(gr$svLen)]
    # Select DEL >= 50 bp
    gr <- gr[abs(gr$svLen) >= 50]
    gr <- remove_blacklist(gr, confidence_regions_gr, sample)
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
  sv_callset_bedpe <- rtracklayer::import(bedpe_file)
  bpgr <- pairs2breakpointgr(sv_callset_bedpe)
  gr <- sort(bpgr)
  gr <- remove_blacklist(gr, confidence_regions_gr, sample)
  gr
}

make_percent <- function(x){
  round(x*100,digits = 1)
}

datasets <- c('NA12878', 'NA24385', 'CHM1_CHM13')

truth_set_file <- list()
truth_set_file[['NA24385']] <-
  '/Users/lsantuari/Documents/Data/germline/NA24385/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz'
truth_set_file[['NA12878']] <-
  '/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe'
truth_set_file[['CHM1']] <-
  '/Users/lsantuari/Documents/Data/germline/CHM/Huddleston2016/structural_variants/CHM1_SVs.annotated.vcf.gz'
truth_set_file[['CHM13']] <-
  '/Users/lsantuari/Documents/Data/germline/CHM/Huddleston2016/structural_variants/CHM13_SVs.annotated.vcf.gz'
truth_set_file[['CHM1_CHM13']] <-
  '/Users/lsantuari/Documents/Data/germline/CHM/Huddleston2016/structural_variants/CHM1_CHM13_pseudodiploid_SVs.vcf.gz'

gr <- list()

for (test_sample in datasets) {
  
  # test_sample <- 'NA12878'

  gr[[test_sample]] <- list()
  
  truth_set <- list()
  
  #Load blacklist
  confidence_regions_gr <- load_confidence_regions(test_sample)
  
  if (test_sample %in% c('NA24385', 'CHM1', 'CHM13', 'CHM1_CHM13'))
  {
    
    truth_set[[test_sample]] <-
      load_truth_set_vcf(truth_set_file[[test_sample]],
                         confidence_regions_gr,
                         test_sample)
    
  } else if (test_sample == 'NA12878') {
    
    truth_set[[test_sample]] <-
      load_bedpe(truth_set_file[[test_sample]], 
                 confidence_regions_gr, 
                 test_sample)
  
    }
  
  for (train_sample in datasets) {
    
    # train_sample <- 'NA12878'
    
    #if (train_sample == test_sample)
    #{
      
      print(paste('Plotting training on',train_sample,'testing on',test_sample,sep = ' '))
      outDir <- paste('/Users/lsantuari/Documents/Processed/Results_DeepSV/complexCNN_v0/',
                      'CNN_CV_train_',train_sample,'/train_',train_sample,'_test_',test_sample,'/predictions/', sep='') 
      print(outDir)
      
      # sample <- 'NA12878'
      # vcf_file <- truth_set_file[[sample]]
      # sv_callset_vcf <-
      #   VariantAnnotation::readVcf(vcf_file)
      # bpgr <- breakpointRanges(sv_callset_vcf)
      # begr <- breakendRanges(sv_callset_vcf)
      
      if(test_sample =='CHM1_CHM13'){
        sv_caller_list <- c('gridss', 'manta','lumpy')
      }else{
        sv_caller_list <- c('gridss', 'manta', 'lumpy', 'delly')
      }
      
      # Load sv_callers results
      for (sv_caller in sv_caller_list)
        #for(sv_caller in c('manta'))
      {
        # sv_caller <- 'gridss'
        print(paste('Loading', sv_caller))
        # NA24385 mapped with novoalign
        #vcf_file <- file.path('/Users/lsantuari/Documents/Data/germline',paste(sample,'novoalign',sep='_'),'SV/Filtered',paste(sv_caller,'.vcf',sep=''))
        # NA24385 mapped with bwa mem
        # vcf_file <-
        #   file.path(
        #     '/Users/lsantuari/Documents/Data/germline',
        #     sample,
        #     'SV/Filtered',
        #     paste(sv_caller, '.vcf', sep = '')
        #   )
        # NA12878 mapped with bwa mem
        #vcf_file <- file.path('/Users/lsantuari/Documents/Data/germline/trio',sample,'SV/Filtered',paste(sv_caller,'.flt.vcf',sep=''))
        
        if(test_sample %in% c('NA12878', 'NA24385'))
        {
        vcf_file <-
          file.path(
            '/Users/lsantuari/Documents/Data/germline/',
            test_sample,
            'SV/Filtered',
            paste(sv_caller, '.vcf', sep = '')
          )
        }else if(test_sample == 'CHM1_CHM13')
        {
          vcf_file <- file.path('/Users/lsantuari/Documents/Data/germline/CHM/SV',test_sample,'Filtered',paste(sv_caller,'.vcf',sep=''))
        }
        
        # CHM1_CHM13 mapped with bwa mem
        #vcf_file <- file.path('/Users/lsantuari/Documents/Data/germline/CHM/SV',sample,'Filtered',paste(sv_caller,'.vcf',sep=''))
        print(paste('Loading ',vcf_file,sep=''))
        gr[[test_sample]][[sv_caller]] <- load_sv_caller_vcf(vcf_file, confidence_regions_gr, test_sample, sv_caller)
      }
      
      for(c in sv_caller_list)
      {
        print(paste(c, length(gr[[test_sample]][[c]]), 'SVs'))
      }
      
      sv_caller_list <- c(sv_caller_list, 'deepsv')
      bedpe_file <-
        file.path(
          #paste('/Users/lsantuari/Documents/Processed/channel_maker_output/',sample,'/cnn/CNN_CV_100919/DeepSV_DEL.svLen.bedpe',sep="")
          paste(
            outDir, 'DeepSV_DEL.svLen.bedpe',
            sep = ""
          )
        )
      gr[[test_sample]][['deepsv']] <-
        load_bedpe(bedpe_file, confidence_regions_gr, test_sample)
      gr[[test_sample]][['deepsv']] <-
        gr[[test_sample]][['deepsv']][gr[[test_sample]][['deepsv']]$NA. >= (50)]
      
      
      # Only consider Chr1 to ChrX
      for (sv_caller in sv_caller_list)
      {
        print(sv_caller)
        print(length(gr[[test_sample]][[sv_caller]]))
        print(length(gr[[test_sample]][[sv_caller]][seqnames(gr[[test_sample]][[sv_caller]]) !=
                                                 'Y']))
        # gr[[sample]][[sv_caller]] <- gr[[sample]][[sv_caller]][seqnames(gr[[sample]][[sv_caller]])!='Y']
      }
      
      truth_svgr <- truth_set[[test_sample]]
      for (sv_caller in sv_caller_list)
      {
        gr[[test_sample]][[sv_caller]]$caller <- sv_caller
      }
      
      for (svcaller in sv_caller_list)
      {
        gr[[test_sample]][[svcaller]]$truth_matches <-
          countBreakpointOverlaps(
            gr[[test_sample]][[svcaller]],
            truth_svgr,
            # read pair based callers make imprecise calls.
            # A margin around the call position is required when matching with the truth set
            maxgap = 200,
            # Since we added a maxgap, we also need to restrict the mismatch between the
            # size of the events. We don't want to match a 100bp deletion with a
            # 5bp duplication. This will happen if we have a 100bp margin but don't also
            # require an approximate size match as well
            sizemargin = 0.25,
            ignore.strand = TRUE,
            # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
            # by restricting the margin based on the size of the event, we can make sure
            # that simple events actually do overlap
            restrictMarginToSizeMultiple =
              0.5,
            # Some callers make duplicate calls and will sometimes report a variant multiple
            # times with slightly different bounds. countOnlyBest prevents these being
            # double-counted as multiple true positives.
            countOnlyBest = TRUE
          )
      }
      
      if(test_sample == 'CHM1_CHM13')
      {
        svgr <- c(gr[[test_sample]][['gridss']],
                  gr[[test_sample]][['manta']],
                  gr[[test_sample]][['lumpy']],
                  gr[[test_sample]][['deepsv']])
      }else{
        svgr <- c(gr[[test_sample]][['gridss']],
                  gr[[test_sample]][['manta']],
                  gr[[test_sample]][['lumpy']],
                  gr[[test_sample]][['delly']],
                  gr[[test_sample]][['deepsv']])
      }
      
      # Plotting Precision and Recall, from StructuralVariantAnnotation vignette:
      # https://bioconductor.org/packages/devel/bioc/vignettes/StructuralVariantAnnotation/inst/doc/vignettes.html
      
      main.title <- c(
        "NA24385\nNIST_SVs_Integration_v0.6 truth set",
        "NA12878\nsv_classify truth set",
        "CHM1_CHM13\nHuddleston2016 truth set"
      )
      names(main.title) <- c("NA24385", "NA12878", "CHM1_CHM13")
      
      suppressPackageStartupMessages(require(dplyr))
      suppressPackageStartupMessages(require(ggplot2))
      ggplot(
        as.data.frame(svgr) %>%
          dplyr::select(QUAL, caller, truth_matches) %>%
          dplyr::group_by(caller, QUAL) %>%
          dplyr::summarise(calls = n(),
                           tp = sum(truth_matches > 0)) %>%
          dplyr::group_by(caller) %>%
          dplyr::arrange(dplyr::desc(QUAL)) %>%
          dplyr::mutate(
            cum_tp = cumsum(tp),
            cum_n = cumsum(calls),
            cum_fp = cum_n - cum_tp,
            precision = cum_tp / cum_n,
            recall = cum_tp / length(truth_svgr)
          )
      ) +
        aes(x = recall,
            y = precision,
            colour = caller) +
        geom_point() +
        geom_line() +
        scale_y_continuous(labels = scales::percent) +
        scale_x_continuous(
          labels = scales::percent,
          sec.axis = sec_axis(~ (.) * length(truth_svgr), name = "true positives")
        ) +
        labs(title = main.title[test_sample])
      ggsave(filename = paste(outDir,'train_',train_sample,'_test_',test_sample,'.png',sep=''), width = 6, height = 5, units='in')
      res.df <- as.data.frame(svgr) %>%
        dplyr::select(caller, truth_matches) %>%
        dplyr::group_by(caller) %>%
        dplyr::summarise(calls = n(),
                         tp = sum(truth_matches > 0)) %>%
        dplyr::group_by(caller) %>%
        dplyr::mutate(
          cum_tp = cumsum(tp),
          cum_n = cumsum(calls),
          cum_fp = cum_n - cum_tp,
          precision = round(cum_tp / cum_n,digits = 1),
          recall = round(cum_tp / length(truth_svgr), digits = 1)
        )
      res.df$F1 = with(res.df, 2 * (precision * recall) / (precision + recall))
      
      res.df$precision <- make_percent(res.df$precision)
      res.df$recall <- make_percent(res.df$recall)
      res.df$F1 <- make_percent(res.df$F1)
      
      write.table(res.df, file=paste(outDir, 'performance_results.csv',sep=''), quote=F, row.names = F)
      
    #}
  }
}