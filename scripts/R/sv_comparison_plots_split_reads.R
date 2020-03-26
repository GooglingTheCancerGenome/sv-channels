source("sv_comparison_functions.R")

datasets <- c('NA12878', 'NA24385', 'CHM1_CHM13')

truth_set_file <- list()

truth_set_file[['NA24385']] <-
  'nstd167.GRCh37.variant_call.fixed.vcf.gz'
truth_set_file[['NA12878']] <-
  'Personalis_1000_Genomes_deduplicated_deletions.bedpe'
truth_set_file[['CHM1_CHM13']] <-
  'CHM1_CHM13_pseudodiploid_SVs.vcf.gz'

gr <- list()

for (test_sample in datasets) {
  #test_sample <- 'NA24385'
  
  print(paste(
    'Plotting testing on',
    test_sample,
    sep = ' '
  ))
  
  outDir <-''
  
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
  
  ts <- truth_set[[test_sample]]
  save(ts, file="NA24385_truthset.rda",
      compress = TRUE)
  load("NA24385_truthset.rda")
  dir <- './split_reads_only/NA24385/GRIDSS/'
  write.table(breakpointgr2bedpe(ts), file = file.path(dir, 'NA24385_truthset.bedpe'), sep = '\t', quote=F, row.names = FALSE,
              col.names = FALSE)
  
  ########################################
  
  require(rjson)
  
  #Write split read positions to BEDPE file
  dir <- './split_reads_only/NA24385/IAP'
  json_file <- file.path(dir, 'split_reads.json.gz')
  json_data <- fromJSON(file=json_file)
  my_data <- data.frame()
  
  idx_sr_pos <- 3
  get_value <- function(json, i){
    sapply(json, function(x){x[[i]]})
  }
    my_data <- data.frame(chr1=get_value(json_data[[idx_sr_pos]],1),
                     pos1_s=get_value(json_data[[idx_sr_pos]],2),
                     pos1_e=as.integer(as.vector(get_value(json_data[[idx_sr_pos]],2)))+1,
                     chr1=get_value(json_data[[idx_sr_pos]],3),
                     pos1_s=get_value(json_data[[idx_sr_pos]],4),
                     pos1_e=as.integer(as.vector(get_value(json_data[[idx_sr_pos]],4)))+1
                    )
  write.table(my_data, file = file.path(dir, 'split_reads.bedpe'), sep = '\t', quote=F, row.names = FALSE,
              col.names = FALSE)
  
  bpgr <- pairs2breakpointgr(rtracklayer::import(file.path(dir, 'split_reads.bedpe')))
  split_reads_gr <- sort(bpgr)
  
  overlaps <- findBreakpointOverlaps(
    truth_set[[test_sample]],
    split_reads_gr,
    maxgap = 0,
    ignore.strand = TRUE
  )
  
  length(truth_set[[test_sample]])
  length(truth_set[[test_sample]][unique(queryHits(overlaps))])
  truth_set[[test_sample]] <- truth_set[[test_sample]][unique(queryHits(overlaps))]
  #truth_set[[test_sample]] <- truth_set[[test_sample]][unique(queryHits(overlaps))]
  
  ########################################
  
  gr <- load_callsets(test_sample)

  mode_name <- 'manta_targeted'
  
  modes <- c('complexCNN_v0')
  names(modes) <- c('complexCNN_v0')
             
  for (m in modes) {
    #m <- modes[1]
  for (train_sample in datasets) {
    #train_sample <- 'NA24385'
    
    if( test_sample == train_sample ){
      suf <- ''
      deepsv_name <- paste('DeepSV',train_sample, suf, names(modes)[modes==m], sep='_')
    }else{
      deepsv_name <- paste('DeepSV_trained_on',train_sample, names(modes)[modes==m], sep='_')
    }
    
    print(paste('Loading', deepsv_name))
    
    #if (train_sample == test_sample)
    #{
    
    
    dataDir <-
      paste(
        '/Users/lsantuari/Documents/Processed/Results_DeepSV/channels_33/',m,'/',
        'CNN_CV_train_',
        train_sample,
        '/train_',
        train_sample,
        '_test_',
        test_sample,
        '/predictions/',
        sep = ''
      )
    print(dataDir)
    
    sv_caller_list <- c('gridss', 'manta', 'lumpy', 'delly')
    sv_caller_list <- c(sv_caller_list, deepsv_name)
    
    # DeepSV BEDPE
    bedpe_file <-
      file.path(
        paste(dataDir, 'DeepSV_DEL.svLen.bedpe',
              sep = ""))
    gr[[test_sample]][[deepsv_name]] <-
      load_bedpe(bedpe_file, confidence_regions_gr, test_sample)
    gr[[test_sample]][[deepsv_name]] <-
      gr[[test_sample]][[deepsv_name]][gr[[test_sample]][[deepsv_name]]$NA. >= (50)]

    
  }
  }
    
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
          restrictMarginToSizeMultiple = 0.5,
          # Some callers make duplicate calls and will sometimes report a variant multiple
          # times with slightly different bounds. countOnlyBest prevents these being
          # double-counted as multiple true positives.
          countOnlyBest = TRUE
        )
      
    }
    
    svgr <- unlist(GRangesList(gr[[test_sample]]))
     
        # Plotting Precision and Recall, from StructuralVariantAnnotation vignette:
    # https://bioconductor.org/packages/devel/bioc/vignettes/StructuralVariantAnnotation/inst/doc/vignettes.html
    
    main.title <- c(
      #"NA24385\nNIST_SVs_Integration_v0.6 truth set",
      "NA24385\nWegner2019 truth set",
      "NA12878\nsvclassify truth set",
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
        sec.axis = sec_axis( ~ (.) * length(truth_svgr), name = "true positives")
      ) +
      labs(title = main.title[test_sample])
    
    ggsave(
      filename = paste(
        outDir,
        mode_name,
        '_test_',
        test_sample,
        '.png',
        sep = ''
      ),
      width = 10,
      height = 7,
      units = 'in'
    )
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
        precision = signif(cum_tp / cum_n, digits = 4),
        recall = signif(cum_tp / length(truth_svgr), digits = 4)
      )
    res.df$F1 = with(res.df, 2 * (precision * recall) / (precision + recall))
    
    res.df$precision <- make_percent(res.df$precision)
    res.df$recall <- make_percent(res.df$recall)
    res.df$F1 <- make_percent(res.df$F1)
    
    write.table(
      res.df,
      file = paste(outDir, mode_name, '_test_',
                   test_sample, '_performance_results.csv', sep = ''),
      quote = F,
      row.names = F
    )
    
    #}
}
