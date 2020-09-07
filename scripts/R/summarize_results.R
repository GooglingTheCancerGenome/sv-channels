source('/Users/lsantuari/Documents/Local_GitHub/sv-channels/scripts/R/summarize_results_fun.R')

parent_dir <-
  '/Users/lsantuari/Documents/Projects/GTCG/paper_data/simulated-data-results'
out_dir <- file.path(parent_dir, 'results')

regions_for_filtering <-
  '/Users/lsantuari/Documents/Local_GitHub/sv-channels/data/ENCFF001TDO.bed'
mode <- 'split_reads'

sweep <- list()
sweep[['cov']] <- c(15, 30, 45, 60, 75, 90)
sweep[['isize']] <- c(150, 200, 250, 300, 400, 500)
#sweep[['isize']] <- c(500)
sweep[['rlen']] <- c(75, 100, 150, 250)
abbrv <- c('cov', 'isize', 'rlen')
names(abbrv) <- c('coverage', 'insert-size', 'read-length')

for (param in c('read-length'))
{
  print(param)
  for (param_val in sweep[[abbrv[param]]])
  {
    fld_val <- paste(abbrv[param], param_val, sep = "")
    print(fld_val)
    for (sample.name in c('hmz-sv', 'htz-sv'))
    {
      print(sample.name)
      
      files.dir <-
        file.path(parent_dir, param, fld_val, sample.name)
      print(files.dir)
      
      truth_set_file <-
        file.path(files.dir, paste(sample.name, '.vcf', sep = ''))
      
      # SV callsets
      callsets <- list()
      callsets[['GRIDSS']] <- file.path(files.dir, 'gridss.vcf')
      callsets[['Manta']] <- file.path(files.dir, 'manta.vcf')
      callsets[['Lumpy']] <- file.path(files.dir, 'lumpy.vcf')
      callsets[['DELLY']] <- file.path(files.dir, 'delly.vcf')
      
      for (svtype in c('DEL', 'INS', 'INV', 'DUP', 'TRA'))
      {
        truth_set_svtype <-
          load_vcf(truth_set_file,
                   svtype,
                   'survivor',
                   regions_for_filtering)
        
        sv_regions <- list()
        sv_regions[['GRIDSS']] <-
          load_vcf(callsets[['GRIDSS']], svtype, 'gridss', regions_for_filtering)
        sv_regions[['Manta']] <-
          load_vcf(callsets[['Manta']], svtype, 'manta', regions_for_filtering)
        sv_regions[['Lumpy']] <-
          load_vcf(callsets[['Lumpy']], svtype, 'lumpy', regions_for_filtering)
        sv_regions[['DELLY']] <-
          load_vcf(callsets[['DELLY']], svtype, 'delly', regions_for_filtering)
        
        callsets[['sv-channels']]  <-
          file.path(files.dir,
                    'sv-channels',
                    paste(mode, '_', svtype, '.bedpe', sep = ''))
        sv_regions[['sv-channels']] <-
          load_bedpe(callsets[['sv-channels']], regions_for_filtering)
        
        seqlengths(truth_set_svtype) <-
          c("10" = 135534747, "12" = 133851895)
        for (c in names(sv_regions))
        {
          seqlengths(sv_regions[[c]]) <- c("10" = 135534747, "12" = 133851895)
        }
        
        for (n in names(sv_regions)) {
          message(paste(n, length(sv_regions[[n]])))
        }
        
        # exclude callsets with zero calls
        sv_regions <- sv_regions[sapply(sv_regions, length) != 0]
        
        for (c in names(sv_regions))
        {
          sv_regions[[c]]$caller <- c
        }
        
        for (c in names(sv_regions))
        {
          if (svtype == 'INS') {
            sv_regions[[c]]$truth_matches <-
              countBreakpointOverlaps(
                sv_regions[[c]],
                truth_set_svtype,
                # using a smaller margin for insertions, insertion location should be precise
                maxgap = 5,
                # sizemargin cannot be used for insertions
                # sizemargin = 0.25,
                ignore.strand = TRUE,
                restrictMarginToSizeMultiple = 0.5
                # countOnlyBest cannot be used for insertions
                # countOnlyBest = TRUE
              )
            
          } else{
            sv_regions[[c]]$truth_matches <-
              countBreakpointOverlaps(
                sv_regions[[c]],
                truth_set_svtype,
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
        }
        
        sv_regions_unlisted <- unlist(GRangesList(sv_regions))
        
        main.title <-
          paste(param, param_val, sample.name, mode, svtype, sep = ' ')
        print(main.title)
        
        ggplot(
          as.data.frame(sv_regions_unlisted) %>%
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
              recall = cum_tp / length(truth_set)
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
            sec.axis = sec_axis(~ (.) * length(truth_set_svtype), name = "true positives")
          ) +
          labs(title = main.title)
        
        output_path <-
          file.path(out_dir, param, fld_val, sample.name, svtype, mode)
        dir.create(output_path, recursive = TRUE)
        filename <-
          file.path(output_path, 'precision_recall_plot.png')
        print(filename)
        ggsave(file = filename)
        
        res.df <- as.data.frame(sv_regions_unlisted) %>%
          dplyr::select(caller, truth_matches) %>%
          dplyr::group_by(caller) %>%
          dplyr::summarise(calls = n(),
                           TP = sum(truth_matches > 0)) %>%
          dplyr::group_by(caller) %>%
          dplyr::mutate(
            FP = calls - TP,
            precision = signif(TP / calls, digits = 4),
            recall = signif(TP / length(truth_set), digits = 4)
          )
        res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
        res.df$precision <- make_percent(res.df$precision)
        res.df$recall <- make_percent(res.df$recall)
        res.df$F1_score <- make_percent(res.df$F1_score)
        
        filename <-
          file.path(output_path, 'precision_recall_plot.csv')
        print(filename)
        write.csv(file = filename,
                  res.df,
                  quote = FALSE,
                  row.names = FALSE)
        
      }
    }
  }
}

############


#Summarize results

wide <- data.frame()

for (param in names(abbrv))
{
  print(param)
  for (param_val in sweep[[abbrv[param]]])
  {
    fld_val <- paste(abbrv[param], param_val, sep = "")
    print(fld_val)
    for (sample.name in c('hmz-sv', 'htz-sv'))
    {
      print(sample.name)
      
      for (svtype in c('DEL', 'INS', 'INV', 'DUP', 'TRA'))
      {
        files.dir <-
          file.path(parent_dir,
                    'results',
                    param,
                    fld_val,
                    sample.name,
                    svtype,
                    mode)
        print(files.dir)
        
        filename <- 'precision_recall_plot.csv'
        df <- read.csv(file.path(files.dir, filename))
        
        df$caller <- as.vector(df$caller)
        
        wide <-
          rbind(
            wide,
            cbind(
              df,
              param = param,
              param_val = param_val,
              sample.name = sample.name,
              svtype = svtype,
              mode = mode
            )
          )
      }
    }
  }
}
# wide


my_labeller <- function(variable, value) {
  return(paste(value, ' (N=', table(truth_set_bedpe$sourceId)[value], ')', sep =
                 ''))
}

for (param in names(abbrv))
{
  print(param)
  for (sample.name in c('hmz-sv', 'htz-sv'))
  {
    print(sample.name)
    
    wide_param <- wide[wide$param == param &
                         wide$sample.name == sample.name, ]
    plot <-
      ggplot(wide_param, aes(caller, F1_score)) + geom_col(aes(fill = caller)) +
      facet_grid(svtype ~ param_val) + theme(
        text = element_text(size = 15),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    plot <-
      plot + ggtitle(paste(
        sample.name,':',param,
        #"TP:green",
        #"FP:red",
        sep = ''
      ))  +
      geom_text(
        aes(label = FP),
        position = position_dodge(width = 0.9),
        vjust = -0.25,
        color = 'darkred',
        size = 2.5
      ) +
      geom_text(
        aes(label = TP),
        position = position_dodge(width = 0.9),
        vjust = -1.50,
        color = 'darkgreen',
        size = 2.5
      )
    
    ggsave(
      plot,
      file = file.path(out_dir, paste(param, sample.name, '_F1_score.png', sep = '_')),
      dpi = 600,
      h = 7,
      w = 12
    )
  }
}