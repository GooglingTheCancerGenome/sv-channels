source('/Users/lsantuari/Documents/Local_GitHub/sv-channels/scripts/R/summarize_results_fun.R')

parent_dir <-
  '/Users/lsantuari/Documents/Projects/GTCG/paper_data/real-data-results_noCR2SR_enh/'
out_dir <- file.path(parent_dir, 'results')

regions_for_filtering <-
  '/Users/lsantuari/Documents/Local_GitHub/sv-channels/data/ENCFF001TDO.bed'
N_regions_for_filtering <-
  '/Users/lsantuari/Documents/Local_GitHub/sv-channels/data/Ns.bed'

modes <- c('split_reads', 'gridss', 'manta', 'lumpy', 'delly')
samples.vec <- c('NA12878', 'NA24385', 'CHM1_CHM13')
svtype <-'DEL'

truth_set_list <- c(
  '~/Documents/Processed/SPIDER/notebook-data/NA12878/in/Personalis_1000_Genomes_deduplicated_deletions.bedpe',
  '~/Documents/Processed/SPIDER/notebook-data/NA24385/in/nstd167.GRCh37.variant_call.vcf',
  '~/Documents/Processed/SPIDER/notebook-data/CHM1_CHM13/in/nstd137.GRCh37.variant_call.vcf'
)
names(truth_set_list) <- samples.vec

for (sample.name in samples.vec)
{
  sample.name <- 'NA24385'
  print(sample.name)
  
  if(sample.name %in% c('NA24385', 'CHM1_CHM13'))
  {
  truth_set <- load_vcf(truth_set_list[[sample.name]], svtype, sample.name, regions_for_filtering, N_regions_for_filtering)
  }else if (sample.name == 'NA12878')
  {
  truth_set <- load_bedpe(truth_set_list[[sample.name]], regions_for_filtering, N_regions_for_filtering)
  }
  
  files.dir <-
    file.path(parent_dir, sample.name)
  print(files.dir)
  
  # SV callsets
  callsets <- list()
  callsets[['GRIDSS']] <- file.path(files.dir, 'gridss.vcf')
  callsets[['Manta']] <- file.path(files.dir, 'manta.vcf')
  callsets[['Lumpy']] <- file.path(files.dir, 'lumpy.vcf')
  callsets[['DELLY']] <- file.path(files.dir, 'delly.vcf')
  #callsets[['ensemble_set']] <- file.path(files.dir, 'all.vcf')
  
  sv_regions <- list()
  sv_regions[['GRIDSS']] <-
    load_vcf(callsets[['GRIDSS']], svtype, 'gridss', regions_for_filtering, N_regions_for_filtering)
  sv_regions[['Manta']] <-
    load_vcf(callsets[['Manta']], svtype, 'manta', regions_for_filtering, N_regions_for_filtering)
  sv_regions[['Lumpy']] <-
    load_vcf(callsets[['Lumpy']], svtype, 'lumpy', regions_for_filtering, N_regions_for_filtering)
  sv_regions[['DELLY']] <-
    load_vcf(callsets[['DELLY']], svtype, 'delly', regions_for_filtering, N_regions_for_filtering)
  #sv_regions[['ensemble_set']] <-
  #  load_vcf(callsets[['ensemble_set']], svtype, 'ensemble_set', regions_for_filtering, N_regions_for_filtering)
  
  for (mode in modes)
  {
      # mode <- 'split_reads'
      files.dir <-
        file.path(parent_dir, sample.name)
      print(files.dir)
        mode.name <- paste('sv-channels',mode,sep=':')
        callsets[[mode.name]]  <-
          file.path(files.dir,
                    'sv-channels', mode,
                    paste('win50', '_', svtype, '.bedpe', sep = ''))
        sv_regions[[mode.name]] <-
          load_bedpe(callsets[[mode.name]], regions_for_filtering, N_regions_for_filtering)
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
                truth_set,
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
                truth_set,
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
          paste(sample.name, svtype, 'CV10', sep = ' ')
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
          theme(text = element_text(size=30)) +
          scale_y_continuous(labels = scales::percent) +
          scale_x_continuous(
            labels = scales::percent,
            sec.axis = sec_axis(~ (.) * length(truth_set), name = "true positives")
          ) +
          labs(title = main.title) + scale_fill_manual(values=cbbPalette)
        
        output_path <-
          file.path(out_dir, sample.name, svtype)
        dir.create(output_path, recursive = TRUE)
        filename <-
          file.path(output_path, 'precision_recall_plot.png')
        print(filename)
        ggsave(file = filename, width = 10, height = 5)
        
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
