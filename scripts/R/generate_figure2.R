# remove everything in the environment
rm(list = ls())

setwd('~/Documents/Local_GitHub/sv-channels/scripts/R')

source('./aux_functions.R')

DEBUG <- FALSE

default_param <- 'coverage'
default_param_val <- 60
default_sample <- 'hmz-sv'
default_svtype <- 'INV'

get_default_param_val <- function(parameters.vec, param) {
  if (DEBUG) {
    return(c(default_param_val))
  } else{
    return(parameters.vec[[param]])
  }
}

input.dir <-
  '~/Documents/Projects/GTCG/sv-channels/manuscript/Figures/F2/simulated-data-results/'
encode.list <- '../../data/ENCFF001TDO.bed'
n.regions <- '../../data/reference_N_regions.bed'
default.mode <- 'split_reads'


# create a parser and add command-line arguments
p <-
  arg_parser(
    "Generate figure with the performance (F1-score) of sv-channels
    versus [GRIDSS, Manta, Lumpy, DELLY] for the SV type [DEL, INS, INV, DUP, BND] using SVs simulated with sv-gen"
  )

p <- add_argument(p,
                  "-inputpath",
                  help = "path to the results of the simulated data",
                  type = "character",
                  default = input.dir)
p <- add_argument(
  p,
  "-outputpath",
  help = "path to the output folder",
  type = "character",
  default = file.path(input.dir, '../plots')
)
p <- add_argument(p,
                  "-encodelist",
                  help = "path to the ENCODE blacklist",
                  type = "character",
                  default = encode.list)
p <- add_argument(p,
                  "-nregions",
                  help = "path to the ENCODE blacklist",
                  type = "character",
                  default = n.regions)
p <- add_argument(p,
                  "-mode",
                  help = "sv-channels mode: [split_reads, gridss, manta, delly, lumpy]",
                  type = "character",
                  default = default.mode)

argv <- parse_args(p)

parent_dir <- argv$inputpath
out_dir <- argv$outputpath
encode.blacklist <- argv$encodelist
regions.with.Ns <- argv$nregions
mode <- argv$mode

# TODO: refactor to load from JSON

# list of parameters
parameters <- list()
parameters[['cov']] <- c(5, 10, 15, 30, 45, 60, 75, 90)
parameters[['isize']] <- c(200, 250, 300, 400, 500, 600)
parameters[['rlen']] <- c(36, 50, 75, 100, 150, 250)

# abbreviations
abbrv <- c('cov', 'isize', 'rlen')
names(abbrv) <- c('coverage', 'insert-size', 'read-length')

if (DEBUG) {
  params.list <- c(default_param)
  samples.list <- c(default_sample)
  svtype.list <- c(default_svtype)
} else{
  params.list <- names(abbrv)
  samples.list <- c('hmz-sv', 'htz-sv')
  svtype.list <- c('DEL', 'INS', 'INV', 'DUP', 'CTX')
}

# Loop over three parameter vectors: coverage, insert size and read length
for (param in params.list)
{
  # Loop over the parameter values
  for (param_val in get_default_param_val(parameters, abbrv[param]))
  {
    fld_val <- paste(abbrv[param], param_val, sep = "")
    print(fld_val)
    for (sample.name in samples.list)
    {
      print(sample.name)
      
      files.dir <-
        file.path(parent_dir, param, fld_val, sample.name)
      print(files.dir)
      
      truth_set_file <-
        file.path(files.dir, paste(sample.name, '.proper.vcf', sep = ''))
      
      # SV callsets
      callsets <- list()
      callsets[['GRIDSS']] <- file.path(files.dir, 'gridss.vcf')
      callsets[['Manta']] <- file.path(files.dir, 'manta.vcf')
      callsets[['Lumpy']] <- file.path(files.dir, 'lumpy.vcf')
      callsets[['DELLY']] <- file.path(files.dir, 'delly.vcf')
      #callsets[['ensemble']] <- file.path(files.dir, 'all.vcf')
      
      for (svtype in svtype.list)
      {
        truth_set_svtype <-
          load_vcf(truth_set_file,
                   svtype,
                   sample.name,
                   encode.blacklist,
                   regions.with.Ns)
        
        sv_regions <- list()
        for (caller in c('GRIDSS', 'Manta', 'Lumpy', 'DELLY'))
        {
          if (file.exists(callsets[[caller]]) &
              count_sv_lines(callsets[[caller]]))
          {
            sv_regions[[caller]] <-
              load_vcf(callsets[[caller]],
                       svtype,
                       tolower(caller),
                       encode.blacklist,
                       regions.with.Ns)
          }
        }
        
        callsets[['sv-channels']] <-
          file.path(files.dir,
                    'sv-channels.split_reads.vcf')
        if (file.exists(callsets[['sv-channels']]))
        {
          sv_regions[['sv-channels']] <-
            load_vcf(callsets[['sv-channels']],
                     svtype,
                     'sv-channels',
                     encode.blacklist,
                     regions.with.Ns)
        }
        
        # TODO: use seqs.bed as input file
        # chr.lengths <- c("10" = 135534747, "12" = 133851895)
        #
        # seqlengths(truth_set_svtype) <-
        #   chr.lengths[names(seqlengths(truth_set_svtype))]
        # for (c in names(sv_regions))
        # {
        #   seqlengths(sv_regions[[c]]) <-
        #     chr.lengths[names(seqlengths(sv_regions[[c]]))]
        # }
        #
        # for (n in names(sv_regions)) {
        #   message(paste(n, length(sv_regions[[n]])))
        # }
        
        # exclude callsets with zero calls
        sv_regions <-
          sv_regions[sapply(sv_regions, length) != 0]
        
        # add caller name
        for (c in names(sv_regions))
        {
          sv_regions[[c]]$caller <- c
        }
        
        # collapse breakpoints by merging
        for (c in names(sv_regions))
        {
          # insertions are a special case
          if (svtype == 'INS') {
            sv_regions[[c]]$truth_matches <-
              countBreakpointOverlaps(
                sv_regions[[c]],
                truth_set_svtype,
                # using a smaller margin for insertions
                maxgap = 5,
                # sizemargin cannot be used for insertions
                #sizemargin = 0.25,
                ignore.strand = FALSE,
                #restrictMarginToSizeMultiple = 0.5,
                # countOnlyBest cannot be used for insertions
                countOnlyBest = TRUE
              )
            
          } else {
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
        
        sv_regions_unlisted <-
          unlist(GRangesList(sv_regions))
        
        if (length(sv_regions_unlisted) > 0)
        {
          main.title <-
            paste(param, param_val, sample.name, mode, svtype, sep = ' ')
          print(main.title)
          
          p <- ggplot(
            as.data.frame(sv_regions_unlisted) %>%
              dplyr::select(QUAL, caller, truth_matches) %>%
              dplyr::group_by(caller, QUAL) %>%
              dplyr::summarise(
                calls = n(),
                tp = sum(truth_matches > 0)
              ) %>%
              dplyr::group_by(caller) %>%
              dplyr::arrange(dplyr::desc(QUAL)) %>%
              dplyr::mutate(
                cum_tp = cumsum(tp),
                cum_n = cumsum(calls),
                cum_fp = cum_n - cum_tp,
                precision = cum_tp / cum_n,
                recall = cum_tp / length(truth_set_svtype)
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
              recall = signif(TP / length(truth_set_svtype), digits = 4)
            )
          res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
          res.df$precision <-
            make_percent(res.df$precision)
          res.df$recall <- make_percent(res.df$recall)
          res.df$F1_score <- make_percent(res.df$F1_score)
          
          filename <-
            file.path(output_path, 'precision_recall_plot.csv')
          print(filename)
          write.csv(
            file = filename,
            res.df,
            quote = FALSE,
            row.names = FALSE
          )
        }
      }
    }
  }
}

############

#
# #Summarize results
#
# wide <- data.frame()
#
# for (param in names(abbrv))
# {
#   print(param)
#   for (param_val in sweep[[abbrv[param]]])
#   {
#     fld_val <- paste(abbrv[param], param_val, sep = "")
#     print(fld_val)
#     for (sample.name in c('hmz-sv', 'htz-sv'))
#     {
#       print(sample.name)
#
#       for (svtype in c('DEL', 'INS', 'INV', 'DUP', 'TRA'))
#       {
#         files.dir <-
#           file.path(parent_dir,
#                     'results',
#                     param,
#                     fld_val,
#                     sample.name,
#                     svtype,
#                     mode)
#         print(files.dir)
#
#         filename <- 'precision_recall_plot.csv'
#         if (file.exists(file.path(files.dir, filename)))
#         {
#           df <- read.csv(file.path(files.dir, filename))
#
#           df$caller <- as.vector(df$caller)
#
#           wide <-
#             rbind(
#               wide,
#               cbind(
#                 df,
#                 param = param,
#                 param_val = param_val,
#                 sample.name = sample.name,
#                 svtype = svtype,
#                 mode = mode
#               )
#             )
#         }
#       }
#     }
#   }
# }
# # wide
#
#
# my_labeller <- function(variable, value) {
#   return(paste(value, ' (N=', table(truth_set_bedpe$sourceId)[value], ')', sep =
#                  ''))
# }
#
# for (param in names(abbrv))
# {
#   print(param)
#   for (sample.name in c('hmz-sv', 'htz-sv'))
#   {
#     print(sample.name)
#
#     wide_param <- wide[wide$param == param &
#                          wide$sample.name == sample.name, ]
#     plot <-
#       ggplot(wide_param, aes(caller, F1_score)) +
#       geom_col(aes(fill = caller)) +
#       facet_grid(svtype ~ param_val) +
#       theme(
#         text = element_text(size = 30),
#         axis.text.y = element_text(size = 15),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank()
#       ) +
#       scale_fill_manual(values = cbbPalette)
#     plot <-
#       plot + ggtitle(paste(sample.name, ':', param,
#                            #"TP:green",
#                            #"FP:red",
#                            sep = ''))  # +
#     # geom_text(
#     #   aes(label = FP),
#     #   position = position_dodge(width = 0.9),
#     #   vjust = -0.25,
#     #   color = 'darkred',
#     #   size = 2
#     # ) +
#     # geom_text(
#     #   aes(label = TP),
#     #   position = position_dodge(width = 0.9),
#     #   vjust = -1.50,
#     #   color = 'darkgreen',
#     #   size = 2
#     # )
#     #+ theme(
#     #    panel.grid.major.y = element_line(colour = "grey", linetype="solid")
#     #  )
#
#     ggsave(
#       plot,
#       file = file.path(
#         out_dir,
#         paste(param, sample.name, '_F1_score.png', sep = '_')
#       ),
#       dpi = 600,
#       h = 7,
#       w = 12
#     )
#   }
# }
# 
