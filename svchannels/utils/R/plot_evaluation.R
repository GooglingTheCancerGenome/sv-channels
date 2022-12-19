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
#default_sample <- 'HG00420'
#default_sample <- 'HG00736'

sample_list <- c('HG02924','HG01114','HG01881','HG00420','HG03992','HG01053','HG02018','NA06991')

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

#for(default_sample in sample_list)
#{
default_sample <- 'HG00420'
print(default_sample)
  
res.df <- data.frame()

# create a parser and add command-line arguments
p <-
  arg_parser("Plot benchmark results")
p <-
  add_argument(
    p,
    "-svchannels",
    help = paste("path to the sv-channels VCF file"),
    type = "character",
    #default = '/Users/lsantuari/Documents/GitHub/sv-channels/svchannels/manta_out.vcf',
    default = paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/UMCU_hpc/nested_locso-cv/',default_sample,'.chr1.sv-channels.vcf',sep='')
  )

p <-
  add_argument(
    p,
    "-manta",
    help = paste("path to the Manta VCF file"),
    type = "character",
    default = paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/Manta/',default_sample,'/manta.vcf',sep='')
  )

p <-
  add_argument(
    p,
    "-gridss",
    help = paste("path to the GRIDSS VCF file"),
    type = "character",
    default = paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/GRIDSS2/',default_sample,'.vcf.gz',sep='')
  )

p <-
  add_argument(
    p,
    "-filter",
    help = paste("path to the ENCODE bl + Ns"),
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/reference/exclude.N-and-encodebl.bed'
  )


p <-
  add_argument(
    p,
    "-t",
    help = "path to the truth set",
    type = "character",
    default = paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/truth_sets_DEL/',default_sample,'.DEL.vcf',sep='')
  )

p <-
  add_argument(p,
               "-s",
               help = "sample name",
               type = "character",
               default = default_sample)

p <-
  add_argument(p,
               "-chr",
               help = "validation chromosome",
               type = "character",
               default = 'chr22')

p <-
  add_argument(
    p,
    "-f",
    help = "Filename",
    type = "character",
    default = paste('/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/UMCU_hpc/nested_locso-cv/sv-channels.',default_sample,sep='')
  )

# parse the command line arguments
argv <- parse_args(p)

sampleID <- argv$s
svtype <- 'DEL'

callsets <- data.frame(cbind(
  c('sv-channels', 'manta', 'gridss'),
  c(argv$svchannels, argv$manta, argv$gridss)
))
  
sv_regions <- list()
for (i in 1:nrow(callsets))
{
  if (tools::file_ext(callsets[i, 2]) == 'vcf' | tools::file_ext(callsets[i, 2]) == 'gz')
  {
    sv_regions[[callsets[i, 1]]] <-
      load_vcf(callsets[i, 2], svtype, callsets[i, 1], argv$filter)
  } else if (tools::file_ext(callsets[i, 2]) == 'bedpe')
  {
    sv_regions[[callsets[i, 1]]] <-
      pairs2breakpointgr(rtracklayer::import(callsets[i, 2]))
  }
}
# exclude callsets with zero calls
sv_regions <- sv_regions[sapply(sv_regions, length) != 0]
for (c in names(sv_regions))
{
  sv_regions[[c]]$caller <- c
}

if (tools::file_ext(argv$t) == 'vcf')
{
  truth_set <- load_vcf(argv$t, svtype, 'truth_set', argv$filter)
  ignore.strand.flag <- FALSE
} else if (tools::file_ext(argv$t) == 'bedpe')
{
  truth_set <- pairs2breakpointgr(rtracklayer::import(argv$t))
  ignore.strand.flag <- TRUE
}

#remove chr22, the validation chromosome
truth_set <- truth_set[seqnames(truth_set) == "chr1"]
truth_set <- truth_set[hasPartner(truth_set)]

for (c in names(sv_regions))
{
  sv_regions[[c]] <- sv_regions[[c]][seqnames(sv_regions[[c]]) == "chr1"]
  sv_regions[[c]] <- sv_regions[[c]][hasPartner(sv_regions[[c]])]
}

for (c in names(sv_regions))
{
  if (svtype == 'INS') {
    sv_regions[[c]]$truth_matches <-
      countBreakpointOverlaps(
        sv_regions[[c]],
        truth_set,
        # using a smaller margin for insertions, insertion location should be precise
        maxgap = 100,
        # sizemargin cannot be used for insertions
        # sizemargin = 0.25,
        ignore.strand = ignore.strand.flag,
        restrictMarginToSizeMultiple = 0.5,
        # countOnlyBest cannot be used for insertions
        countOnlyBest = TRUE
      )
    
  } else{
    sv_regions[[c]]$truth_matches <-
      countBreakpointOverlaps(
        sv_regions[[c]],
        truth_set,
        # read pair based callers make imprecise calls.
        # A margin around the call position is required when matching with the truth set
        maxgap = 100,
        # Since we added a maxgap, we also need to restrict the mismatch between the
        # size of the events. We don't want to match a 100bp deletion with a
        # 5bp duplication. This will happen if we have a 100bp margin but don't also
        # require an approximate size match as well
        sizemargin = 0.25,
        ignore.strand = ignore.strand.flag,
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
  paste(argv$s, sep = ' ')
print(main.title)

my.data <- as.data.frame(sv_regions_unlisted) %>%
  dplyr::select(QUAL, caller, truth_matches) %>%
  dplyr::group_by(caller, QUAL) %>%
  dplyr::summarise(calls = n(),
                   tp = sum(truth_matches > 0))
my.data <- as.data.frame(my.data)

summary(my.data$QUAL[my.data$caller=='gridss'])
summary(my.data$QUAL[my.data$caller=='sv-channels'])
summary(my.data$QUAL[my.data$caller=='manta'])

my.data <- as.data.frame(sv_regions_unlisted) %>%
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

my.data %>% dplyr::select(caller, recall) %>% dplyr::group_by(caller) %>%
  dplyr::summarise(max_recall = max(recall))

min.y = min(my.data$precision[my.data$caller=='manta'])
max.y = max(my.data$precision[my.data$caller=='manta'])

hist(my.data$QUAL[my.data$caller=='sv-channels'])

my.x.intercept <- my.data$recall[my.data$caller=='sv-channels'][
  which.min(abs(my.data$QUAL[my.data$caller=='sv-channels']-0.5))]

my.y.intercept <- my.data$precision[my.data$caller=='sv-channels'][
  which.min(abs(my.data$QUAL[my.data$caller=='sv-channels']-0.5))]

eval.plot <- ggplot(
  #my.data[my.data$caller!='gridss',]
  my.data
) +
  aes(x = recall,
      y = precision,
      colour = caller) +
  geom_point() +
  geom_line() +
  theme(text = element_text(size=8),
        plot.title = element_text(size=18),
        legend.text = element_text(size=12),
        legend.title = element_blank(),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"),
        panel.grid = element_line(color = "grey",
                                  size = 0.5,
                                  linetype = 1)
        ) +
  ylim(0,1) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(
    labels = scales::percent,
    sec.axis = sec_axis(~ (.) * length(truth_set), name = "true positives")
  ) +
  # geom_ribbon(data=subset(my.data[my.data$caller!='gridss',], 
  #                         recall > min(my.data$recall[my.data$caller=='manta']) & 
  #                           recall < max(my.data$recall[my.data$caller=='sv-channels'])),
  #             aes(ymin=max.y,
  #                 ymax=precision), ymin=0,
  #             fill = "green", alpha=0.5) +
  labs(title = main.title) + scale_fill_manual(values=cbbPalette) +
  # geom_vline(xintercept = my.data$recall[my.data$caller=='sv-channels'][
  #   which.min(abs(my.data$QUAL[my.data$caller=='sv-channels']))])
  geom_vline(xintercept = my.x.intercept) +
  geom_hline(yintercept = my.y.intercept)
  #annotate(geom="text", x=-0.1, y=my.y.intercept, label=paste(as.integer(my.y.intercept*100),'%',sep='')) +
  #annotate(geom="text", x=my.x.intercept, y=-0.1, label=paste(as.integer(my.x.intercept*100),'%',sep=''))

eval.plot
ggsave(eval.plot, filename = paste(argv$f, '.plot.png', sep=''))

# Filter good quality DELs
c <- 'sv-channels'
sv_regions[[c]] <-sv_regions[[c]][sv_regions[[c]]$QUAL>0.5,] 
sv_regions[[c]] <- sv_regions[[c]][hasPartner(sv_regions[[c]])]

sv_regions_unlisted <- unlist(GRangesList(sv_regions))
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
    recall = signif(TP / length(truth_set), digits = 4)
  )
res.df$F1_score = with(res.df, 2 * (precision * recall) / (precision + recall))
res.df$precision <- make_percent(res.df$precision)
res.df$recall <- make_percent(res.df$recall)
res.df$F1_score <- make_percent(res.df$F1_score)

write.csv(res.df, file = paste(argv$f, '.table.csv', sep=''))

#}