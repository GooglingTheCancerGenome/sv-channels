rm(list=ls())
suppressPackageStartupMessages(require(argparser))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
  vcf_gr <- vcf_gr[seqnames(vcf_gr) %in% c(1:22)]
  vcf_gr <- vcf_gr[vcf_gr$partner %in% names(vcf_gr)]

  return(vcf_gr)
}

# create a parser and add command-line arguments
p <-
  arg_parser("Plot benchmark results")
p <-
  add_argument(
    p,
    "-i",
    help = paste("path to the CSV file with SV callsets (VCF)"),
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/Figures/F3/callsets_NA24385.csv'
  )
p <-
  add_argument(p,
               "-t",
               help = "path to the truth set",
               type = "character",
               default = '/Users/lsantuari/Documents/Processed/SPIDER/split_reads/NA24385/nstd167.GRCh37.variant_call.vcf')
p <-
  add_argument(p,
               "-s",
               help = "sample name",
               type = "character",
               default = 'NA24385')

p <-
  add_argument(p,
               "-sr",
               help = "path to split reads in BEDPE format",
               type = "character",
               default = '/Users/lsantuari/Documents/Processed/SPIDER/split_reads/NA24385/split_reads_real.bedpe')

p <-
  add_argument(p,
               "-b",
               help = "consider only split read positions?",
               type = "boolean",
               default = FALSE)

# parse the command line arguments
argv <- parse_args(p)
svtype <- 'DEL'
callsets <- read.csv(argv$i, header=FALSE)
sv_regions <- list()
for(i in 1:nrow(callsets))
{
  if(tools::file_ext(callsets[i,2]) == 'vcf')
  {
  sv_regions[[callsets[i,1]]] <- load_vcf(callsets[i,2], 'DEL', callsets[i,1])
  }else if(tools::file_ext(callsets[i,2]) == 'bedpe')
  {
  sv_regions[[callsets[i,1]]] <- pairs2breakpointgr(rtracklayer::import(callsets[i,2]))
  }
}
# exclude callsets with zero calls
sv_regions <- sv_regions[sapply(sv_regions, length) != 0]
for (c in names(sv_regions))
{
  sv_regions[[c]]$caller <- c
}

if(tools::file_ext(argv$t) == 'vcf')
{
  truth_set <- load_vcf(argv$t, 'DEL', 'truth_set')
}else if(tools::file_ext(argv$t) == 'bedpe')
{
  truth_set <- pairs2breakpointgr(rtracklayer::import(argv$t))
}

if(argv$b)
{
  
# Consider only SVs with split read positions?
split_reads <- pairs2breakpointgr(rtracklayer::import(argv$sr))

for (c in names(sv_regions))
{
print(c)
myHits <- findBreakpointOverlaps(
  sv_regions[[c]],
  split_reads,
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = TRUE,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)
print(length(unique(queryHits(myHits))))
sv_regions[[c]] <- sv_regions[[c]][unique(queryHits(myHits)),]
sv_regions[[c]] <- sv_regions[[c]][sv_regions[[c]]$partner %in% names(sv_regions[[c]])]
}

myHits <- findBreakpointOverlaps(
  truth_set,
  split_reads,
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = TRUE,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)
truth_set <- truth_set[unique(queryHits(myHits)),]

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
        maxgap = 100,
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
  paste(argv$s, svtype, sep = ' ')
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

filename <- 'precision_recall_plot.png'
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

filename <- 'precision_recall_plot.csv'

write.csv(file = filename,
          res.df,
          quote = FALSE,
          row.names = FALSE)
