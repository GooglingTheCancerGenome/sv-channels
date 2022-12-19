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
  vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("PASS", "../..")]

  #if(caller=='truth_set')
  #{
    #vcf_gr <- vcf_gr[which(unlist(info(vcf_gr)$ALGORITHMS)!='manta')]
  #}

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
  print(table(vcf_gr$svtype))
  vcf_gr <- vcf_gr[which(vcf_gr$svtype == svtype)]
  message(length(vcf_gr))
  # Select SVs >= 50 bp
  if (!svtype %in% c('TRA', 'INS'))
  {
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]
  }
  # select only chromosomes 1 to 22
  vcf_gr <-
    vcf_gr[seqnames(vcf_gr) %in% paste('chr', c(1:22), sep = '')]
  vcf_gr <- vcf_gr[vcf_gr$partner %in% names(vcf_gr)]

  return(vcf_gr)
}

res.df <- data.frame()

# create a parser and add command-line arguments
p <-
  arg_parser("Plot benchmark results")
p <-
  add_argument(
    p,
    "--svchannels",
    help = paste("path to the sv-channels VCF file"),
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/validation-220322/validation_HG01053_test_HG01114/sv-channels.DEL.vcf'
  )

p <-
  add_argument(
    p,
    "--manta",
    help = paste("path to the Manta VCF file"),
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/manta_loocv/HG01114/results/manta.vcf'
  )

p <-
  add_argument(
    p,
    "--gridss",
    help = paste("path to the GRIDSS VCF file"),
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/manta_loocv/HG01114/results/gridss.vcf'
  )

p <-
  add_argument(
    p,
    "-t",
    help = "path to the truth set",
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/truth_sets_DEL/HG01114.DEL.vcf'
  )

p <-
  add_argument(p,
               "-n",
               help = "sample name",
               type = "character",
               default = 'HG01114')

p <-
  add_argument(
    p,
    "-f",
    help = "Filename",
    type = "character",
    default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/validation-220322/validation_HG01053_test_HG01114/sv-channels'
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
  if (tools::file_ext(callsets[i, 2]) == 'vcf')
  {
    sv_regions[[callsets[i, 1]]] <-
      load_vcf(callsets[i, 2], svtype, callsets[i, 1])
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
  truth_set <- load_vcf(argv$t, svtype, 'truth_set')
  ignore.strand.flag <- FALSE
} else if (tools::file_ext(argv$t) == 'bedpe')
{
  truth_set <- pairs2breakpointgr(rtracklayer::import(argv$t))
  ignore.strand.flag <- TRUE
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

eval.plot <- ggplot(
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
  theme(text = element_text(size=8)) +
  xlim(0.5, 1) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(
    labels = scales::percent,
    sec.axis = sec_axis(~ (.) * length(truth_set), name = "true positives")
  ) +
  labs(title = main.title) + scale_fill_manual(values=cbbPalette)

ggsave(eval.plot, filename = paste(argv$f, '.plot.png', sep=''))

res.df <- data.frame(
    sample = sampleID,
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
  )


