
source(benchmark_functions.R)

library(argparser, quietly=TRUE)

script.name <- basename(sub(".*=", "", commandArgs()[4])) # script name

# create a parser and add command-line arguments
p <- arg_parser("Generate benchmark plot and table")

p <- add_argument(p, "-i", help="Input CSV file with filenames and paths", type="character")
p <- add_argument(p, "-s", help="Sample name", type="character")
p <- add_argument(p, "-ts", help="Truth set", type="character")
p <- add_argument(p, "-svtype", help="SV type", type="character")
p <- add_argument(p, "-ib", help="BED regions to include", type="character")
p <- add_argument(p, "-eb", help="BED regions to exclude", type="character")
p <- add_argument(p, "-op", help="Output file for plot", type="character")
p <- add_argument(p, "-ot", help="Output file for table", type="character")

# StructuralVariantAnnotation::countBreakpointOverlaps arguments
p <- add_argument(p, "-p", help="SVTYPE=INS if insertion length >= SV length * p",
  type="numeric", default=0.7)
p <- add_argument(p, "-l", help="Minimum SV length to consider", type="int",
  default=50)

# parse the command line arguments
argv <- parse_args(p)
if (is.na(argv$i))
{
    print(p)
    q(status=1)
}

# Read SV callsets from input CSV file
svcallset_table <- read.csv(file = argv$i)
callsets <- list()
names(callsets) <- svcallset_table[,1]
callsets <- svcallset_table[,2]

#Read the truth set
truth_set <- load_vcf(argv$ts, argv$svtype, argv$s, argv$eb, argv$ib)

#Load callsets
sv_regions <- list()
for( c in names(callsets) )
{
    sv_regions[[c]] <- load_vcf(callsets[[c]], argv$svtype, c, argv$eb, argv$ib)
    sv_regions[[c]]$caller <- c
}

# Sanity check on seqlengths
for( c in names(sv_regions) )
{
    if(seqlengths(truth_set)!=seqlengths(sv_regions[[c]]))
    {
        print(paste("seqlengths differ between truth_set and sv_regions for caller", c))
    }
}

#Count overlaps on both breakpoints
for (c in names(sv_regions))
{

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
sv_regions <- unlist(GRangesList(sv_regions))

#Create benchmark plot
main.title <- argv$s

benchmark_plot <- ggplot(
  as.data.frame(sv_regions) %>%
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
    sec.axis = sec_axis( ~ (.) * length(truth_set), name = "true positives")
  ) +
  labs(title = main.title)

ggsave(filename=argv$op, plot=benchmark_plot)

#Create table with performance metrics
res.df <- as.data.frame(sv_regions) %>%
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

write.csv(res.df, file=argv$ot, quote = FALSE, row.names = FALSE)
