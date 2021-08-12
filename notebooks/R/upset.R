source('~/Documents/Local_GitHub/sv-channels-latest/notebooks/R/aux_functions.R')

sampleID <- 'NA06991'

truth_set_file <- paste(
  '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/truth_sets/',
  sampleID,
  '.DEL.vcf',
  sep = ''
)

callsets_files <- paste(
  '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/results/manta_loocv/',
  sampleID,
  '/callsets_',
  sampleID,
  '.csv',
  sep = ''
)

callsets <- read.csv(callsets_files, header = FALSE)
sv_regions <- list()
for (i in 1:nrow(callsets))
{
  if (tools::file_ext(callsets[i, 2]) == 'vcf')
  {
    sv_regions[[callsets[i, 1]]] <-
      load_vcf(callsets[i, 2], 'DEL', callsets[i, 1])
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

if (tools::file_ext(truth_set_file) == 'vcf')
{
  truth_set <- load_vcf(truth_set_file, 'DEL', 'truth_set')
  ignore.strand.flag <- FALSE
} else if (tools::file_ext(truth_set_file) == 'bedpe')
{
  truth_set <- pairs2breakpointgr(rtracklayer::import(truth_set_file))
  ignore.strand.flag <- TRUE
}

overlaps <- list()
for(c in names(sv_regions))
{
overlaps[[c]] <-
  queryHits(findBreakpointOverlaps(
    truth_set,
    sv_regions[[c]],
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
    restrictMarginToSizeMultiple = 0.5
  ))
}
overlaps[['truth_set']] <- c(1:length(truth_set))
upset(fromList(overlaps), order.by = "freq")
