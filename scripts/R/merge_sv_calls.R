#!/usr/bin/env Rscript

options(scipen=999)

suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(argparser))

script.name <-
  basename(sub(".*=", "", commandArgs()[4])) # script name

script.dir <- dirname(commandArgs()[2])
source(file.path(script.dir, 'R/benchmark_functions.R'))

#source('~/Documents/Local_GitHub/sv-channels/scripts/R/benchmark_functions.R')

# create a parser and add command-line arguments
p <-
  arg_parser("Convert sv-channels BEDPE predictions into VCF format.")
p <-
  add_argument(
    p,
    "-i",
    help = paste("a path containing multiple BEDPE files (CV mode"),
    type = "character"
  )
p <-
  add_argument(p, "-f", help = "ENCODE blacklist", type = "character")
p <-
  add_argument(p, "-m", help = "mode: 'split_reads', 'gridss', 'manta', 'delly', 'lumpy'", type = "character")
p <-
  add_argument(p, "-o", help = "Output in BEDPE", type = "character")

# parse the command line arguments
argv <- parse_args(p)
if (is.na(argv$i))
{
  message('Input path missing')
  print(p)
  q(status = 1)
} else if (is.na(argv$o))
{
  message('Output file missing')
  print(p)
  q(status = 1)
}

input_path <- argv$i
regions_for_filtering <- argv$f

# For each svtype, collect all predictions in BEDPE format and merge the SVs using SVA
sv_regions <- list()
sv_types <- c('DEL', 'INS', 'INV', 'DUP', 'TRA')
#for (svtype in )

for (svtype in sv_types)
{
  if(file.exists(file.path(input_path, svtype)))
  {
  # svtype <- 'DEL'
  print(svtype)
  filenames <-
    list.files(
      path = file.path(input_path,
                       svtype,
                       mode),
      pattern = "cnn_predictions.bedpe$",
      recursive = TRUE,
      full.names = TRUE
    )
  
  print(filenames)
  bedpe.file <-
    file.path(input_path, paste(mode, '_', svtype, '.bedpe', sep = ''))
  print(bedpe.file)
  if (file.exists(bedpe.file)) {
    file.remove(bedpe.file)
  }
  
  for (f in filenames)
  {
    cmd <-
      paste(
        "awk '{if($1==$4){print $0\"\t*\t*\t\" $5-$2}else{print $0\"\t*\t*\t\"0}}'",
        f,
        '>>',
        bedpe.file
      )
    system(cmd)
  }
  
  print('Loading predictions...')
  sv_regions[[svtype]] <-
    load_bedpe(bedpe.file, regions_for_filtering)
  
  # rename NA. column of CNN into svLen
  names(mcols(sv_regions[[svtype]]))[names(mcols(sv_regions[[svtype]])) ==
                                       'NA.'] <- 'svLen'
  
  # Merge SVs
  if (svtype != 'INS')
  {
    hits <- findBreakpointOverlaps(
      sv_regions[[svtype]],
      sv_regions[[svtype]],
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
    data <- as.data.frame(hits)
    getfirst <- function(x) {
      x[1]
    }
    get_best_qual <-
      function(x) {
        x[which.max(sv_regions[[svtype]]$QUAL[x])]
      }
    aggdata <-
      aggregate(data$queryHits,
                by = list(data$subjectHits),
                FUN = get_best_qual)
    sv_regions[[svtype]] <- sv_regions[[svtype]][unique(aggdata$x)]
    sv_regions[[svtype]] <-
      sv_regions[[svtype]][sv_regions[[svtype]]$partner %in% names(sv_regions[[svtype]])]
  }
  print(length(sv_regions[[svtype]]))
  }
}

# print(breakpointgr2pairs(sv_regions[[svtype]]))
#callset <- unlist(GRangesList(sv_regions))

# Export to BEDPE
for (t in sv_types)
{
  if(file.exists(file.path(input_path, t)))
  {
  rtracklayer::export(breakpointgr2pairs(sv_regions[[t]]),
                      con = paste("results_", t, ".bedpe", sep = ""))
  }
}
