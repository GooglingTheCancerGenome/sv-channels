#!/usr/bin/env Rscript

source('./aux_functions.R')

# create a parser and add command-line arguments
p <-
  arg_parser("Convert sv-channels BEDPE predictions into VCF format.")
p <-
  add_argument(
    p,
    "-i",
    help = paste("a path containing multiple BEDPE files (CV mode"),
    type = "character",
    default = '../genome_wide/results'
  )
p <-
  add_argument(p, "-f", help = "ENCODE blacklist", type = "character",
  default = '../../data/ENCFF001TDO.bed'
  )
p <-
  add_argument(p, "-n", help = "BED file with regions containing Ns", type = "character",
  default = '../../data/reference_N_regions.bed'
  )
p <-
  add_argument(p, "-m", help = "mode: 'split_reads', 'gridss', 'manta', 'delly', 'lumpy'",
  type = "character",
  default="split_reads"
  )
p <-
  add_argument(p, "-o", help = "Output in BEDPE", type = "character",
  default="../genome_wide/results/sv-channels"
  )

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
ref_regions <- argv$n
mode <- argv$m
output_fn <- argv$o

# For each svtype, collect all predictions in BEDPE format and merge the SVs using SVA
sv_regions <- GRangesList()
sv_types <- c('DEL', 'INS', 'INV', 'DUP', 'CTX')
#for (svtype in )

if(file.exists(input_path))
{
  filenames <-
    list.files(
      path = input_path,
      pattern = "correct.bedpe$",
      recursive = TRUE,
      full.names = TRUE
    )
    print(filenames)
}

f <- filenames[1]

bedpe.file <- file.path(input_path, 'results.bedpe')
print(bedpe.file)

cmd <-
  paste(
    "awk '{if($1==$4){print $0\"\t*\t*\t\" $5-$2}else{print $0\"\t*\t*\t\"0}}'",
    f,
    '>',
    bedpe.file
  )
system(cmd)

  print('Loading predictions...')

  svtype <- 'DEL'

  if(file.exists(bedpe.file))
  {

  sv_regions[[svtype]] <-
    load_bedpe(bedpe.file, regions_for_filtering, ref_regions)

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
  }
  print(length(sv_regions[[svtype]]))


# Export to BEDPE
#for (svtype in sv_types)
#{
if(length(sv_regions[[svtype]])>0){

    bp_pairs <- breakpointgr2pairs(sv_regions[[svtype]])
    out_file = paste(output_fn, 'DEL', 'bedpe', sep='.')
    rtracklayer::export(bp_pairs, con = out_file)

}
#}
#concatenate the BEDPE files
#system(paste("cat ", output_fn, ".*.bedpe", " > ", output_fn, ".bedpe", sep=""))


