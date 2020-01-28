#!/usr/bin/env Rscript

suppressPackageStartupMessages(require(StructuralVariantAnnotation))
library(tools)
library(argparser, quietly=TRUE)

script.name<-basename(sub(".*=", "", commandArgs()[4])) # script name

# create a parser and add command-line arguments
p <- arg_parser("Convert VCF output of Manta, DELLY, LUMPY or GRIDSS to BEDPE format.")
p <- add_argument(p, "-i", help="Input in VCF", type="character")
p <- add_argument(p, "-o", help="Output in BEDPE", type="character")
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

if (is.na(argv$o))
{
  argv$o = paste(file_path_sans_ext(argv$i), "bedpe", sep='.')
}

# SV types assigned according to
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
apply_svtype <- function(gr, p_inslen)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "BP",
      ifelse(
        gr$insLen >= abs(gr$svLen) * p_inslen,
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
  gr
}

# t est if there is at least one argument: if not, return an error
#if (length(args) <= 1) {
#  print_usage()
#if (! "o" %in% names(argv))
#{
#  # default output file
#  argv$o = paste(file_path_sans_ext(argv$i), "bedpe", sep='.')
#}

sv_callset_vcf <- VariantAnnotation::readVcf(argv$i)

# Not including breakends (unpaired breakpoints)
# bpgr <- breakpointRanges(sv_callset_vcf)
# breakends are excluded
# begr <- breakendRanges(sv_callset_vcf)
# gr <- sort(c(bpgr, begr))

gr <- breakpointRanges(sv_callset_vcf)

gr <- apply_svtype(gr, p_inslen=argv$p)
# select SVs >= 50 bp. svLen==NA for svtype=='BP'
gr <- gr[abs(gr$svLen) >= argv$l | gr$svtype == 'BP']

bedpe <- breakpointgr2bedpe(gr)

# create a vector with mappings: sourceId -> svtype
svtype_vec <- gr$svtype
names(svtype_vec) <- names(gr)

# typeof(svtype_vec) must be equal to typeof(bedpe_keys)
bedpe_keys <- as.vector(bedpe[,7])
bedpe_svtype <- cbind(bedpe[,1:6], svtype_vec[bedpe_keys])

# check that all SVs with svtype==BP have breakpoints on different chromosomes
if(any(bedpe_svtype[,1]==bedpe_svtype[,4]&bedpe_svtype[,7]=='BP'))
{
  stop("Some SVs with svtype BP contain breakpoints that are both on the same chromosomes")
}

# write output in BEDPE
write.table(file=argv$o, bedpe_svtype, quote=FALSE, row.names=FALSE,
            col.names=FALSE, sep='\t')
