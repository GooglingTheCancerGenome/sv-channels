# to run: Rscript --vanilla vcf2bedpe.R -i [svcaller].vcf -o [svcaller].bedpe
# require StructuralVariantAnnotation (Bioconductor) and argparser (CRAN)

suppressPackageStartupMessages(require(StructuralVariantAnnotation))
library(tools)
library(argparser, quietly=TRUE)

# Create a parser
p <- arg_parser("Convert GRIDSS/manta/LUMPY/DELLY VCF output to BEDPE")

# Add command line arguments
p <- add_argument(p, "-i", help="VCF input file", type="character")
p <- add_argument(p, "-o", help="BEDPE output file", type="character")
p <- add_argument(p, "-p", help="SVTYPE=INS if insertion length greater or equal to SV_length * p", type="numeric", default=0.7)
p <- add_argument(p, "-l", help="minimum SV length to consider", type="int", default=50)

# Parse the command line arguments
argv <- parse_args(p)

#SV type assignment based on
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

# test if there is at least one argument: if not, return an error
if (length(argv)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (! "p" %in% names(argv)) {
  # default output file
  argv$o = paste(file_path_sans_ext(argv$i), "bedpe", sep='.')
}

sv_callset_vcf <-
  VariantAnnotation::readVcf(argv$i)

# Not including breakends (unpaired breakpoints)
# bpgr <- breakpointRanges(sv_callset_vcf)
# breakends are excluded
# begr <- breakendRanges(sv_callset_vcf)
# gr <- sort(c(bpgr, begr))

gr <- breakpointRanges(sv_callset_vcf)

gr <- apply_svtype(gr, p_inslen=argv$p)
# Select SVs >= 50 bp. svLen==NA for svtype=='BP'
gr <- gr[abs(gr$svLen) >= argv$l | gr$svtype == 'BP']

bedpe <- breakpointgr2bedpe(gr)

# create vector with mapping sourceId -> svtype
svtype_vec <- gr$svtype
names(svtype_vec) <- names(gr)

# typeof(svtype_vec) must be equal to typeof(bedpe_keys)
bedpe_keys <- as.vector(bedpe[,7])
bedpe_svtype <- cbind(bedpe[,1:6], svtype_vec[bedpe_keys])

#check that all SVs with svtype==BP have breakpoints on different chromosomes
if(any(bedpe_svtype[,1]==bedpe_svtype[,4]&bedpe_svtype[,7]=='BP'))
{
  stop("Some SVs with svtype BP contain breakpoints that are both on the same chromosomes")
}

# write output
write.table(file=argv$o, bedpe_svtype, quote=FALSE, 
            row.names = FALSE, col.names = FALSE, sep='\t')
