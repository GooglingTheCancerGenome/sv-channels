# to run: Rscript --vanilla vcf2bedpe.R [svcaller].vcf [svcaller].bedpe

suppressPackageStartupMessages(require(StructuralVariantAnnotation))
library(tools)

#SV type assignment based on
# https://github.com/PapenfussLab/gridss/blob/7b1fedfed32af9e03ed5c6863d368a821a4c699f/example/simple-event-annotation.R#L9
apply_svtype <- function(gr)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "BP",
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
  gr
}

args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = paste(file_path_sans_ext(args[1]), "bedpe", sep='.')
}

sv_callset_vcf <-
  VariantAnnotation::readVcf(args[1])

# bpgr <- breakpointRanges(sv_callset_vcf)
# breakends are excluded
# begr <- breakendRanges(sv_callset_vcf)
# gr <- sort(c(bpgr, begr))

gr <- breakpointRanges(sv_callset_vcf)
gr <- apply_svtype(gr)
# Select SVs >= 50 bp
gr <- gr[abs(gr$svLen) >= 50]

bedpe <- breakpointgr2bedpe(gr)

write.table(file=args[2], bedpe, quote=FALSE, row.names = FALSE, col.names = FALSE, sep='\t')