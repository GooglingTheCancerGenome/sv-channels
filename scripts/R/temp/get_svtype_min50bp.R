# run as Rscript --vanilla get_svype.R -i input.vcf -o output.vcf

# infer SVTYPE and filter SVs <= 50 bp

library("optparse")
suppressPackageStartupMessages(require(StructuralVariantAnnotation))

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

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input vcf file", metavar="character"),
  make_option(c("-o", "--output"), type="character", default='tmp.vcf', 
              help="output vcf file", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

input_vcf <- opt$input
output_vcf <- opt$output

sv_callset_vcf <-
  VariantAnnotation::readVcf(input_vcf)

bpgr <- breakpointRanges(sv_callset_vcf)
begr <- breakendRanges(sv_callset_vcf)
gr <- sort(c(bpgr, begr))
gr <- apply_svtype(gr)

svtype_vec <- gr$svtype
names(svtype_vec) <- gr$sourceId

#table(info(sv_callset_vcf)$SVTYPE)
info(sv_callset_vcf)$SVTYPE <- as.vector(svtype_vec[rownames(info(sv_callset_vcf))])
#table(info(sv_callset_vcf)$SVTYPE)

gt50 <- gr$sourceId[abs(gr$svLen) >= 50]
sv_callset_vcf <- sv_callset_vcf[which(rownames(info(sv_callset_vcf))%in%gt50),]

summary(abs(unlist(info(sv_callset_vcf)$SVLEN)))

VariantAnnotation::writeVcf(sv_callset_vcf, output_vcf)
