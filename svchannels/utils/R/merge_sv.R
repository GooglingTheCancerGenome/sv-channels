rm(list = ls())
suppressPackageStartupMessages(require(argparser))
suppressPackageStartupMessages(require(ggplot2))
suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))


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

load_bed <- function(bed_file)
{
  bed_regions <- rtracklayer::import(bed_file)
  # set NCBI seqlevels
  seqlevelsStyle(bed_regions) <- "UCSC"
  return(bed_regions)
}

load_vcf <- function(vcf_file, svtype, caller, filter.regions)
{
  # Load VCF file
  vcf_gr <-
    VariantAnnotation::readVcf(vcf_file)
  
  # Keep only SVs that passed the filtering (PASS or .)
  vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("PASS", ".")]
  
  if (caller == 'lumpy')
  {
    # Read evidence support as a proxy for QUAL
    fixed(vcf_gr)$QUAL <- unlist(info(vcf_gr)$SU)
    
  } else if (caller == 'delly')
  {
    # Split-read support plus Paired-end read support as a proxy for QUAL
    sr_support <- info(vcf_gr)$SR
    sr_support[is.na(vcf_gr)] <- 0
    fixed(vcf_gr)$QUAL <- sr_support + info(vcf_gr)$PE
  }
  
  vcf_gr <- StructuralVariantAnnotation::breakpointRanges(vcf_gr)
  vcf_gr <- infer_svtype(vcf_gr)
  
  # Select only one SV type
  vcf_gr <- vcf_gr[which(vcf_gr$svtype == svtype)]
  
  message(paste('Loading', length(vcf_gr), svtype, 'calls for', caller, sep=" "))
  
  # Select SVs >= 50 bp
  if (! svtype %in% c('CTX', 'INS'))
  {
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]
  }
  
  # Filter regions
  vcf_gr <- filter_regions('ENCODE blacklist', vcf_gr, load_bed(filter.regions), mode = 'remove')
  
  # select only chromosomes 1 to 22
  vcf_gr <-
    vcf_gr[seqnames(vcf_gr) %in% paste('chr', c(1:22), sep = '')]
  vcf_gr <- vcf_gr[vcf_gr$partner %in% names(vcf_gr)]
  
  return(vcf_gr)
}

filter_regions <-
  function(filter.name, regions_to_filter, ref_regions, mode = 'remove')
  {
    print(length(regions_to_filter))
    if (mode == 'keep')
    {
      result <-
        regions_to_filter[overlapsAny(regions_to_filter, ref_regions) &
                            overlapsAny(partner(regions_to_filter), ref_regions),]
    } else if (mode == 'remove') {
      result <- regions_to_filter[! (
        overlapsAny(regions_to_filter, ref_regions) |
          overlapsAny(partner(regions_to_filter), ref_regions)
      ),]
    }
    message(paste(length(result), 'calls after filtering for', filter.name, sep=" "))
    return(result)
  }

p <-
  arg_parser("Merge Manta and GRIDSS2 SVs in VCF format. Output in BEDPE format.")
p <-
  add_argument(p, "-m", help = "Manta input in VCF", type = "character",
               default = '/Users/lsantuari/Documents/GitHub/sv-channels/notebooks/SV/manta.vcf')
p <-
  add_argument(p, "-g", help = "GRIDSS2 input in VCF", type = "character",
               default = '/Users/lsantuari/Documents/GitHub/sv-channels/notebooks/SV/HG00420.vcf.gz')

p <-
  add_argument(p, "-filter", help = "Filter regions in BED format", type = "character",
               default = '/Users/lsantuari/Documents/Projects/GTCG/sv-channels/sv-channels_manuscript/1KG_trios/reference/exclude.N-and-encodebl.bed')
p <-
  add_argument(p, "-o", help = "Output in BEDPE", type = "character",
               default = '/Users/lsantuari/Documents/GitHub/sv-channels/notebooks/SV/Manta_union_GRIDSS2.bedpe')

argv <- parse_args(p)
manta_vcf <- argv$m
gridss_vcf <- argv$g
filter <- argv$filter

sv_regions <- list()
sv_regions['manta'] <-
  load_vcf(manta_vcf, 'DEL', 'manta', filter)
sv_regions['gridss'] <-
  load_vcf(gridss_vcf, 'DEL', 'gridss', filter)

hits <- findBreakpointOverlaps(
  sv_regions[['manta']],
  sv_regions[['gridss']],
  # using a smaller margin for insertions, insertion location should be precise
  maxgap = 100,
  # sizemargin cannot be used for insertions
  # sizemargin = 0.25,
  ignore.strand = TRUE,
  restrictMarginToSizeMultiple = 0.5
)

all_sv <- c(sv_regions[['manta']],
  sv_regions[['gridss']][-subjectHits(hits)])
bedpe <- StructuralVariantAnnotation::breakpointgr2bedpe(all_sv)
svtype_vec <- all_sv$svtype
names(svtype_vec) <- names(all_sv)

bedpe_keys <- as.vector(bedpe[, 7])
bedpe_svtype <- cbind(bedpe[, 1:6], SVTYPE = svtype_vec[bedpe_keys], name = bedpe$name)
bedpe_svtype <- bedpe_svtype[,-which(names(bedpe_svtype)=='name')]

print('SVTYPEs in BEDPE output:')
print(table(bedpe_svtype$SVTYPE))

#only write deletions
bedpe_svtype <- bedpe_svtype[which(bedpe_svtype$SVTYPE=='DEL'),]

# write output in BEDPE
write.table(
  file = argv$o,
  bedpe_svtype,
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE,
  sep = '\t'
)
