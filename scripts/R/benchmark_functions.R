
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))

# SV type inference
infer_svtype <- function(gr)
{
  gr$svtype <-
    ifelse(
      seqnames(gr) != seqnames(partner(gr)),
      "TRA", # Using TRA instead of ITX or BP
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

# load a BED file and returns a GRanges object
load_bed <- function(bed_file)
{
    bed_regions <- rtracklayer::import(bed_file)
    # set NCBI seqlevels
    seqlevelsStyle(bed_regions) <- "NCBI"
    return(bed_regions)
}

# load a BEDPE file and returns a Pairs object
load_bedpe <- function(bedpe_file, filter_regions)
{
    bedpe_gr <- pairs2breakpointgr(rtracklayer::import(bedpe_file))
    bedpe_gr <- filter_regions(bedpe_gr, load_bed(filter_regions), mode = 'remove')
    return(bedpe_gr)
}

# Filter regions using a BED file either for inclusion (keep) or exclusion (remove)
filter_regions <- function(regions_to_filter, ref_regions, mode='remove')
{
    print(length(regions_to_filter))
    if (mode == 'keep')
    {
        result <- regions_to_filter[overlapsAny(regions_to_filter, ref_regions) &
        overlapsAny(partner(regions_to_filter), ref_regions),]
    } else if (mode == 'remove') {
        result <- regions_to_filter[! (
        overlapsAny(regions_to_filter, ref_regions) |
        overlapsAny(partner(regions_to_filter), ref_regions)
        ),]
    }
    print(length(result))
    return(result)
}

# load a VCF file and returns a Pairs object
load_vcf <- function(vcf_file, svtype, caller, exclude_regions, include_regions)
{
    # Load VCF file
    vcf_gr <-
    VariantAnnotation::readVcf(vcf_file)

    # set NCBI seqlevels
    seqlevelsStyle(vcf_gr) <- 'NCBI'

    if (caller == 'survivor')
    {
        # SURVIVOR simSV assigns LowQual to all artificial SVs
        vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("LowQual")]

        # handle SURVIVOR simSV SVTYPE <TRA> as TIGRA <CTX>
        if ("TRA" %in% info(vcf_gr)$SVTYPE & ! ("CT" %in% names(info(vcf_gr))))
        {
            info(vcf_gr)$SVTYPE[info(vcf_gr)$SVTYPE == "TRA"] <- "CTX"
        }
    }else {
        # Keep only SVs that passed the filtering (PASS or .)
        vcf_gr <- vcf_gr[rowRanges(vcf_gr)$FILTER %in% c("PASS", ".")]
    }

    if (caller == 'lumpy')
    {
        # Read evidence support as a proxy for QUAL
        fixed(vcf_gr)$QUAL <- unlist(info(vcf_gr)$SU)
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
    vcf_gr <- vcf_gr[which(vcf_gr$svtype == svtype)]

    # Select SVs >= 50 bp
    vcf_gr <- vcf_gr[abs(vcf_gr$svLen) >= 50]

    #Filter regions
    if(file.exists(exclude_regions))
    {
        vcf_gr <- filter_regions(vcf_gr, load_bed(exclude_regions), mode = 'remove')
    }
    if(file.exists(include_regions))
    {
        vcf_gr <- filter_regions(vcf_gr, load_bed(include_regions), mode = 'keep')
    }

    return(vcf_gr)
}