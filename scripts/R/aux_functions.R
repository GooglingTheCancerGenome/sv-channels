suppressPackageStartupMessages(require(tools))
suppressPackageStartupMessages(require(StructuralVariantAnnotation))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(argparser))

# this option avoid the scientific notation for numbers (such as 1e2)
options(scipen = 999)

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

make_percent <- function(x) {
    signif(x * 100, digits = 4)
}

getRegionsFromTruthSet <-
function(ts.file) {
    # handle BED(PE) and VCF files
    if (file_ext(ts.file) %in% c('vcf', 'gz')) {
        vcf <- VariantAnnotation::readVcf(ts.file)
        # fix: SVLEN type: CharacterList->IntegerList
        info(vcf)$SVLEN <- IntegerList(info(vcf)$SVLEN)
        return(StructuralVariantAnnotation::breakpointRanges(vcf))
    }

    if (file_ext(ts.file) == 'bed') {
        bedpe.file <- paste0(file_path_sans_ext(ts.file), '.bedpe')
        cmd <-
        paste(
        "awk 'BEGIN {a=0; OFS=\"\t\"} NR>1 {print $1,$2,$2+1,$1,$3,\
        $3+1,\"DEL_\" a,-1,\"+\",\"+\",\"DEL\"; a+=1}'",
        ts.file,
        '>',
        bedpe.file
        )
        system(cmd)
        ts.file <- bedpe.file
    }
    return(pairs2breakpointgr(rtracklayer::import(ts.file)))
}

infer_svtype <- function(gr)
{
    gr$svtype <-
    ifelse(
    seqnames(gr) != seqnames(partner(gr)),
    "CTX",
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
    seqlevelsStyle(bed_regions) <- "NCBI"
    return(bed_regions)
}

load_bedpe <- function(bedpe_file, encode.blacklist, n.regions)
{
    bedpe_gr <- pairs2breakpointgr(rtracklayer::import(bedpe_file))

    if(!is.na(encode.blacklist) && encode.blacklist != ""){
        bedpe_gr <- filter_regions('ENCODE blacklist', bedpe_gr, load_bed(encode.blacklist), mode = 'remove')
    } 
    if(!is.na(n.regions) && n.regions != ""){
        bedpe_gr <- filter_regions('regions with Ns', bedpe_gr, load_bed(n.regions), mode = 'remove')
    }
    return(bedpe_gr)
}

load_vcf <- function(vcf_file, svtype, caller, encode.blacklist, n.regions)
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
    vcf_gr <- filter_regions('ENCODE blacklist', vcf_gr, load_bed(encode.blacklist), mode = 'remove')
    vcf_gr <- filter_regions('regions with Ns', vcf_gr, load_bed(n.regions), mode = 'remove')

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

count_sv_lines <- function(filename){
    length(grep('^#', readLines(filename), perl = TRUE, invert = TRUE)) > 0
}

############
