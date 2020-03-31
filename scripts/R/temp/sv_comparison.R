

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#
# # The following initializes usage of Bioc devel
# BiocManager::install(version='devel')
#
# BiocManager::install("StructuralVariantAnnotation")

require(vcfR)
require(StructuralVariantAnnotation)


# file_paths <-
#   c(
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/delly.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/gridss.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/last_nanosv.sorted.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/lumpy.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/manta.test_set.vcf",
#     "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2011-call-set.test_set.sorted.vcf",
#     "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.vcf",
#     "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"
#   )

file_paths <-
  c(
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/delly.flt.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/gridss.flt.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/last_nanosv.sorted.flt.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/lumpy.flt.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/manta.flt.vcf",
    "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2011-call-set.bed",
    "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.vcf",
    "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"#,
    #"/Users/lsantuari/Documents/Data/germline/N"
  )

# file_paths <-
#   c(
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Unfiltered/delly.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Unfiltered/gridss.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Unfiltered/last_nanosv.sorted.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Unfiltered/lumpy.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Unfiltered/manta.vcf",
#     "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2011-call-set.bed",
#     "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.vcf",
#     "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"
#   )

# file_paths <-
#   c(
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/delly.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/gridss.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/last_nanosv.sorted.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/lumpy.test_set.vcf",
#     "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/manta.test_set.vcf",
#     "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2011-call-set.test_set.sorted.vcf",
#     "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.vcf",
#     "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"#,
#     #"/Users/lsantuari/Documents/Data/germline/N"
#   )

names(file_paths) <-
  c('delly',
    'gridss',
    'nanosv',
    'lumpy',
    'manta',
    'mills',
    'deepsv',
    'sv_classify')

gridss.vcf <- read.vcfR(file_paths['gridss'])

gr <- list()

for (sv_caller in names(file_paths)[1:5])
{
  #sv_caller <- 'gridss'
  print(paste('considering ', sv_caller, '...', sep = ""))
  sv_callset_vcf <-
    VariantAnnotation::readVcf(as.vector(file_paths[sv_caller]))
  bpgr <- breakpointRanges(sv_callset_vcf)
  begr <- breakendRanges(sv_callset_vcf)
  gr[[sv_caller]] <- sort(c(bpgr, begr))
}

#gr$gridss <- sort(c(bpgr, begr))

bedpe.file <-
  "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe"
bedpe.pairs <- rtracklayer::import(bedpe.file)
gr$mills <- pairs2breakpointgr(bedpe.pairs)
seqlevelsStyle(gr$mills) <- "NCBI"
gr$mills <- sort(gr$mills)


# Load the NA24385 SVs
vcf_file <- 'Documents/Data/germline/NA24385/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz'
sv_callset_vcf <-
  VariantAnnotation::readVcf(vcf_file)
bpgr <- breakpointRanges(sv_callset_vcf)
begr <- breakendRanges(sv_callset_vcf)
gr[['NA24385_SV']] <- sort(c(bpgr, begr))
gr[['NA24385_SV']] <- gr[['NA24385_SV']][which(gr[['NA24385_SV']]$svtype=="DEL")]

# Load the CHM1_CHM13 SVs
vcf_file <- 'Documents/Data/germline/CHM/Huddleston2016/structural_variants/CHM1_CHM13_pseudodiploid_SVs.vcf.gz'
sv_callset_vcf <-
  VariantAnnotation::readVcf(vcf_file)
bpgr <- breakpointRanges(sv_callset_vcf)
begr <- breakendRanges(sv_callset_vcf)
gr[['CHM1_CHM13_SV']] <- sort(c(bpgr, begr))
gr[['CHM1_CHM13_SV']] <- gr[['CHM1_CHM13_SV']][which(gr[['CHM1_CHM13_SV']]$svtype=="DEL")]
seqlevelsStyle(gr$CHM1_CHM13_SV) <- "NCBI"

# p <-
#  "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.corrected.bedpe"
sample_name <- "CHM1_CHM13"
bedpe.file <-
  paste(
    "/Users/lsantuari/Documents/Processed/channel_maker_output/",
    sample_name,
    "/candidate_positions.bedpe.gz",
    sep = ""
  )
bedpe.pairs <- rtracklayer::import(bedpe.file)
gr$deepsv <-
  pairs2breakpointgr(bedpe.pairs, placeholderName = "deepsv")
gr$deepsv <- sort(gr$deepsv)


vcf_file <- 'Documents/Data/germline/NA24385/SV/Filtered/gridss.vcf'
sv_callset_vcf <-
  VariantAnnotation::readVcf(vcf_file)
bpgr <- breakpointRanges(sv_callset_vcf)
begr <- breakendRanges(sv_callset_vcf)
gr[['gridss']] <- sort(c(bpgr, begr))



bedpe.file <-
  "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"
bedpe.pairs <- rtracklayer::import(bedpe.file)
gr$sv_classify <-
  pairs2breakpointgr(bedpe.pairs, placeholderName = "sv_classify")
gr$sv_classify <- sort(gr$sv_classify)

# for (n in names(gr))
# {
#   gr[[n]] <- gr[[n]][seqnames(gr[[n]]) %in% c("1", "2", "3")]
# }

df <- data.frame()
for (gt_dataset in c('mills', 'sv_classify'))
{
  ground_truth <- gr[[gt_dataset]]
  for (filter in c(TRUE, FALSE))
  {
    filtered <- if (filter)
      "_filtered"
    else
      ""
    
    for (sv_caller in names(gr)[names(gr) %in% c('delly', 'gridss', 'nanosv', 'lumpy', 'manta', 'deepsv')])
    {
      print(sv_caller)
      # sv_caller <- 'gridss'
      sv_callset <- gr[[sv_caller]]
      if (sv_caller != 'deepsv' & filter == TRUE) {
        sv_callset <- sv_callset[sv_callset$FILTER %in% c('.', 'PASS')]
      }
      overlaps <-
        countBreakpointOverlaps(sv_callset,
                                ground_truth,
                                maxgap = 200,
                                ignore.strand = TRUE)
      precision <-
        round(length(overlaps[overlaps > 0]) / length(sv_callset), digits = 2)
      recall <-
        round(length(overlaps[overlaps > 0]) / length(ground_truth), digits = 2)
      f1_score <-
        round(2 * (precision * recall) / (precision + recall), digits = 2)
      df <-
        rbind(
          df,
          data.frame(
            sv_caller = sv_caller,
            filtered = if (filter)
              "PASS"
            else
              "all",
            ground_truth = gt_dataset,
            precision = precision,
            recall = recall,
            f1_score = f1_score,
            TP = length(overlaps[overlaps > 0]),
            FP = length(sv_callset) - length(overlaps[overlaps > 0]),
            FN = length(ground_truth) - length(overlaps[overlaps > 0])
          )
        )
    }
  }
}
df <- df[order(df$f1_score, decreasing = TRUE),]

# write.table(
#   file = paste(
#     "~/Documents/Processed/Window_pairs_split_reads/060519/SVAnno_results",
#     ".csv",
#     sep = ""
#   ),
#   df,
#   row.names = FALSE
# )

df <- data.frame()
sv_callset_ref <-  gr$gridss
gap <- 200

# for (n in names(gr)[names(gr)!='deepsv'&names(gr)!='gridss'])
for (n in names(gr)[names(gr) != 'gridss'])
{
  print(n)
  sv_callset <-  gr[[n]]
  overlaps <-
    countBreakpointOverlaps(sv_callset,
                            sv_callset_ref,
                            maxgap = gap,
                            ignore.strand = TRUE)
  
  overlaps_mills <-
    countBreakpointOverlaps(gr$mills,
                            sv_callset[overlaps > 0],
                            maxgap = gap,
                            ignore.strand = TRUE)
  
  overlaps_sv_classify <-
    countBreakpointOverlaps(gr$sv_classify,
                            sv_callset[overlaps > 0],
                            maxgap = gap,
                            ignore.strand = TRUE)
  
  df <-
    rbind(
      df,
      data.frame(
        svcaller = n,
        overlap = round(length(overlaps[overlaps > 0]) / length(overlaps), digits = 2),
        sv_callset_size = length(overlaps),
        number_of_overlaps = length(overlaps[overlaps > 0]),
        overlap_mills = round(length(overlaps_mills[overlaps_mills > 0]) / length(overlaps_mills), digits = 2),
        sv_callset_size_mills = length(overlaps_mills),
        number_of_overlaps_mills = length(overlaps_mills[overlaps_mills > 0]),
        overlap_sv_classify = round(
          length(overlaps_sv_classify[overlaps_sv_classify > 0]) / length(overlaps_sv_classify),
          digits = 2
        ),
        sv_callset_size_sv_classify = length(overlaps_sv_classify),
        number_of_overlaps_sv_classify = length(overlaps_sv_classify[overlaps_sv_classify > 0])
      )
    )
}
df

write.table(
  file = paste(
    "~/Documents/Processed/Window_pairs_split_reads/060519/Candidate_positions_overlaps_GRIDSS",
    ".csv",
    sep = ""
  ),
  df,
  row.names = FALSE
)

overlaps <-
  countBreakpointOverlaps(gr$sv_classify,
                          gr$gridss,
                          maxgap = 100,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)
overlaps <-
  countBreakpointOverlaps(gr$sv_classify,
                          gr$deepsv,
                          maxgap = 100,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)



overlaps <-
  countBreakpointOverlaps(gr$sv_classify,
                          gr$deepsv,
                          maxgap = 200,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)

overlaps <-
  countBreakpointOverlaps(gr$sv_classify,
                          gr$gridss,
                          maxgap = 200,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)


gr$deepsv[seqnames(gr$deepsv)%in%c(1:22,'X')]
overlaps <-
  countBreakpointOverlaps(gr$CHM1_CHM13_SV,
                          gr$deepsv,
                          maxgap = 200,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)


rtracklayer::export(
  gr$sv_classify[which(overlaps == 0)],
  con = paste(
    "/Users/lsantuari/Documents/Processed/channel_maker_output/",
    sample_name,
    "sv_classify.not_in_DeepSV_cpos.bedpe",
    sep = ""
  )
)

ground_truth <-
  gr$sv_classify[which(seqnames(gr$sv_classify) %in% c('1', '2', '3'))]

gridss_svcallset <- gr$gridss[which(gr$gridss$svtype=="DEL")]
gridss_svcallset <-
  gridss_svcallset[which(seqnames(gridss_svcallset) %in% c('1', '2', '3'))]


overlaps <-
  countBreakpointOverlaps(ground_truth,
                          gridss_svcallset,
                          maxgap = 200,
                          ignore.strand = TRUE)
length(overlaps[overlaps > 0]) / length(overlaps)
length(overlaps[overlaps > 0]) / length(gridss_svcallset)

gr$gridss$svtype <-
  ifelse(
    seqnames(gr$gridss) != seqnames(partner(gr$gridss)),
    "BP",
    ifelse(
      gr$gridss$insLen >= abs(gr$gridss$svLen) * 0.7,
      "INS",
      ifelse(
        strand(gr$gridss) == strand(partner(gr$gridss)),
        "INV",
        ifelse(xor(
          start(gr$gridss) < start(partner(gr$gridss)), strand(gr$gridss) == "-"
        ), "DEL",
        "DUP")
      )
    )
  )

gr$gridss <- gr$gridss[gr$gridss$svtype=="DEL"]

gridss_svcallset["gridss69_18809h"]
length(gridss_svcallset)

length(names(gr$gridss))
write.table(file="Documents/Data/germline/trio/NA12878/SV/Filtered/gridss.sym.sv_benchmark.ids.txt", names(gr$gridss), row.names = F, col.names = F, quote=F)

unique(seqnames(gr$gridss))

gridss.vcf.del <- gridss.vcf[which(getID(gridss.vcf)%in%names(gr$gridss))]
head(getINFO(gridss.vcf.del))

mine_id <- read.table(file="Documents/Data/germline/trio/NA12878/SV/Filtered/gridss.sym.mine.ids.sorted.txt")[,1]
mine_id <- c(paste(mine_id,'h', sep=''), paste(mine_id,'o', sep=''))

gridss_svcallset <- gr$gridss[which(gr$gridss$svtype=="DEL")]
length(gridss_svcallset[which(names(gridss_svcallset)%in%mine_id)])/length(gridss_svcallset)

length(which(gr$gridss[which(names(gr$gridss)%in%mine_id)]$svtype=="INS"))
