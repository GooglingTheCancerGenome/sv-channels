require(StructuralVariantAnnotation)

file_paths <-
  c(
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/delly.test_set.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/gridss.test_set.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/last_nanosv.sorted.test_set.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/lumpy.test_set.vcf",
    "/Users/lsantuari/Documents/Data/germline/trio/NA12878/SV/Filtered/Test_set/manta.test_set.vcf",
    "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2011-call-set.test_set.sorted.vcf",
    "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.vcf",
    "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"
  )
names(file_paths) <-
  c('delly',
    'gridss',
    'nanosv',
    'lumpy',
    'manta',
    'mills',
    'deepsv',
    'sv_classify')

gr <- list()

for (sv_caller in names(file_paths)[1:5])
{
  print(paste('considering ', sv_caller, '...', sep = ""))
  sv_callset_vcf <-
    VariantAnnotation::readVcf(as.vector(file_paths[sv_caller]))
  bpgr <- breakpointRanges(sv_callset_vcf)
  begr <- breakendRanges(sv_callset_vcf)
  gr[[sv_caller]] <- sort(c(bpgr, begr))
}

#gr$gridss <- sort(c(bpgr, begr))

p <-
  "/Users/lsantuari/Documents/External_GitHub/sv_benchmark/input.na12878/lumpy-Mills2012-call-set.bedpe"
gr$mills <- bedpe2breakpointgr(p)
seqlevelsStyle(gr$mills) <- "NCBI"
gr$mills <- sort(gr$mills)


p <-
  "/Users/lsantuari/Documents/Processed/Window_pairs_split_reads/060519/predictions/CV_results_DEL_predicted_7.sorted.corrected.bedpe"
gr$deepsv <- bedpe2breakpointgr(p, placeholderName = "deepsv")
gr$deepsv <- sort(gr$deepsv)

p <-
  "/Users/lsantuari/Documents/Data/svclassify/Personalis_1000_Genomes_deduplicated_deletions.bedpe"
gr$sv_classify <-
  bedpe2breakpointgr(p, placeholderName = "sv_classify")
gr$sv_classify <- sort(gr$sv_classify)

for (n in names(gr))
{
  gr[[n]] <- gr[[n]][seqnames(gr[[n]]) %in% c("1", "2", "3")]
}

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
write.table(
  file = paste(
    "~/Documents/Processed/Window_pairs_split_reads/060519/SVAnno_results",
    ".csv",
    sep = ""
  ),
  df,
  row.names = FALSE
)