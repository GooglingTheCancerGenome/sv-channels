hits <- findBreakpointOverlaps(
  sv_regions[['manta']],
  sv_regions[['sv-channels']],
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = ignore.strand.flag,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)
length(sv_regions[['manta']][-queryHits(hits),])


hits.m.t <- findBreakpointOverlaps(
  sv_regions[['manta']],
  truth_set,
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = ignore.strand.flag,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)
hits.s.t <- findBreakpointOverlaps(
  sv_regions[['sv-channels']],
  truth_set,
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = ignore.strand.flag,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)

m.no.s <- subjectHits(hits.m.t)[-which(subjectHits(hits.m.t)%in%subjectHits(hits.s.t))]
sv_regions[['manta']][m.no.s]


m.in.s <- findBreakpointOverlaps(
  sv_regions[['manta']][unique(queryHits(hits))],
  truth_set,
  # read pair based callers make imprecise calls.
  # A margin around the call position is required when matching with the truth set
  maxgap = 100,
  # Since we added a maxgap, we also need to restrict the mismatch between the
  # size of the events. We don't want to match a 100bp deletion with a
  # 5bp duplication. This will happen if we have a 100bp margin but don't also
  # require an approximate size match as well
  sizemargin = 0.25,
  ignore.strand = ignore.strand.flag,
  # We also don't want to match a 20bp deletion with a 20bp deletion 80bp away
  # by restricting the margin based on the size of the event, we can make sure
  # that simple events actually do overlap
  restrictMarginToSizeMultiple = 0.5
)
length(m.in.s)
length(hits.m.t)
length(hits.s.t)
