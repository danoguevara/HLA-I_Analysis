longest_orf_peptide <- function(dna) {
  best <- ""
  for (f in 0:2) {
    aa <- as.character(translate(subseq(dna, start = 1 + f), if.fuzzy.codon = "solve"))
    segs <- strsplit(aa, "\\*")[[1]]  # split at stops
    for (seg in segs) {
      mpos <- regexpr("M", seg, fixed = TRUE)[1]
      if (mpos > 0) {
        pep <- substr(seg, mpos, nchar(seg))
        if (nchar(pep) > nchar(best)) best <- pep
      }
    }
  }
  if (nchar(best) == 0) NA_character_ else best
}