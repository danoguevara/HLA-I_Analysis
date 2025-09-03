fix_consecutive_rank <- function(gr, rank_col = c("exon_rank","cds_rank")) {
  if (length(gr) == 0) return(gr)
  rank_col <- match.arg(rank_col)
  if (rank_col %in% colnames(mcols(gr))) {
    # order by current rank, then renumber 1..n
    ord <- order(mcols(gr)[[rank_col]])
    gr  <- gr[ord]
    mcols(gr)[[rank_col]] <- seq_len(length(gr))
  } else {
    # no rank present: sort 5'→3' respecting strand for determinism
    s <- as.character(unique(strand(gr)))
    if (length(s) == 1 && s %in% c("+","-")) {
      ord <- if (s == "+") order(start(gr), end(gr)) else order(start(gr), end(gr), decreasing = TRUE)
      gr  <- gr[ord]
    }
  }
  gr
}