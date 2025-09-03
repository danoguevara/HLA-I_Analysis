make_skip_iso_peptides <- function(ev_row) {
  gene   <- as.character(ev_row[["GeneID"]])
  txlist <- txs_by_gene[[gene]]
  if (is.null(txlist) || length(txlist) == 0L) return(NULL)
  
  chr    <- as.character(ev_row[["chr"]])
  strand <- as.character(ev_row[["strand"]])
  
  # rMATS starts are 0-based; add +1. Ends are 1-based already.
  up   <- GRanges(chr, IRanges(start = as.integer(ev_row[["upstreamES"]]) + 1L,
                               end   = as.integer(ev_row[["upstreamEE"]])),
                  strand = strand)
  down <- GRanges(chr, IRanges(start = as.integer(ev_row[["downstreamES"]]) + 1L,
                               end   = as.integer(ev_row[["downstreamEE"]])),
                  strand = strand)
  skip <- GRanges(chr, IRanges(start = as.integer(ev_row[["exonStart_0base"]]) + 1L,
                               end   = as.integer(ev_row[["exonEnd"]])),
                  strand = strand)
  
  label <- tolower(as.character(ev_row[["novel"]]))
  out <- list(); k <- 0L
  
  # Only consider transcripts that contain BOTH flanking exons
  tx_with_flanks <- Filter(function(tx) {
    exs <- exons_by_tx[[tx]]
    !is.null(exs) && length(exs) > 0L &&
      length(findOverlaps(exs, up,   ignore.strand = FALSE)) > 0L &&
      length(findOverlaps(exs, down, ignore.strand = FALSE)) > 0L
  }, txlist)
  
  if (identical(label, "canonical")) {
    if (length(tx_with_flanks) == 0L) return(NULL)
    for (tx in tx_with_flanks) {
      exs <- exons_by_tx[[tx]]
      # remove cassette if it overlaps any exon(s)
      ov_skip <- findOverlaps(exs, skip, ignore.strand = FALSE)
      exs2 <- if (length(ov_skip) > 0L) exs[-unique(queryHits(ov_skip))] else exs
      exs2 <- fix_consecutive_rank(exs2, "exon_rank")
      if (length(exs2) == 0L) next
      
      cdna <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38,
                                    GRangesList("tx" = exs2))[[1]]
      
      pep <- NA_character_; frameshift <- NA
      cds <- cds_by_tx[[tx]]
      if (!is.null(cds) && length(cds) > 0L) {
        ov_cds <- findOverlaps(cds, skip, ignore.strand = FALSE)
        cds2 <- if (length(ov_cds) > 0L) cds[-unique(queryHits(ov_cds))] else cds
        cds2 <- fix_consecutive_rank(cds2, "cds_rank")
        
        if (length(cds2) > 0L) {
          cds_seq <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38,
                                           GRangesList("tx" = cds2))[[1]]
          pep <- as.character(translate(cds_seq, if.fuzzy.codon = "solve"))
          frameshift <- (sum(width(cds)) - sum(width(cds2))) %% 3L != 0L
        }
      }
      
      k <- k + 1L
      out[[k]] <- data.frame(
        ID         = as.character(ev_row[["ID"]]),
        GeneID     = gene,
        geneSymbol = as.character(ev_row[["geneSymbol"]]),
        tx         = tx,
        chr        = chr,
        strand     = strand,
        novel      = label,
        cdna       = as.character(cdna),
        peptide    = pep,
        frameshift = frameshift,
        stringsAsFactors = FALSE
      )
    }
  } else if (identical(label, "novel")) {
    if (length(tx_with_flanks) > 0L) {
      for (tx in tx_with_flanks) {
        exs <- exons_by_tx[[tx]]
        ov_skip <- findOverlaps(exs, skip, ignore.strand = FALSE)
        exs2 <- if (length(ov_skip) > 0L) exs[-unique(queryHits(ov_skip))] else exs
        
        exs2 <- fix_consecutive_rank(exs2, "exon_rank")
        if (length(exs2) == 0L) next
        
        cdna <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38,
                                      GRangesList("tx" = exs2))[[1]]
        pep <- longest_orf_peptide(cdna)
        
        k <- k + 1L
        out[[k]] <- data.frame(
          ID         = as.character(ev_row[["ID"]]),
          GeneID     = gene,
          geneSymbol = as.character(ev_row[["geneSymbol"]]),
          tx         = tx,
          chr        = chr,
          strand     = strand,
          novel      = label,
          cdna       = as.character(cdna),
          peptide    = pep,
          frameshift = NA,   # undefined in ORF mode
          stringsAsFactors = FALSE
        )
      }
    } else {
      # Minimal synthetic: just join up+down when no transcript hosts both flanks
      exs2 <- c(up, down)
      exs2 <- fix_consecutive_rank(exs2, "exon_rank")
      cdna <- extractTranscriptSeqs(BSgenome.Hsapiens.UCSC.hg38,
                                    GRangesList("tx" = exs2))[[1]]  
      
      pep <- longest_orf_peptide(cdna)
      out[[1]] <- data.frame(
        ID         = as.character(ev_row[["ID"]]),
        GeneID     = gene,
        geneSymbol = as.character(ev_row[["geneSymbol"]]),
        tx         = NA_character_,
        chr        = chr,
        strand     = strand,
        novel      = label,
        cdna       = as.character(cdna),
        peptide    = pep,
        frameshift = NA,
        stringsAsFactors = FALSE
      )
    }
  } else {
    return(NULL)
  }
  
  do.call(rbind, out)
}