get_skipping_counts <- function(df, annot, preffix = "SJC") {
  
  if (preffix == "IJC") {
    SJC_counts <- df[, c("IJC_SAMPLE_1", "IJC_SAMPLE_2")]
    SJC_counts <- SJC_counts %>%
      separate(IJC_SAMPLE_1, into = paste0("IJC_SAMPLE_1_", 1:5), sep = ",", fill = "right", convert = TRUE) %>%
      separate(IJC_SAMPLE_2, into = paste0("IJC_SAMPLE_2_", 1:15), sep = ",", fill = "right", convert = TRUE)
    
  } else {
    SJC_counts <- df[, c("SJC_SAMPLE_1", "SJC_SAMPLE_2")]
    SJC_counts <- SJC_counts %>%
      separate(SJC_SAMPLE_1, into = paste0("SJC_SAMPLE_1_", 1:5), sep = ",", fill = "right", convert = TRUE) %>%
      separate(SJC_SAMPLE_2, into = paste0("SJC_SAMPLE_2_", 1:15), sep = ",", fill = "right", convert = TRUE)
  }
  
  colnames(SJC_counts) <- annot$id
  
  SJC_CPM <- SJC_counts
  for (sample in colnames(SJC_CPM)[1:20]) {
    SJC_CPM[[sample]] <- SJC_CPM[[sample]]*(10**6)/annot$lib_size[annot$id == sample]
  }
  
  # Counts distribution
  SJC_CPM$Mean.N <- rowMeans(SJC_CPM[, annot$id[annot$tumor == "N"]])
  SJC_CPM$Mean.T <- rowMeans(SJC_CPM[, annot$id[annot$tumor == "T"]])
  SJC_CPM$Mean <- rowMeans(SJC_CPM[, annot$id])
  
  SJC_CPM <- cbind(SJC_CPM, df[, c("ID", "geneSymbol", "GeneID")])
  
  return(SJC_CPM)
}