get_recurrent_events <- function(counts, res, annot, event) {
  df <- counts
  df$event <- event
  
  #if (event == "RI") {
  #  df <- df[res$IncLevelDifference < 0, ] 
  #} else {
  #  df <- df[res$IncLevelDifference > 0, ] 
  #}
  
  normal_ids <- annot$id[annot$tumor == "N"]
  tumor_ids <- annot$id[annot$tumor == "T"]
  
  normal_max <- do.call(pmax, c(df[normal_ids], na.rm = TRUE))
  tumor_min  <- do.call(pmin, c(df[tumor_ids], na.rm = TRUE))
  
  boolean_vector <- tumor_min > normal_max
  
  df <- df[boolean_vector, ]
  
  return(df)
}