mean_preprocessing <- function(df, columns, cutoff) {
  
  if (length(columns) == 2) {
    normal_list <- strsplit(df[[columns[1]]], ",", fixed = TRUE)
    tumor_list  <- strsplit(df[[columns[2]]], ",", fixed = TRUE)
    
    safe_mean <- function(x) {
      x <- x[!(x %in% c("", "NA", "NaN"))]
      nums <- as.numeric(x)
      if (length(nums) == 0 || all(is.na(nums))) return(0)
      return(mean(nums, na.rm = TRUE))
    }
    
    mean_PSI_N <- vapply(normal_list, safe_mean, numeric(1))
    mean_PSI_T <- vapply(tumor_list,  safe_mean, numeric(1))
    
    df$mean_PSI_N <- mean_PSI_N
    df$mean_PSI_T <- mean_PSI_T
    
    keep <- !((mean_PSI_N < cutoff | mean_PSI_N > 1 - cutoff) & (mean_PSI_T < cutoff | mean_PSI_T > 1 - cutoff))
    return(df[keep, ])
    
  } else if (length(columns) == 4) {
    normal_1_list <- strsplit(df[[columns[1]]], ",", fixed = TRUE)
    normal_2_list <- strsplit(df[[columns[2]]], ",", fixed = TRUE)
    tumor_1_list  <- strsplit(df[[columns[3]]], ",", fixed = TRUE)
    tumor_2_list  <- strsplit(df[[columns[4]]], ",", fixed = TRUE)
    
    mean_JC_N <- vapply(seq_along(normal_1_list), function(i) {
      mean(as.numeric(normal_1_list[[i]]) + as.numeric(normal_2_list[[i]]))
    }, numeric(1))
    
    mean_JC_T <- vapply(seq_along(tumor_1_list), function(i) {
      mean(as.numeric(tumor_1_list[[i]]) + as.numeric(tumor_2_list[[i]]))
    }, numeric(1))
    
    df$mean_JC_N <- mean_JC_N
    df$mean_JC_T <- mean_JC_T
    
    keep <- (mean_JC_N >= cutoff) & (mean_JC_T >= cutoff)
    return(df[keep, ])
  }
}
