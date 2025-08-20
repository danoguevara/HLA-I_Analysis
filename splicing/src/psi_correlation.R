psi_correlation <- function(df1, df2, columns = c("ID", "geneSymbol", "IncLevelDifference"),
                            title, p_size = 1, p_shape = 20, label = FALSE, cutoff = NULL)
{
  # Integrate data frames
  merged_df <- merge(df1[, columns], df2[, columns],
                     by = columns[1:2], suffixes = c("_JC", "_JCEC"))
  
  # Annotation position
  x_rng <- range(merged_df$IncLevelDifference_JC, na.rm = TRUE)
  y_rng <- range(merged_df$IncLevelDifference_JCEC, na.rm = TRUE)
  x_pos <- x_rng[2] - 0.15 * diff(x_rng)   # 5% in from left
  y_pos <- y_rng[1] + 0.08 * diff(y_rng)   # 5% down from top
  
  # Pearson correlation
  cor_test <- cor.test(merged_df$IncLevelDifference_JC, merged_df$IncLevelDifference_JCEC, method = "pearson")
  r_squared <- round(cor_test$estimate^2, 3)
  
  # Plotting
  plot <- ggplot(merged_df, aes(x = IncLevelDifference_JC, y = IncLevelDifference_JCEC)) +
    geom_point(alpha = 0.7, size = p_size, shape = p_shape) + theme_bw() +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    annotate(geom = "point", x = 0, y = 0, color = "red", fill = "white", size = p_size+1, shape = 21) +
    annotate("text", x = x_pos, y = y_pos, label = paste0("R² = ", r_squared), size = 4, hjust = 0) +
    labs(title = title,
         x = expression(Delta*" PSI"*" - Junction Count"),
         y = expression(Delta*" PSI"*" - Junction & Exon Count")) +
    scale_x_continuous(breaks = seq(-0.8, 0.8, 0.4)) +
    scale_y_continuous(breaks = seq(-0.8, 0.8, 0.4)) +
    coord_cartesian(xlim = c(min(merged_df$IncLevelDifference_JC) - 0.1,
                             max(merged_df$IncLevelDifference_JC)   + 0.1),
                    ylim = c(min(merged_df$IncLevelDifference_JCEC) - 0.1,
                             max(merged_df$IncLevelDifference_JCEC) + 0.1))
  
  if (label) {
    label_df <- subset(merged_df, abs(IncLevelDifference_JC) >= cutoff | abs(IncLevelDifference_JCEC) >= cutoff)
    label_df$nudge_y <- ifelse(label_df$IncLevelDifference_JCEC > label_df$IncLevelDifference_JC, 0.08, -0.08)
    
    plot <- plot + geom_label_repel(data = label_df, aes(label = geneSymbol), size = 3,
                                    nudge_y = label_df$nudge_y, segment.alpha = 0.5, segment.size = 0.3,
                                    box.padding = unit(0.4, "lines"), point.padding = unit(0.3, "lines"),
                                    max.overlaps = Inf, show.legend = FALSE)
  }
  
  return(plot)
}