plot_gene_bars <- function(agg_df, annot, gene, free_y = TRUE) {
  df <- dplyr::filter(agg_df, geneSymbol == gene)
  if (nrow(df) == 0) return(NULL)
  
  # Add annotation info
  df <- left_join(df, annot[, c("id", "tumor")],
                  by = c("sample" = "id"))
  
  samp_levels <- annot$id
  df$sample <- factor(df$sample, levels = samp_levels)
  
  ggplot(df, aes(x = sample, y = mean_CPM, fill = tumor)) +
    geom_col(color = "black") + theme_bw() +
    facet_wrap(~ event, nrow = 1, scales = if (free_y) "free_y" else "fixed") +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = c("N" = "lightblue", "T" = "lightcoral"), name = "Sample type") +
    labs(title = paste0(gene), x = "Sample", y = "Counts Per Million") +
    theme(axis.text.x = element_text(size = 12, angle = 45, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(), axis.title.y = element_text(size = 14),
          strip.text = element_text(size = 16), strip.background = element_blank(),
          plot.title = element_text(size = 16, hjust = 0))
}
