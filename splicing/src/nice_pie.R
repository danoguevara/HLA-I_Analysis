nice_pie <- function(df, selected_cols = c("typeCol", "freqCol"), colors, order = FALSE, levels = NULL) {
  # Manage data
  pie.df <- df[, selected_cols]
  pie.df$percent <- round(pie.df[, 2]*100/sum(pie.df[, 2]), 1)
  pie.df$label <- paste0(pie.df[, 2], "\n", "(", pie.df$percent, "%)")
  
  if (order) {
    pie.df[, 1] <- factor(pie.df[, 1], levels = levels)
    pie.df <- pie.df[order(pie.df[, 1]), ]
    #pie.df$label <- paste0(pie.df[, 1], "\n", pie.df[, 2], "\n", "(", pie.df$percent, "%)")
  }
  
  # Plot
  p.pie <- ggplot(pie.df, aes(x = "", y = .data[[selected_cols[2]]], fill = .data[[selected_cols[1]]])) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    theme_void() +
    geom_text(aes(label = label),
              position = position_stack(vjust = 0.5),
              size = 4,
              color = "black") +
    scale_fill_manual(values = colors,
                      labels = c("Alt. 3' Splice Sites (A3SS)", "Alt. 5' Splice Sites (A5SS)",
                                 "Mutually Exclusive Exons (MXE)", "Retained Intron (RI)", "Skipped Exon (SE)")) +
    guides(fill = guide_legend(title = "Splicing Event")) +
    theme(legend.title = element_text(size = 16, face = "bold"),
          legend.text = element_text(size = 14))
  
  return(p.pie)
}