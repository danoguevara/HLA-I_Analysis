nice_Hist <- function(counts, mean = "Mean", thresholds = c(seq(0, 2, 0.5), 3, 5), xlim = 50, bin = 1) {
  levels_vector <- c(paste0("= ", thresholds[1]), paste0("<= ", thresholds[2]), paste0("<= ", thresholds[3]),
                     paste0("<= ", thresholds[4]), paste0("<= ", thresholds[5]), paste0("<= ", thresholds[6]),
                     paste0("<= ", thresholds[7]), paste0("> ", thresholds[7]))
  
  counts <- counts %>% mutate(category = case_when(.data[[mean]] == thresholds[1] ~ levels_vector[1],
                                                   .data[[mean]] <= thresholds[2] ~ levels_vector[2],
                                                   .data[[mean]] <= thresholds[3] ~ levels_vector[3],
                                                   .data[[mean]] <= thresholds[4] ~ levels_vector[4],
                                                   .data[[mean]] <= thresholds[5] ~ levels_vector[5],
                                                   .data[[mean]] <= thresholds[6] ~ levels_vector[6],
                                                   .data[[mean]] <= thresholds[7] ~ levels_vector[7],
                                                   TRUE ~ levels_vector[8]))
  
  counts$category <- factor(counts$category, levels = levels_vector)
  
  #fill_colors <- setNames(c("#00004F", "#6000C6", "#B3A3FF", "steelblue", "#608FE6", "#85C0FF", "grey70", "white"), levels_vector)
  #fill_colors <- setNames(c("#08306b", "#08519c", "#2171b5", "#4292c6", "#6baed6", "#9ecae1", "grey70", "white"), levels_vector)
  fill_colors <- setNames(c("#440154", "#443A83", "#31688E", "#21918C", "#35B779", "#FDE725", "grey70", "white"), levels_vector)
  
  plot <- ggplot(counts[counts[[mean]] <= xlim, ], aes(x = .data[[mean]], fill = category)) +
    geom_histogram(binwidth = bin, boundary = 0, closed = "left", color = "black", position = "stack") + theme_bw() +
    xlab("Mean Normalized Counts") + ylab("Frequency") +
    scale_x_continuous(breaks = seq(0, xlim, by = 5)) +
    scale_fill_manual(values = fill_colors, name = "Counts") +
    theme(strip.text = element_text(size = 14), axis.text = element_text(size = 12),
          axis.title = element_text(size = 16), legend.title = element_text(size = 14),
          legend.text = element_text(size = 12))
  
  return(plot)
}
