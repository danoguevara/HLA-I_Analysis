intersect_SE <- function(df_list, indexes, id_col) {
  if (length(indexes) != length(df_list)) {
    stop("`indexes` must have the same amount of elements as `df_list`")
  }
  
  names(df_list) <- indexes
  
  ids_1 <- df_list[[indexes[1]]][, id_col]
  ids_2 <- df_list[[indexes[2]]][, id_col]
  ids_3 <- df_list[[indexes[3]]][, id_col]
  ids_4 <- df_list[[indexes[4]]][, id_col]
  ids_5 <- df_list[[indexes[5]]][, id_col]
  
  # Pairwise intersections
  Reduce(intersect, list(ids_1, ids_2))
  Reduce(intersect, list(ids_1, ids_3))
  Reduce(intersect, list(ids_1, ids_4))
  Reduce(intersect, list(ids_1, ids_5))
  
  Reduce(intersect, list(ids_2, ids_3))
  Reduce(intersect, list(ids_2, ids_4))
  Reduce(intersect, list(ids_2, ids_5))
  
  Reduce(intersect, list(ids_3, ids_4))
  Reduce(intersect, list(ids_3, ids_5))
  
  Reduce(intersect, list(ids_4, ids_5))
  
  # 3-way comparisons
  Reduce(intersect, list(ids_1, ids_2, ids_3))
  Reduce(intersect, list(ids_1, ids_2, ids_4))
  Reduce(intersect, list(ids_1, ids_2, ids_5))
  
  Reduce(intersect, list(ids_2, ids_3, ids_4))
  Reduce(intersect, list(ids_2, ids_3, ids_5))
  
  Reduce(intersect, list(ids_3, ids_4, ids_5))
  
  # 4-way comparisons
  Reduce(intersect, list(ids_1, ids_2, ids_3, ids_4))
  Reduce(intersect, list(ids_2, ids_3, ids_4, ids_5))
  
  # Full set comparisons
  Reduce(intersect, list(ids_1, ids_2, ids_3, ids_4, ids_5))
}


