library(dplyr)
library(stringr)

mat <- read.csv('raw_data/kidney_p38/matrix_kidney.csv', check.names = FALSE)
meta <- read.csv('raw_data/kidney_p38/metadata_kidney.csv')
colnames(meta)[1] <- 'cell'

# Get names of cells in each cluster
lower_segment_cells <- meta$cell[meta$cluster == 'lower segment PC' | meta$cluster == 'lower tubule PC']
bar_shaped_stellate_cells <- meta$cell[meta$cluster == 'bar-shape stellate cell']
init_trans_pc_cells <- meta$cell[meta$cluster == 'initial and transitional PC']
lower_uterer_pc_cells <- meta$cell[meta$cluster == 'lower ureter PC']
main_segment_pc_cells <- meta$cell[meta$cluster == 'main segment PC']
main_segment_stellate_cells <- meta$cell[meta$cluster == 'main segment stellate cell']
stem_cell_cells <- meta$cell[meta$cluster == 'stem cell']
upper_uterer_pc_cells <- meta$cell[meta$cluster == 'upper ureter PC']

cells <- list(lower_segment_cells,bar_shaped_stellate_cells,init_trans_pc_cells,lower_uterer_pc_cells,
              main_segment_pc_cells,main_segment_stellate_cells,stem_cell_cells,upper_uterer_pc_cells)
names <- list('lower_segment', 'bar_shaped_stellate', 'init_trans_pc', 'lower_uterer_pc',
              'main_segment_pc', 'main_segment_stellate', 'stem_cell', 'upper_uterer_pc')
mats <- list()

# Make a subset for each cluster
for (i in 1:8) {
  sub_mat <- mat[, cells[[i]]]
  rownames(sub_mat) <- mat$geneID
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:8) {
  matrix <- mats[[i]]
  
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression 
  # CPM
  total_cts <- sum(matrix)
  scale_fac <- 1e6 / total_cts
  cpm <- rowSums(matrix) * scale_fac
  
  df <- data.frame(gene = rownames(matrix),
                   average_count = avg,
                   pct_expression = pct,
                   CPM = cpm)
  
  dfs[[i]] <- df
}

# Save all to csv
for (i in 1:8) {
  write.csv(dfs[[i]], paste0('data/kidney_subsets_agg/', names[[i]], '.csv'), row.names = FALSE)
}
