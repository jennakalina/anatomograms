library(dplyr)
library(tidyr)
library(Seurat)

#seuObj <- readRDS('raw_data/trachea/ctrl_filtered.rds')
#condition <- 'control'
seuObj <- readRDS('raw_data/trachea/HSD_filtered.rds')
condition <- 'hsd'

# ASP - air sac primordium
# DB  - dorsal branch
# DT  - dorsal trunk
# SB  - spiracular branch
# PC  - progenitor cells
# TC  - transverse connective
# VB  - visceral branch
# LT  - lateral trunk
# GB  - ganglionic branches

mat <- seuObj@assays$RNA$data %>% as.matrix()
meta <- seuObj@meta.data
meta <- meta %>% tibble::rownames_to_column('cell') %>% rename(cluster = celltype)

# Get names of cells in each cluster
asp_cells <- meta$cell[meta$cluster == 'ASP']
db_cells <- meta$cell[meta$cluster == 'DB']
dt_cells <- meta$cell[meta$cluster == 'DT']
sb_cells <- meta$cell[meta$cluster == 'SB']
pc_cells <- meta$cell[meta$cluster == 'PC']
tc_cells <- meta$cell[meta$cluster == 'TC']
vb_cells <- meta$cell[meta$cluster == 'VB']
lt_cells <- meta$cell[meta$cluster == 'LT']
gb_cells <- meta$cell[meta$cluster == 'GB']

cells <- list(asp_cells, db_cells, dt_cells, sb_cells, pc_cells,
              tc_cells, vb_cells, lt_cells, gb_cells)
names <- list('air_sac_primordium', 'dorsal_branch', 'dorsal_trunk', 'spiracular_branch', 'progenitor_cells',
              'transverse_connective', 'visceral_branch', 'lateral_trunk', 'ganglionic_branches')
mats <- list()

# Make a subset for each cluster
for (i in 1:length(cells)) {
  sub_mat <- mat[, cells[[i]]]
  rownames(sub_mat) <- rownames(mat)
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:length(cells)) {
  matrix <- mats[[i]]
  
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression (percent of cells expressing the gene)
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
for (i in 1:length(cells)) {
  write.csv(dfs[[i]], paste0('results/trachea/', condition, '/', names[[i]], '.csv'), row.names = FALSE)
}
