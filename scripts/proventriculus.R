library(dplyr)
library(stringr)

mat <- read.csv('raw_data/proventriculus_p50/FC_08073.mat.csv', check.names = FALSE)
meta <- read.csv('raw_data/proventriculus_p50/FC_08073.meta.csv')
colnames(meta) <- c('cell', 'celltype', 'condition')
mat <- mat %>% tibble::column_to_rownames('X')

# Separate into wild type and Yki
wt_cells <- meta$cell[meta$condition == 'wt']
wt_mat <- mat[, wt_cells]
wt_meta <- meta[meta$condition == 'wt', ]
wt_meta <- wt_meta %>% mutate(celltype = gsub('^[^[:alpha:]]+', '', celltype))

yki_cells <- meta$cell[meta$condition == 'yki']
yki_mat <- mat[, yki_cells]
yki_meta <- meta[meta$condition == 'yki', ]
yki_meta <- yki_meta %>% mutate(celltype = gsub('^[^[:alpha:]]+', '', celltype))

### First wt
# Get names of cells in each cluster
esophagus_cells <- wt_meta$cell[wt_meta$celltype == 'esophagus']
muscle_cells <- wt_meta$cell[wt_meta$celltype == 'muscles']
zone6_cells <- wt_meta$cell[wt_meta$celltype == 'zone6']
zone56border_cells <- wt_meta$cell[wt_meta$celltype == 'zone5_6_border_cells']
zone5_cells <- wt_meta$cell[wt_meta$celltype == 'zone5_dome']
zone4_cells <- wt_meta$cell[wt_meta$celltype == 'zone4']
zone23_cells <- wt_meta$cell[wt_meta$celltype == 'zone2_3']
zone1_cells <- wt_meta$cell[wt_meta$celltype == 'zone1']

cells <- list(esophagus_cells, muscle_cells, zone6_cells, zone56border_cells,
              zone5_cells, zone4_cells, zone23_cells, zone1_cells)
names <- list('esophagus', 'muscles', 'zone6', 'zone56border', 
              'zone5', 'zone4', 'zone23', 'zone1')
mats <- list()

# Make a subset for each cluster
wt_mat <- wt_mat %>% tibble::rownames_to_column('geneID')
for (i in 1:8) {
  sub_mat <- wt_mat[, cells[[i]]]
  rownames(sub_mat) <- wt_mat$geneID
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:8) {
  matrix <- mats[[i]]
  
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression (percent of cells expressing the gene)
  # CPM
  total_cts <- sum(matrix)
  scale_fac <- 1e6 / total_cts
  cpm <- rowSums(matrix) * scale_fac
  
  df <- data.frame(gene = rownames(matrix),
                   fbgn = NA,
                   average_count = avg,
                   pct_expression = pct,
                   CPM = cpm)
  dfs[[i]] <- df
}

# Save all to csv
for (i in 1:8) {
  write.csv(dfs[[i]], paste0('results/proventriculus_p50/wt/', names[[i]], '.csv'), row.names = FALSE)
}

### Now Yki
# Get names of cells in each cluster
esophagus_cells <- yki_meta$cell[yki_meta$celltype == 'esophagus']
muscle_cells <- yki_meta$cell[yki_meta$celltype == 'muscles']
zone6_cells <- yki_meta$cell[yki_meta$celltype == 'zone6']
zone56border_cells <- yki_meta$cell[yki_meta$celltype == 'zone5_6_border_cells']
zone5_cells <- yki_meta$cell[yki_meta$celltype == 'zone5_dome']
zone4_cells <- yki_meta$cell[yki_meta$celltype == 'zone4']
zone23_cells <- yki_meta$cell[yki_meta$celltype == 'zone2_3']
zone1_cells <- yki_meta$cell[yki_meta$celltype == 'zone1']

cells <- list(esophagus_cells, muscle_cells, zone6_cells, zone56border_cells,
              zone5_cells, zone4_cells, zone23_cells, zone1_cells)
names <- list('esophagus', 'muscles', 'zone6', 'zone56border', 
              'zone5', 'zone4', 'zone23', 'zone1')
mats <- list()

# Make a subset for each cluster
yki_mat <- yki_mat %>% tibble::rownames_to_column('geneID')
for (i in 1:8) {
  sub_mat <- yki_mat[, cells[[i]]]
  rownames(sub_mat) <- yki_mat$geneID
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:8) {
  matrix <- mats[[i]]
  
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # pct expression
  # CPM
  total_cts <- sum(matrix)
  scale_fac <- 1e6 / total_cts
  cpm <- rowSums(matrix) * scale_fac
  
  df <- data.frame(gene = rownames(matrix),
                   fbgn = NA,
                   average_count = avg,
                   pct_expression = pct,
                   CPM = cpm)
  dfs[[i]] <- df
}

# Save all to csv
for (i in 1:8) {
  write.csv(dfs[[i]], paste0('results/proventriculus_p50/yki/', names[[i]], '.csv'), row.names = FALSE)
}
