library(dplyr)
library(stringr)

# Read in data
mat <- read.csv('raw_data/kidney_p49/normalized_matrix.csv', check.names = FALSE)
meta <- read.csv('raw_data/kidney_p49/meta.csv')
colnames(meta)[1] <- 'cell'
mat <- mat %>% tibble::column_to_rownames('X')

# Remove clusters we won't use
clus_to_remove <- c('13 adult tracheal cells', '15 visceral muscle', '16 fat body', '9 principal cells-like')

meta_filt <- meta %>% filter(!annotation %in% clus_to_remove)

# Separate into wild type and Yki
wt_cells <- meta_filt$cell[meta_filt$sample == 'EGT>w1118 for Yki']
wt_mat <- mat[, wt_cells]
wt_meta <- meta_filt[meta_filt$sample == 'EGT>w1118 for Yki', ]
wt_meta <- wt_meta %>% mutate(annotation = gsub('^[^[:alpha:]]+', '', annotation))

yki_cells <- meta_filt$cell[meta_filt$sample == 'Yki']
yki_mat <- mat[, yki_cells]
yki_meta <- meta_filt[meta_filt$sample == 'Yki', ]
yki_meta <- yki_meta %>% mutate(annotation = gsub('^[^[:alpha:]]+', '', annotation))

### First wt
# Get names of cells in each cluster
lower_segment_cells <- wt_meta$cell[wt_meta$annotation == 'lower tubule principal cells' | wt_meta$annotation == 'lower tubule principal cells (Pvr-act)' |
                                      wt_meta$annotation == 'lower segment principal cells (Pvr-act)' | wt_meta$annotation == 'lower segment principal cells 2 (Pvr-act)']
init_trans_pc_cells <- wt_meta$cell[wt_meta$annotation == 'initial and transitional principal cells']
lower_uterer_pc_cells <- wt_meta$cell[wt_meta$annotation == 'lower ureter principal cells']
main_segment_pc_cells <- wt_meta$cell[wt_meta$annotation == 'main segment principal cells' | wt_meta$annotation == 'main segment principal cells 2']
main_segment_stellate_cells <- wt_meta$cell[wt_meta$annotation == 'stellate cells']
stem_cell_cells <- wt_meta$cell[wt_meta$annotation == 'renal stem cells' | wt_meta$annotation == 'renal stem cells (Pvr-act)' |
                                  wt_meta$annotation == 'renal stem cells 2 (Pvr-act)']
upper_uterer_pc_cells <- wt_meta$cell[wt_meta$annotation == 'upper ureter principal cells']

cells <- list(lower_segment_cells,init_trans_pc_cells,lower_uterer_pc_cells,main_segment_pc_cells,
              main_segment_stellate_cells,stem_cell_cells,upper_uterer_pc_cells)
names <- list('lower_segment', 'init_trans_pc', 'lower_uterer_pc', 'main_segment_pc', 
              'main_segment_stellate', 'stem_cell', 'upper_uterer_pc')
mats <- list()

# Make a subset for each cluster
wt_mat <- wt_mat %>% tibble::rownames_to_column('geneID')
for (i in 1:7) {
  sub_mat <- wt_mat[, cells[[i]]]
  rownames(sub_mat) <- wt_mat$geneID
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:7) {
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
for (i in 1:7) {
  write.csv(dfs[[i]], paste0('results/kidney_p49/wt_egt/', names[[i]], '.csv'), row.names = FALSE)
}

### Now Yki
# Get names of cells in each cluster
lower_segment_cells <- yki_meta$cell[yki_meta$annotation == 'lower tubule principal cells' | yki_meta$annotation == 'lower tubule principal cells (Pvr-act)' |
                                       yki_meta$annotation == 'lower segment principal cells (Pvr-act)' | yki_meta$annotation == 'lower segment principal cells 2 (Pvr-act)']
init_trans_pc_cells <- yki_meta$cell[yki_meta$annotation == 'initial and transitional principal cells']
lower_uterer_pc_cells <- yki_meta$cell[yki_meta$annotation == 'lower ureter principal cells']
main_segment_pc_cells <- yki_meta$cell[yki_meta$annotation == 'main segment principal cells' | yki_meta$annotation == 'main segment principal cells 2']
main_segment_stellate_cells <- yki_meta$cell[yki_meta$annotation == 'stellate cells']
stem_cell_cells <- yki_meta$cell[yki_meta$annotation == 'renal stem cells' | yki_meta$annotation == 'renal stem cells (Pvr-act)' |
                                   yki_meta$annotation == 'renal stem cells 2 (Pvr-act)']
upper_uterer_pc_cells <- yki_meta$cell[yki_meta$annotation == 'upper ureter principal cells']

cells <- list(lower_segment_cells,init_trans_pc_cells,lower_uterer_pc_cells,main_segment_pc_cells,
              main_segment_stellate_cells,stem_cell_cells,upper_uterer_pc_cells)
names <- list('lower_segment', 'init_trans_pc', 'lower_uterer_pc', 'main_segment_pc', 
              'main_segment_stellate', 'stem_cell', 'upper_uterer_pc')
mats <- list()

# Make a subset for each cluster
yki_mat <- yki_mat %>% tibble::rownames_to_column('geneID')
for (i in 1:7) {
  sub_mat <- yki_mat[, cells[[i]]]
  rownames(sub_mat) <- yki_mat$geneID
  mats[[i]] <- sub_mat
}

dfs <- list()
for (i in 1:7) {
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
for (i in 1:7) {
  write.csv(dfs[[i]], paste0('results/kidney_p49/yki/', names[[i]], '.csv'), row.names = FALSE)
}

### Make loop for last three
conditions <- c('w1118 for Pvr', 'Pvr Activation', 'Pvr RNAi')
cond_names <- list('wt_pvr', 'pvr_activated', 'pvr_rnai')

for (j in 1:length(conditions)) {
  condition <- conditions[[j]]
  
  # Filter down to condition
  cond_cells <- meta_filt$cell[meta_filt$sample == condition]
  cond_mat <- mat[, cond_cells]
  cond_meta <- meta_filt[meta_filt$sample == condition, ]
  cond_meta <- cond_meta %>% mutate(annotation = gsub('^[^[:alpha:]]+', '', annotation))
  
  # Get names of cells in each cluster
  lower_segment_cells <- cond_meta$cell[cond_meta$annotation == 'lower tubule principal cells' | cond_meta$annotation == 'lower tubule principal cells (Pvr-act)' |
                                          cond_meta$annotation == 'lower segment principal cells (Pvr-act)' | cond_meta$annotation == 'lower segment principal cells 2 (Pvr-act)']
  init_trans_pc_cells <- cond_meta$cell[cond_meta$annotation == 'initial and transitional principal cells']
  lower_uterer_pc_cells <- cond_meta$cell[cond_meta$annotation == 'lower ureter principal cells']
  main_segment_pc_cells <- cond_meta$cell[cond_meta$annotation == 'main segment principal cells' | cond_meta$annotation == 'main segment principal cells 2']
  main_segment_stellate_cells <- cond_meta$cell[cond_meta$annotation == 'stellate cells']
  stem_cell_cells <- cond_meta$cell[cond_meta$annotation == 'renal stem cells' | cond_meta$annotation == 'renal stem cells (Pvr-act)' |
                                      cond_meta$annotation == 'renal stem cells 2 (Pvr-act)']
  upper_uterer_pc_cells <- cond_meta$cell[cond_meta$annotation == 'upper ureter principal cells']
  
  cells <- list(lower_segment_cells,init_trans_pc_cells,lower_uterer_pc_cells,main_segment_pc_cells,
                main_segment_stellate_cells,stem_cell_cells,upper_uterer_pc_cells)
  names <- list('lower_segment', 'init_trans_pc', 'lower_uterer_pc', 'main_segment_pc', 
                'main_segment_stellate', 'stem_cell', 'upper_uterer_pc')
  mats <- list()
  
  # Make a subset for each cluster
  cond_mat <- cond_mat %>% tibble::rownames_to_column('geneID')
  for (i in 1:7) {
    if (length(cells[[i]]) == 0) {
      message(paste('Skipping', names[i], 'in condition', condition))
      next
    }
    
    sub_mat <- cond_mat[, c("geneID", cells[[i]]), drop = FALSE]
    rownames(sub_mat) <- sub_mat$geneID
    sub_mat$geneID <- NULL
    
    mats[[i]] <- sub_mat
  }
  
  dfs <- list()
  for (i in 1:7) {
    if (is.null(mats[[i]])) next
    
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
  for (i in 1:7) {
    if (is.null(dfs[[i]])) next
    
    write.csv(dfs[[i]], paste0('results/kidney_p49/', cond_names[[j]], '/', names[[i]], '.csv'), row.names = FALSE)
  }
}

