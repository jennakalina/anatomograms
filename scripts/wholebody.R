library(dplyr)
library(stringr)
library(tidyr)

mat <- readRDS('raw_data/p61/matrix_fixed_mapping.rds')
meta <- read.delim('raw_data/p61/metadata.tsv')
colnames(meta) <- c('cell', 'UMAP_1', 'UMAP_2', 'celltype', 'condition')

# Separate into conditions
conditions <- c('5D_Wt', '5D_Yki', '8D_Wt', '8D_Yki')
mats <- list()
metas <- list()

for (i in 1:length(conditions)) {
  cond <- conditions[[i]]
  
  cells <- meta$cell[meta$condition == cond]
  sub_mat <- mat[, cells]
  sub_meta <- meta[meta$condition == cond, ] %>%
    mutate(celltype = gsub('^[^[:alpha:]]+', '', celltype))
  
  mats[[cond]] <- sub_mat
  metas[[cond]] <- sub_meta
}

# Loop to run through each condition and get each cluster
for (i in 1:length(conditions)) {
  cond <- conditions[[i]]
  sub_mat <- mats[[cond]]
  sub_meta <- metas[[cond]]
  
  # Get names of cells in each cluster
  epithelial_cells <- sub_meta$cell[sub_meta$celltype == 'epithelial_cell']
  muscle_cells <- sub_meta$cell[sub_meta$celltype == 'muscle_cell' | sub_meta$celltype == 'indirect_flight_muscle']
  heart_cells <- sub_meta$cell[sub_meta$celltype == 'heart_muscle']
  visceral_muscle_cells <- sub_meta$cell[sub_meta$celltype == 'visceral_muscle_of_midgut']
  intestine_stem_cells <- sub_meta$cell[sub_meta$celltype == 'intestinal_stem_cell']
  gut_ec_cells <- sub_meta$cell[sub_meta$celltype == 'enterocyte' | sub_meta$celltype == 'enterocyte-like']
  gut_ee_cells <- sub_meta$cell[sub_meta$celltype == 'enteroendocrine_cell']
  cardia_cells <- sub_meta$cell[sub_meta$celltype == 'cardia']
  hindgut_cells <- sub_meta$cell[sub_meta$celltype == 'hindgut_A' | 
                                   sub_meta$celltype == 'hindgut_B' | 
                                   sub_meta$celltype == 'hindgut_C']
  fat_body_cells <- sub_meta$cell[sub_meta$celltype == 'fat_body']
  oenocyte_cells <- sub_meta$cell[sub_meta$celltype == 'oenocyte']
  hemocyte_cells <- sub_meta$cell[sub_meta$celltype == 'hemocyte']
  trachea_cells <- sub_meta$cell[sub_meta$celltype == 'tracheal_cell']
  malpighian_tubule_cells <- sub_meta$cell[sub_meta$celltype == 'malpighian_tubule_principal_cell']
  salivary_gland_cells <- sub_meta$cell[sub_meta$celltype == 'salivary_gland']
  ventral_nervous_system_cells <- sub_meta$cell[sub_meta$celltype == 'ventral_nervous_system']
  glial_cells <- sub_meta$cell[sub_meta$celltype == 'cell_body_glial_cell']
  female_repro_cells  <- sub_meta$cell[sub_meta$celltype == 'cyst_stem_cells' | sub_meta$celltype == 'female_reproductive_system_A' |
                                         sub_meta$celltype == 'female_reproductive_system_B' | sub_meta$celltype == 'follicle_cell' |
                                         sub_meta$celltype == 'follicle_cell_stage 9+*' | sub_meta$celltype == 'germline_cell_A_unknown' |
                                         sub_meta$celltype == 'germline_cell_B_germline' | sub_meta$celltype == 'germline_cell_C_unknown stage' |
                                         sub_meta$celltype == 'stalk_follicle_cell' | sub_meta$celltype == 'stretch_follicle_cell' |
                                         sub_meta$celltype == 'oviduct']
  neuron_cells  <- sub_meta$cell[sub_meta$celltype == 'neuron_A' | sub_meta$celltype == 'neuron_B']
  
  # Make list of cluster cells and names
  cluster_cells <- list(epithelial_cells, muscle_cells, heart_cells, visceral_muscle_cells, intestine_stem_cells,
                     gut_ec_cells, gut_ee_cells, cardia_cells, hindgut_cells, fat_body_cells, oenocyte_cells,
                     hemocyte_cells, trachea_cells, malpighian_tubule_cells, salivary_gland_cells, 
                     ventral_nervous_system_cells, glial_cells, female_repro_cells, neuron_cells)
  cluster_names <- c('epithelial_cell', 'muscle', 'heart', 'visceral_muscle_midgut', 'intestine_stem_cell',
                     'gut_EC_cell', 'gut_EE_cell', 'cardia', 'hindgut', 'fat_body', 'oenocyte', 'hemocyte',
                     'trachea', 'malpighian_tubule', 'salivary_gland', 'ventral_nervous_system', 'glial',
                     'female_reproduction', 'neuron')
  
  # Loop through clusters, make df for each
  for (j in 1:length(cluster_cells)) {
    cells <- cluster_cells[[j]]
    name <- cluster_names[[j]]
    
    # Skip if no cells in cluster
    if (length(cells) == 0) next
    
    matrix <- sub_mat[, cells]
    
    avg <- rowMeans(matrix) # average count
    pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression (percent of cells expressing the gene)
    # CPM
    total_cts <- sum(matrix)
    scale_fac <- 1e6 / total_cts
    cpm <- rowSums(matrix) * scale_fac
    
    df <- data.frame(fbgn = rownames(matrix),
                     average_count = avg,
                     pct_expression = pct,
                     CPM = cpm)
    
    write.csv(df, paste0('results/p61/', cond, '/', name, '.csv'), row.names = FALSE, quote = FALSE)
  }
}

