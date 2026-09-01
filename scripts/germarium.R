library(dplyr)
library(stringr)
library(tidyr)
library(Matrix)

# Read in data
meta <- read.csv('raw_data/ovary_CA/metadata.csv') %>% rename(cell = barcode)
# Put matrix together
mat <- readMM('raw_data/ovary_CA/normalised_counts.mtx')
rows <- read.table('raw_data/ovary_CA/normalised_counts.mtx_rows', stringsAsFactors = FALSE)[,1]
cols <- readLines('raw_data/ovary_CA/normalised_counts.mtx_cols')

rownames(mat) <- rows
colnames(mat) <- cols
mat_dense <- as.matrix(mat)

# Get names of cells in each cluster
ant_escort_cells <- meta$cell[meta$cluster == 'anterior escort cell']
cent_escort_cells <- meta$cell[meta$cluster == 'central escort cell']
post_escort_cells <- meta$cell[meta$cluster == 'posterior escort cell']
follicle_stem_cells <- meta$cell[meta$cluster == 'follicle stem cell']
germline_stem_cells <- meta$cell[meta$cluster == 'germline stem cell']
germarium_cap_cells <- meta$cell[meta$cluster == 'cap cell']
term_filament_cells <- meta$cell[meta$cluster == 'terminal filament cell']
prefollicle_cells <- meta$cell[meta$cluster == 'early prefollicle cell' | meta$cluster == 'late prefollicle cell']
follicle_cells <- meta$cell[meta$cluster == 'stage 2-5 main body follicle cell' | 
                              meta$cluster == 'stage 5-6 anterior/central main body follicle cell' |
                              meta$cluster == 'stage 6-7 central main body follicle cell' |
                              meta$cluster == 'stage 8 central main body follicle cell' |
                              meta$cluster == 'stage 9+ main body follicle cell' |
                              meta$cluster == 'stage 8 posterior main body follicle cell' |
                              meta$cluster == 'stage 5-6 posterior main body follicle cell' |
                              meta$cluster == 'stage 7 posterior main body follicle cell']
cytoblast_cells <- meta$cell[meta$cluster == 'undifferentiated germ cell']
fg2_cells <- meta$cell[meta$cluster == 'undifferentiated germ cell']
fg4_cells <- meta$cell[meta$cluster == 'undifferentiated germ cell']
fg8_cells <- meta$cell[meta$cluster == 'undifferentiated germ cell']
fg16_cells <- meta$cell[meta$cluster == 'older germ cell']


# Make list of cluster cells and names
cluster_cells <- list(ant_escort_cells, cent_escort_cells, post_escort_cells, follicle_stem_cells, 
                      germline_stem_cells, germarium_cap_cells, term_filament_cells, prefollicle_cells, 
                      follicle_cells, cytoblast_cells, fg2_cells, fg4_cells, fg8_cells, fg16_cells)
cluster_names <- c('anterior_escort_cells', 'central_escort_cells', 'posterior_escort_cells', 
                   'follicle_stem_cells', 'germline_stem_cells', 'germarium_cap_cell', 
                   'terminal_filament_cells', 'prefollicle_cells', 'follicle_cells', 'cytoblast', 
                   'f_germline_2_cell_cyst', 'f_germline_4_cell_cyst', 'f_germline_8_cell_cyst',
                   'f_germline_16_cell_cyst')
  
# Loop through clusters, make df for each
for (j in 1:length(cluster_cells)) {
  cells <- cluster_cells[[j]]
  name <- cluster_names[[j]]
  
  matrix <- mat_dense[, cells]
    
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression 
  # CPM
  total_cts <- sum(matrix)
  scale_fac <- 1e6 / total_cts
  cpm <- rowSums(matrix) * scale_fac
    
  df <- data.frame(fbgn = rownames(matrix),
                   average_count = avg,
                   pct_expression = pct,
                   CPM = cpm)
    
  write.csv(df, paste0('results/ovary_ca/', name, '.csv'), row.names = FALSE, quote = FALSE)
}
