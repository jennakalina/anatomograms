library(loomR)
library(hdf5r)

lfile <- connect(filename = "r_fca_biohub_all_wo_blood_10x.loom", mode = "r+", skip.validate = TRUE)
lfile

#table(lfile$col.attrs$annotation[]) # Check different objects to make sure it's the right annotation

clusters_keep <- c('adult heart', 'adult abdominal pericardial cell', 
                   'adult ostium', 'cardiomyocyte, working adult heart muscle (non-ostia)')

# Get cells that belong to just these clusters
annotations <- lfile$col.attrs$annotation[]
cells_keep <- annotations %in% clusters_keep
# Check it's the right number
table(cells_keep)

# Subset the matrix
mat_subset <- lfile[['matrix']][which(cells_keep), ]
# Get gene col names
gene_names <- lfile$row.attrs$Gene[]
colnames(mat_subset) <- gene_names
# Get cell row names
cell_names <- lfile$col.attrs$CellID[]
cell_names_subset <- cell_names[which(cells_keep)]
rownames(mat_subset) <- cell_names_subset
mat <- t(mat_subset)

# Separate further into four clusters
annotations_subset <- annotations[which(cells_keep)]

mats <- lapply(clusters_keep, function(clust){
  mat[, annotations_subset == clust, drop = FALSE]
})
names(mats) <- clusters_keep

dfs <- list()
for (i in 1:4) {
  matrix <- mats[[i]]
  
  avg <- rowMeans(matrix) # average count
  pct <- rowSums(matrix != 0) / ncol(matrix) # percent expression 
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
names <- c('adult_heart', 'adult_apc', 'adult_ostium', 'cardiomyocyte')
for (i in 1:4) {
  write.csv(dfs[[i]], paste0('../../results/fca/heart/', names[[i]], '.csv'), row.names = FALSE)
}


