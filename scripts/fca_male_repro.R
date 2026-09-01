library(loomR)
library(hdf5r)

lfile <- connect(filename = 'r_fca_biohub_all_wo_blood_10x.loom', mode = 'r+', skip.validate = TRUE)
lfile

table(lfile$col.attrs$annotation[]) # Check counts

# Get annotations, gene names, cell names
annotations <- lfile$col.attrs$annotation[]
cell_names <- lfile$col.attrs$CellID[]
gene_names <- lfile$row.attrs$Gene[]

### Testis file
testis_groups <- list(
  testis = 'testis',
  spermatozoon = 'mid-late elongation-stage spermatid',
  mid_late_prolif = 'mid-late proliferating spermatogonia',
  spermatid_cyst = c('early elongation stage spermatid','early-mid elongation-stage spermatid','mid-late proliferating spermatogonia'),
  spermatogonium = 'spermatogonium',
  cyst_progenitor = c('cyst stem cell', 'early cyst cell 1', 'early cyst cell 2'),
  male_germline = c('spermatogonium', 'mid-late proliferating spermatogonia'),
  spermatocyte_cyst = c('spermatocyte cyst cell branch a', 'spermatocyte cyst cell branch b'),
  prim_spermatocyte = c('late primary spermatocyte', 'spermatocyte 1', 'spermatocyte 2', 'spermatocyte 3', 'spermatocyte 4'),
  sec_spermatocyte = c('spermatocyte 5', 'spermatocyte 6', 'spermatocyte 7a'),
  spermatid = 'spermatid',
  tail_cyst = c('late cyst cell branch a', 'late cyst cell branch b'),
  head_cyst = 'head cyst cell',
  cyst_testis = c('cyst cell intermediate', 'cyst cell branch a', 'cyst cell branch b')
)

### Male Reproductive System file
male_rep_groups <- list(
  male_rep_sys = c('secretory cell of the male reproductive tract', 'male gonad associated epithelium'),
  ejac_duct = 'anterior ejaculatory duct',
  ejac_bulb = 'ejaculatory bulb',
  access_gland = 'male accessory gland',
  sem_vesicle = 'seminal vesicle',
  testis = 'testis'
)

# Function to calculate stats for each group
process_groups <- function(lfile, annotations, group_list, gene_names, outdir) {
  
  for (group_name in names(group_list)) {
    
    labels <- group_list[[group_name]]
    
    # Get cells to keep
    cells_idx <- which(annotations %in% labels)
    
    if (length(cells_idx) == 0) {
      message(paste("Skipping", group_name, "- no cells found"))
      next
    }
    
    message(paste("Processing:", group_name))
    
    # Subset from matrix
    sub_mat <- lfile[['matrix']][cells_idx, ]
    
    # Matrix is in cells x genes; switch
    sub_mat <- t(sub_mat)
    
    # Calculate stats
    avg <- rowMeans(sub_mat)
    pct <- rowSums(sub_mat != 0) / ncol(sub_mat)
    
    total_cts <- sum(sub_mat)
    scale_fac <- 1e6 / total_cts
    cpm <- rowSums(sub_mat) * scale_fac
    
    df <- data.frame(gene = gene_names,
                     fbgn = NA,
                     average_count = avg,
                     pct_expression = pct,
                     CPM = cpm)
    
    write.csv(df, paste0('../../results/fca/male_reproduction/', outdir, '/', group_name, '.csv'), 
              row.names = FALSE)
    
    rm(sub_mat)
  }
}

# Run for testis
process_groups(lfile, annotations, testis_groups, gene_names, 'testis')

# Run for male reproductive system
process_groups(lfile, annotations, male_rep_groups, gene_names, 'male_repro_sys')

