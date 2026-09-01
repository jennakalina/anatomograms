library(loomR)
library(hdf5r)
library(dplyr)

lfile <- connect(filename = 'r_fca_biohub_all_wo_blood_10x.loom', mode = 'r+', skip.validate = TRUE)
lfile

table(lfile$col.attrs$annotation[]) # Check counts

# Get annotations, gene names, cell names
annotations <- lfile$col.attrs$annotation[]
cell_names <- lfile$col.attrs$CellID[]
gene_names <- lfile$row.attrs$Gene[]

### Posterior file
posterior_groups <- list(
  enteroblast = 'enteroblast',
  intestinal_stem = 'intestinal stem cell',
  enterocyte_post_epi = 'enterocyte of posterior adult midgut epithelium',
  enteroendocrine = 'enteroendocrine cell',
  visceral_midgut = 'visceral muscle of the midgut'
)

### Middle file
middle_groups <- list(
  visceral_midgut = 'visceral muscle of the midgut',
  enteroblast = 'enteroblast',
  intestinal_stem = 'intestinal stem cell',
  large_flat = 'midgut large flat cell',
  enteroendocrine = 'enteroendocrine cell',
  copper = 'copper cell',
  enterocyte = 'adult midgut enterocyte'
)

### Anterior file
anterior_groups <- list(
  enteroblast = 'enteroblast',
  intestinal_stem = 'intestinal stem cell',
  enteroendocrine = 'enteroendocrine cell',
  visceral_midgut = 'visceral muscle of the midgut',
  enterocyte_ant_epi = 'enterocyte of anterior adult midgut epithelium'
)

# Function to calculate stats for each group
process_groups <- function(lfile, annotations, group_list, gene_names, outdir) {
  
  for (group_name in names(group_list)) {
    
    labels <- group_list[[group_name]]
    
    # Get cells to keep
    cells_idx <- which(annotations %in% labels)
    
    if (length(cells_idx) == 0) {
      message(paste("Skipping", group_name, "; no cells found"))
      next
    }
    
    message(paste("Processing", group_name))
    
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
    
    write.csv(df, paste0('../../results/fca/midgut/', outdir, '/', group_name, '.csv'), 
              row.names = FALSE)
    
    rm(sub_mat)
  }
}

# Run for posterior
process_groups(lfile, annotations, posterior_groups, gene_names, 'posterior')

# Run for middle
process_groups(lfile, annotations, middle_groups, gene_names, 'middle')

# Run for anterior
process_groups(lfile, annotations, anterior_groups, gene_names, 'anterior')
