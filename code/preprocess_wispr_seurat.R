library(Seurat)
library(dplyr)
library(data.table)
library(readr)
library(Matrix)

# Set working directory
setwd("C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject")

# ============================================================================
# FUNCTION: Calculate Gini Coefficient for DEG Selection
# ============================================================================
# This replaces Giotto's findMarkers_one_vs_all with method='gini'

calculate_gini <- function(x) {
  # Gini coefficient measures inequality in gene expression
  # Higher Gini = more unequal = better marker gene
  n <- length(x)
  x_sorted <- sort(x)
  gini <- (2 * sum((1:n) * x_sorted)) / (n * sum(x_sorted)) - (n + 1) / n
  return(gini)
}

find_gini_markers <- function(seurat_obj, cluster_column = "seurat_clusters", 
                              top_n = 100) {
  # Get expression matrix (normalized)
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  
  # Get cluster assignments
  clusters <- seurat_obj@meta.data[[cluster_column]]
  unique_clusters <- unique(clusters)
  
  marker_list <- list()
  
  for (cluster in unique_clusters) {
    cat("Finding markers for cluster:", cluster, "\n")
    
    # Cells in this cluster vs others
    cells_in <- which(clusters == cluster)
    cells_out <- which(clusters != cluster)
    
    # Calculate Gini coefficient for each gene
    gini_scores <- apply(expr_matrix, 1, function(gene_expr) {
      # Expression in cluster
      expr_in <- gene_expr[cells_in]
      # Gini coefficient
      gini <- calculate_gini(expr_in)
      return(gini)
    })
    
    # Rank genes by Gini coefficient
    top_genes <- names(sort(gini_scores, decreasing = TRUE)[1:top_n])
    
    marker_list[[as.character(cluster)]] <- data.frame(
      cluster = cluster,
      genes = top_genes,
      gini = gini_scores[top_genes],
      stringsAsFactors = FALSE
    )
  }
  
  # Combine all markers
  all_markers <- do.call(rbind, marker_list)
  rownames(all_markers) <- NULL
  
  return(all_markers)
}

# ============================================================================
# LOAD DATA
# ============================================================================

# Load spatial data (800 spots)
st_obj <- readRDS("data/st_800spots.rds")

# Load scRNA-seq (your uploaded seurat_counts.rds)
sc_obj <- readRDS("code/seurat_counts.rds")

# Load metadata
sc_meta <- read.table("code/sc_mta_data.tsv")
rownames(sc_meta) <- sc_meta$V1
sc_meta <- sc_meta[, -1]

# ============================================================================
# PROCESS scRNA-seq DATA
# ============================================================================

# Intersect genes between scRNA-seq and spatial
inter <- intersect(rownames(sc_obj), rownames(st_obj))
cat("Shared genes:", length(inter), "\n")

sc_obj <- sc_obj[inter, ]
st_obj <- st_obj[inter, ]

# Add cell types to scRNA-seq
sc_obj@meta.data$cell_type <- as.factor(sc_meta$V2)
Idents(sc_obj) <- "cell_type"

# Normalize scRNA-seq data
sc_obj <- NormalizeData(sc_obj, normalization.method = "LogNormalize", 
                        scale.factor = 10000)

# Find marker genes using Gini coefficient (replaces Giotto)
cat("Finding Gini-based marker genes for scRNA-seq...\n")
gini_markers_sc <- find_gini_markers(sc_obj, cluster_column = "cell_type", 
                                     top_n = 100)

topgenes_sc <- gini_markers_sc$genes

# ============================================================================
# PROCESS SPATIAL DATA
# ============================================================================

# Normalize spatial data
st_obj <- NormalizeData(st_obj, normalization.method = "LogNormalize", 
                        scale.factor = 10000)

# Find highly variable genes for spatial clustering
st_obj <- FindVariableFeatures(st_obj, selection.method = "vst", 
                               nfeatures = 2000)

# Scale data
st_obj <- ScaleData(st_obj, features = VariableFeatures(st_obj))

# Run PCA
st_obj <- RunPCA(st_obj, features = VariableFeatures(st_obj), npcs = 10)

# Find neighbors and clusters (replaces Giotto clustering)
st_obj <- FindNeighbors(st_obj, dims = 1:10, k.param = 10)
st_obj <- FindClusters(st_obj, resolution = 0.2, algorithm = 2)  # Leiden

# Rename cluster column to match Giotto naming
st_obj@meta.data$leiden_clus <- st_obj@meta.data$seurat_clusters

# Find marker genes for spatial zones using Gini coefficient
cat("Finding Gini-based marker genes for spatial zones...\n")
gini_markers_st <- find_gini_markers(st_obj, cluster_column = "leiden_clus", 
                                     top_n = 100)

topgenes_st <- gini_markers_st$genes

# ============================================================================
# COMBINE MARKERS
# ============================================================================

all_markers <- unique(c(topgenes_sc, topgenes_st))
cat("Total marker genes:", length(all_markers), "\n")

# ============================================================================
# CALCULATE AVERAGE EXPRESSION PER CELL TYPE
# ============================================================================

# Subset to marker genes
sc_subset <- sc_obj[all_markers, ]

# Calculate average expression (replaces Giotto approach)
avg_expr <- AverageExpression(
  object = sc_subset,
  features = all_markers,
  return.seurat = FALSE,
  group.by = "cell_type",
  slot = "counts"  # Use raw counts for WISpR
)

# ============================================================================
# MATCH GENES BETWEEN scRNA-seq AND SPATIAL
# ============================================================================

st_subset <- st_obj[all_markers, ]
match_idx <- match(rownames(avg_expr$RNA), rownames(st_subset))
st_final <- st_subset[match_idx, ]

# ============================================================================
# SAVE FOR WISpR INPUT
# ============================================================================

# Save scRNA-seq reference
write_csv(
  cbind(rownames(avg_expr$RNA), as.data.frame(avg_expr$RNA)),
  "data/WISpR_sc_reference.csv"
)

# Save spatial data
write_csv(
  cbind(rownames(st_final@assays$RNA@counts),
        as.data.frame(st_final@assays$RNA@counts)),
  "data/WISpR_st_spatial.csv"
)

# Save spatial coordinates for visualization
write_csv(
  as.data.frame(st_obj@meta.data[, c("X_coor", "Y_coor")]),
  "data/spatial_coordinates.csv", 
  row_names = TRUE
)

cat("\nPreprocessing complete!\n")
cat("Files created:\n")
cat(" - WISpR_sc_reference.csv\n")
cat(" - WISpR_st_spatial.csv\n")
cat(" - spatial_coordinates.csv\n")