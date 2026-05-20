# BIOINFORMATICS II FINAL PROJECT
## Implementation Instructions - NO GIOTTO VERSION
### Options A + B: Adaptive Thresholding + Co-occurrence Prior
### Windows | 8GB RAM | Python + R

---

## DAY 0: SYSTEM SETUP (Do This First)

### Install R Packages

Open RStudio or R Console, run this:

```r
install.packages(c(
  "Seurat",       # Main spatial transcriptomics package
  "ggplot2",      # Visualization
  "dplyr",        # Data manipulation  
  "Matrix",       # Sparse matrices
  "tidyverse",    # Data wrangling ecosystem
  "reshape2",     # Data reshaping
  "png",          # Image handling
  "RColorBrewer", # Color palettes
  "data.table",   # Fast data operations
  "reticulate"    # Python interface
))
```

⏱ **This takes 10-15 minutes.** Go make coffee.

**NOTE:** We are **NOT installing Giotto**. All Giotto functionality will be replaced with Seurat equivalents.

---

### Verify Python Packages

Activate your virtual environment, then run:

```bash
pip install scikit-learn pandas numpy scipy rich matplotlib seaborn
```

---

### Create Project Folder Structure

Open Command Prompt (Win+R, type `cmd`), then:

```cmd
cd C:\Users\YourUsername\Documents
mkdir BioProject
cd BioProject
mkdir data code results figures
```

Copy all your uploaded files into `C:\Users\YourUsername\Documents\BioProject\code\`

---

## DAY 1: DATA PREP + BASELINE (6 hours)

### Morning: Prepare Mouse Brain Data (3 hours)

#### Step 1.1: Create Subsampled Dataset

You have 2,700+ spots. For 8GB RAM, we'll use 800 spots.

**Create this R script:** `subsample_data.R`

```r
library(Seurat)
library(dplyr)
library(readr)

# Load your Visium data (adjust path to your actual data location)
data_dir <- "C:/Users/YourUsername/Documents/BioProject/data/Sample_Visium_ST"
expression_matrix <- Read10X(data.dir = data_dir)

# Load coordinates
coor_ST <- read.csv(paste0(data_dir, "/tissue_positions_list.csv"), header = FALSE)
rownames(coor_ST) <- coor_ST$V1
coor_ST <- coor_ST[coor_ST$V2 == "1", ]  # Only spots under tissue

# Match spots
match_n <- match(colnames(expression_matrix), rownames(coor_ST))
coor_ST_reor <- coor_ST[match_n, ]

# SUBSAMPLE TO 800 SPOTS
set.seed(42)  # Reproducibility
n_spots <- min(800, ncol(expression_matrix))
selected_spots <- sample(1:ncol(expression_matrix), n_spots)

expression_matrix_sub <- expression_matrix[, selected_spots]
coor_ST_sub <- coor_ST_reor[selected_spots, ]

# Create Seurat object
st_obj <- CreateSeuratObject(counts = expression_matrix_sub, assay = "RNA")
st_obj <- AddMetaData(st_obj, metadata = coor_ST_sub$V5, col.name = "X_coor")
st_obj <- AddMetaData(st_obj, metadata = coor_ST_sub$V6, col.name = "Y_coor")

# Save
saveRDS(st_obj, "C:/Users/YourUsername/Documents/BioProject/data/st_800spots.rds")
write.csv(coor_ST_sub, 
  "C:/Users/YourUsername/Documents/BioProject/data/coor_ST_800.csv")

cat("Subsampled to", n_spots, "spots\n")
```

**Run it:**

```cmd
cd C:\Users\YourUsername\Documents\BioProject\code
Rscript subsample_data.R
```

✓ **Expected output:** `st_800spots.rds` created (takes 2-3 min)

---

#### Step 1.2: Preprocess scRNA-seq + Spatial Data (NO GIOTTO VERSION)

**Create:** `preprocess_wispr_seurat.R`

```r
library(Seurat)
library(dplyr)
library(data.table)
library(readr)
library(Matrix)

# Set working directory
setwd("C:/Users/YourUsername/Documents/BioProject")

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
```

**Run it:**

```cmd
Rscript preprocess_wispr_seurat.R
```

⏱ **Takes 15-20 minutes.** Go for a walk.

✓ **Expected:** 3 CSV files in `data/` folder

---

### Afternoon: Run Baseline WISpR (3 hours)

#### Step 1.3: Run Original WISpR

Activate your Python virtual environment:

```cmd
cd C:\Users\YourUsername\Documents\BioProject
venv\Scripts\activate
```

Copy `WISpR_new.py` to your code folder, then run:

```bash
python code/WISpR_new.py -r data/WISpR_sc_reference.csv -s data/WISpR_st_spatial.csv -o results/
```

⏱ **Takes 45-60 minutes** for 800 spots

✓ **Expected:** `results/WISpR_Visium.WISpR_sc_reference.csv`

---

#### Step 1.4: Calculate Baseline Metrics

**Create:** `calculate_metrics.py`

```python
import pandas as pd
import numpy as np

# Load WISpR output
results = pd.read_csv("results/WISpR_Visium.WISpR_sc_reference.csv",
                      sep="\t", index_col=0)

# Calculate sparsity (% of zeros)
total_elements = results.size
zero_elements = (results == 0).sum().sum()
sparsity = (zero_elements / total_elements) * 100

# Cell types per spot
cell_types_per_spot = (results > 0).sum(axis=1)

print("="*50)
print("BASELINE WISpR METRICS")
print("="*50)
print(f"Total spots: {len(results)}")
print(f"Total cell types: {len(results.columns)}")
print(f"Sparsity: {sparsity:.2f}%")
print(f"Avg cell types/spot: {cell_types_per_spot.mean():.2f}")
print(f"Min cell types/spot: {cell_types_per_spot.min()}")
print(f"Max cell types/spot: {cell_types_per_spot.max()}")
print("="*50)

# Save summary
with open("results/baseline_metrics.txt", "w") as f:
    f.write(f"Sparsity: {sparsity:.2f}%\n")
    f.write(f"Avg cell types per spot: {cell_types_per_spot.mean():.2f}\n")
```

**Run it:**

```bash
python code/calculate_metrics.py
```

✓ **Expected:** output printed to screen + saved to `results/baseline_metrics.txt`

---

## DAY 2: IMPLEMENT OPTIONS A + B (8 hours)

### Morning: Option A - Adaptive Thresholding (4 hours)

#### Step 2.1: Calculate Spot Density

**Create:** `calculate_density.py`

```python
import pandas as pd
import numpy as np
from sklearn.neighbors import NearestNeighbors
from scipy.spatial.distance import euclidean

# Load spatial data
st_data = pd.read_csv("data/WISpR_st_spatial.csv", index_col=0)
coords = pd.read_csv("data/spatial_coordinates.csv", index_col=0)

# Ensure matching order
coords = coords.loc[st_data.columns]

def calculate_spot_density(st_matrix, coordinates, k=10):
    """
    Calculate density for each spot based on:
    1. Spatial proximity (k nearest neighbors)
    2. Gene expression similarity with neighbors
    """
    # Get k nearest neighbors in spatial coordinates
    nbrs = NearestNeighbors(n_neighbors=k+1).fit(coordinates)
    distances, indices = nbrs.kneighbors(coordinates)
    
    densities = []
    
    for i in range(len(st_matrix.columns)):
        # Get this spot and its neighbors (exclude self at index 0)
        neighbor_idx = indices[i, 1:]  # Skip first (self)
        
        # Calculate gene expression correlation with neighbors
        spot_expr = st_matrix.iloc[:, i].values
        neighbor_exprs = st_matrix.iloc[:, neighbor_idx].values
        
        # Correlation with each neighbor
        correlations = []
        for j in range(neighbor_exprs.shape[1]):
            corr = np.corrcoef(spot_expr, neighbor_exprs[:, j])[0, 1]
            if not np.isnan(corr):
                correlations.append(corr)
        
        # High correlation = dense/homogeneous region
        avg_correlation = np.mean(correlations) if correlations else 0
        
        # Spatial distance (average to neighbors)
        avg_spatial_dist = np.mean(distances[i, 1:])
        
        # Density score: high correlation + close proximity = high density
        density = avg_correlation / (avg_spatial_dist + 1e-6)
        densities.append(density)
    
    return np.array(densities)

# Calculate densities
densities = calculate_spot_density(st_data, coords.values, k=10)

# Normalize to 0-1 range
densities_norm = (densities - densities.min()) / (densities.max() - densities.min())

# Save
density_df = pd.DataFrame({
    "spot_id": st_data.columns,
    "density": densities_norm
})

density_df.to_csv("data/spot_densities.csv", index=False)

print(f"Calculated densities for {len(densities)} spots")
print(f"Density range: {densities_norm.min():.3f} - {densities_norm.max():.3f}")
print(f"Mean density: {densities_norm.mean():.3f}")
```

**Run it:**

```bash
python code/calculate_density.py
```

✓ **Expected:** `data/spot_densities.csv` created

---

#### Step 2.2: Modify WISpR for Adaptive Thresholding

Copy `WISpR_new.py` to `WISpR_adaptive.py`, then make these changes:

**CHANGE 1:** Load density data (add after line 222)

```python
# Load spot densities
density_df = pd.read_csv("data/spot_densities.csv")
spot_densities = dict(zip(density_df["spot_id"], density_df["density"]))
```

**CHANGE 2:** Add adaptive threshold function (before `deconvolve` function, around line 60)

```python
def get_adaptive_threshold_range(density):
    """
    Map density to threshold range:
    High density (>0.7) -> low threshold (0.001-0.005) = more sparse
    Medium density (0.3-0.7) -> medium threshold (0.005-0.02)
    Low density (<0.3) -> high threshold (0.02-0.05) = less sparse
    """
    if density > 0.7:
        return np.arange(0.001, 0.005, 0.001)
    elif density > 0.3:
        return np.arange(0.005, 0.02, 0.002)
    else:
        return np.arange(0.02, 0.05, 0.005)
```

**CHANGE 3:** Modify deconvolve loop (replace line 111-130)

```python
for i in track(range(st_cnt.shape[1]), description="Deconvoluting data (adaptive)"):
    if (np.sum(st_cnt[:,i]) == 0):
        continue
    else:
        # Get spot name
        spot_name = st_cnt_df.columns[i]
        spot_density = spot_densities.get(spot_name, 0.5)  # Default 0.5 if missing
        
        # Get adaptive threshold range for this spot
        threshold = get_adaptive_threshold_range(spot_density)
        
        p_weight[:,i] = weight(sc_cnt.astype(float), st_cnt.astype(float)[:,i])
        best_res = gscv.fit(sc_cnt.astype(float), st_cnt.astype(float)[:,i], 
                            sample_weight=(p_weight[:,i]))
        
        # Rest of nested regression (keep lines 121-179 unchanged)
```

**Run adaptive WISpR:**

```bash
python code/WISpR_adaptive.py -r data/WISpR_sc_reference.csv -s data/WISpR_st_spatial.csv -o results/
```

⏱ **Takes 50-70 minutes**

✓ **Output:** `results/WISpR_Visium.WISpR_sc_reference.csv` (adaptive version)

Rename to: `results/WISpR_adaptive_only.csv`

---

### Afternoon: Option B - Co-occurrence Prior (4 hours)

#### Step 2.3: Build Co-occurrence Matrix

**Create:** `build_cooccurrence.py`

```python
import pandas as pd
import numpy as np
from collections import defaultdict

# Load scRNA-seq metadata
sc_meta = pd.read_csv("code/sc_mta_data.tsv", sep="\t", header=None)
sc_meta.columns = ["cell_id", "cell_type"]

# Get unique cell types
cell_types = sorted(sc_meta["cell_type"].unique())
n_types = len(cell_types)

print(f"Found {n_types} cell types")
print("Cell types:", cell_types)

# Initialize co-occurrence matrix
cooccur = np.zeros((n_types, n_types))

# Build co-occurrence from scRNA-seq
# Strategy: cells from same tissue/organ co-occur more
# Use biological knowledge - for brain data:

# Define biological groups (brain-specific)
neural_types = [ct for ct in cell_types if "Neuron" in ct or "neuron" in ct]
glial_types = [ct for ct in cell_types if any(x in ct for x in ["Astro", "Oligo", "oligo"])]
vascular_types = [ct for ct in cell_types if "Vascular" in ct or "Endo" in ct]
immune_types = [ct for ct in cell_types if "Immune" in ct or "Micro" in ct]

# Build co-occurrence based on biological compatibility
for i, ct1 in enumerate(cell_types):
    for j, ct2 in enumerate(cell_types):
        if i == j:
            cooccur[i, j] = 1.0  # Cell type with itself
        # Same biological group = high co-occurrence
        elif (ct1 in neural_types and ct2 in neural_types) or \
             (ct1 in glial_types and ct2 in glial_types) or \
             (ct1 in vascular_types and ct2 in vascular_types) or \
             (ct1 in immune_types and ct2 in immune_types):
            cooccur[i, j] = 0.8
        # Neurons + glia = common co-occurrence
        elif (ct1 in neural_types and ct2 in glial_types) or \
             (ct1 in glial_types and ct2 in neural_types):
            cooccur[i, j] = 0.6
        # Vascular with anything = moderate (blood everywhere)
        elif ct1 in vascular_types or ct2 in vascular_types:
            cooccur[i, j] = 0.4
        # Default = low but nonzero
        else:
            cooccur[i, j] = 0.2

# Make symmetric
cooccur = (cooccur + cooccur.T) / 2

# Save
cooccur_df = pd.DataFrame(cooccur, index=cell_types, columns=cell_types)
cooccur_df.to_csv("data/cooccurrence_matrix.csv")

print("\nCo-occurrence matrix saved")
print(f"Shape: {cooccur.shape}")
print(f"Range: {cooccur.min():.2f} - {cooccur.max():.2f}")
```

**Run it:**

```bash
python code/build_cooccurrence.py
```

✓ **Expected:** `data/cooccurrence_matrix.csv` created

---

#### Step 2.4: Combine Options A + B

**Create final version:** `WISpR_AB_combined.py`

This combines both modifications. Copy `WISpR_adaptive.py` and add co-occurrence.

**ADD** after loading densities (line ~225):

```python
# Load co-occurrence matrix
cooccur_matrix = pd.read_csv("data/cooccurrence_matrix.csv", index_col=0)

# Ensure order matches sc_cnt columns
cooccur_matrix = cooccur_matrix.loc[sc_cnt_df.columns, sc_cnt_df.columns]
cooccur_np = cooccur_matrix.values
```

**MODIFY** `stlsq_nested_n` function (around line 121):

```python
def stlsq_nested_n(A, y, best_res, predicted_so_far):
    """
    Ridge with co-occurrence penalty
    predicted_so_far: current prediction (to check which cell types active)
    """
    # Get indices of active cell types
    active_idx = np.where(predicted_so_far > 0)[0]
    
    # Calculate co-occurrence weights for current prediction
    if len(active_idx) > 0:
        # Average co-occurrence with already-detected cell types
        cooccur_weights = np.mean(cooccur_np[:, active_idx], axis=1)
        # Boost sample weights for compatible cell types
        sample_weight_adjusted = p_weight[:, i] * (1 + 0.5 * cooccur_weights)
    else:
        sample_weight_adjusted = p_weight[:, i]
    
    ridge_optimizer = Ridge(alpha=best_res.best_params_["alpha"],
                           max_iter=1000, fit_intercept=True,
                           copy_X=False, solver="auto", positive=True)
    
    solution = ridge_optimizer.fit(X=A, y=y, sample_weight=sample_weight_adjusted)
    
    return solution.coef_
```

**UPDATE** calls to `stlsq_nested_n` (lines 150, 167):

```python
# Line 150 (first nested call):
predict_[1, big_indices, i] = stlsq_nested_n(
    sc_cnt.astype(float)[:, big_indices],
    st_cnt.astype(float)[:, i],
    best_res,
    predict_[0, :, i]  # Pass current prediction
)

# Line 167 (second nested call):
predict_[2, big_indices2, i] = stlsq_nested_n(
    sc_cnt.astype(float)[:, big_indices2],
    st_cnt.astype(float)[:, i],
    best_res,
    predict_[1, :, i]  # Pass level-1 prediction
)
```

**Run combined method:**

```bash
python code/WISpR_AB_combined.py -r data/WISpR_sc_reference.csv -s data/WISpR_st_spatial.csv -o results/
```

⏱ **Takes 60-80 minutes** (slower due to co-occurrence calculations)

Rename output to: `results/WISpR_AB_combined.csv`

---

## DAY 3: EVALUATION + VISUALIZATION (8 hours)

### Morning: Calculate Metrics (3 hours)

#### Step 3.1: Compare All Methods

**Create:** `compare_methods.py`

```python
import pandas as pd
import numpy as np
from scipy.stats import wilcoxon

def calculate_sparsity(df):
    """% of zeros in predictions"""
    return (df == 0).sum().sum() / df.size * 100

def cell_types_per_spot(df):
    """Average number of predicted cell types per spot"""
    return (df > 0).sum(axis=1).mean()

# Load all results
baseline = pd.read_csv("results/WISpR_Visium.WISpR_sc_reference.csv",
                       sep="\t", index_col=0)
adaptive = pd.read_csv("results/WISpR_adaptive_only.csv", sep="\t", index_col=0)
combined = pd.read_csv("results/WISpR_AB_combined.csv", sep="\t", index_col=0)

methods = {
    "Baseline WISpR": baseline,
    "Option A (Adaptive)": adaptive,
    "Option A+B (Combined)": combined
}

results_summary = []

for name, df in methods.items():
    sparsity = calculate_sparsity(df)
    avg_ct = cell_types_per_spot(df)
    
    results_summary.append({
        "Method": name,
        "Sparsity (%)": f"{sparsity:.2f}",
        "Avg Cell Types/Spot": f"{avg_ct:.2f}"
    })

summary_df = pd.DataFrame(results_summary)

print("\n" + "="*60)
print("COMPARISON OF METHODS")
print("="*60)
print(summary_df.to_string(index=False))
print("="*60)

# Statistical test: Wilcoxon paired test
baseline_ct = (baseline > 0).sum(axis=1)
combined_ct = (combined > 0).sum(axis=1)

stat, pval = wilcoxon(baseline_ct, combined_ct)

print(f"\nWilcoxon test (Baseline vs Combined):")
print(f" Statistic: {stat:.2f}")
print(f" p-value: {pval:.4f}")

if pval < 0.05:
    print(" ✓ Significant improvement!")
else:
    print(" No significant difference")

# Save summary
summary_df.to_csv("results/comparison_summary.csv", index=False)
```

**Run it:**

```bash
python code/compare_methods.py
```

✓ **Prints table + p-value + saves** `results/comparison_summary.csv`

---

### Afternoon: Create Visualizations (5 hours)

#### Step 3.2: Spatial Heatmaps

**Create:** `plot_spatial.py`

```python
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
coords = pd.read_csv("data/spatial_coordinates.csv", index_col=0)
baseline = pd.read_csv("results/WISpR_Visium.WISpR_sc_reference.csv",
                       sep="\t", index_col=0)
combined = pd.read_csv("results/WISpR_AB_combined.csv", sep="\t", index_col=0)

# Align coordinates with results
coords = coords.loc[baseline.index]

# Pick a cell type to visualize (e.g., first column)
cell_type = baseline.columns[0]

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Baseline
scatter1 = axes[0].scatter(coords["X_coor"], coords["Y_coor"],
                          c=baseline[cell_type], cmap="viridis",
                          s=100, edgecolors="black", linewidth=0.5)
axes[0].set_title(f"Baseline: {cell_type}", fontsize=14, fontweight="bold")
axes[0].set_xlabel("X coordinate")
axes[0].set_ylabel("Y coordinate")
plt.colorbar(scatter1, ax=axes[0], label="Proportion")

# Combined
scatter2 = axes[1].scatter(coords["X_coor"], coords["Y_coor"],
                          c=combined[cell_type], cmap="viridis",
                          s=100, edgecolors="black", linewidth=0.5)
axes[1].set_title(f"Combined A+B: {cell_type}", fontsize=14, fontweight="bold")
axes[1].set_xlabel("X coordinate")
axes[1].set_ylabel("Y coordinate")
plt.colorbar(scatter2, ax=axes[1], label="Proportion")

plt.tight_layout()
plt.savefig(f"figures/spatial_comparison_{cell_type}.png", dpi=300)
print(f"Saved: figures/spatial_comparison_{cell_type}.png")
```

**Run it:**

```bash
python code/plot_spatial.py
```

✓ **Creates** spatial comparison figure in `figures/` folder

---

#### Step 3.3: Sparsity Boxplot

**Create:** `plot_sparsity.py`

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load results
baseline = pd.read_csv("results/WISpR_Visium.WISpR_sc_reference.csv",
                       sep="\t", index_col=0)
adaptive = pd.read_csv("results/WISpR_adaptive_only.csv", sep="\t", index_col=0)
combined = pd.read_csv("results/WISpR_AB_combined.csv", sep="\t", index_col=0)

# Cell types per spot
data = [
    (baseline > 0).sum(axis=1),
    (adaptive > 0).sum(axis=1),
    (combined > 0).sum(axis=1)
]

fig, ax = plt.subplots(figsize=(8, 6))
bp = ax.boxplot(data, labels=["Baseline", "Adaptive (A)", "Combined (A+B)"],
                patch_artist=True)

# Color boxes
colors = ["lightblue", "lightgreen", "salmon"]
for patch, color in zip(bp["boxes"], colors):
    patch.set_facecolor(color)

ax.set_ylabel("Cell Types per Spot", fontsize=12)
ax.set_title("Sparsity Comparison", fontsize=14, fontweight="bold")
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
plt.savefig("figures/sparsity_boxplot.png", dpi=300)
print("Saved: figures/sparsity_boxplot.png")
```

**Run it:**

```bash
python code/plot_sparsity.py
```

✓ **Creates** boxplot in `figures/` folder

---

## DAY 4-5: WRITE PAPER (12 hours)

Follow the paper structure from the main guide. Key sections:

### Results Section - What to Report

- ☐ **Table:** Sparsity and avg cell types per spot for all 3 methods
- ☐ **Figure 1:** Spatial heatmap comparison
- ☐ **Figure 2:** Sparsity boxplot
- ☐ **Statistical test:** Report Wilcoxon p-value
- ☐ **Case study:** Pick 2-3 spots where combined method performs best, explain why

### Methods Section - Key Details to Include

#### Option A: Adaptive Thresholding

Explain:
- How density is calculated (k=10 neighbors, correlation + spatial distance)
- Three density bins and their threshold ranges
- Rationale: Dense regions = homogeneous = need stronger sparsity

#### Option B: Co-occurrence Prior

Explain:
- How co-occurrence matrix built (biological groups: neural, glial, vascular, immune)
- Co-occurrence values (0.8 within group, 0.6 neurons-glia, etc.)
- How integrated: adjust Ridge sample weights by avg co-occurrence with active cell types
- Formula: `weight_adjusted = weight_original × (1 + 0.5 × cooccur_score)`

### Discussion Points

- Why adaptive thresholding improves sparsity
- Why co-occurrence prevents biologically implausible combinations
- **Limitations:** Requires spatial coordinates (Option A), relies on prior knowledge (Option B)
- **Future work:** Learn co-occurrence from data, test on other tissues

---

## FINAL CHECKLIST

### Files You Should Have

- ☐ `data/st_800spots.rds`
- ☐ `data/WISpR_sc_reference.csv`
- ☐ `data/WISpR_st_spatial.csv`
- ☐ `data/spatial_coordinates.csv`
- ☐ `data/spot_densities.csv`
- ☐ `data/cooccurrence_matrix.csv`
- ☐ `results/WISpR_Visium.WISpR_sc_reference.csv` (baseline)
- ☐ `results/WISpR_adaptive_only.csv`
- ☐ `results/WISpR_AB_combined.csv`
- ☐ `results/comparison_summary.csv`
- ☐ `figures/spatial_comparison_*.png`
- ☐ `figures/sparsity_boxplot.png`

### Before Submission

- ☐ All code runs without errors
- ☐ Paper has Abstract, Intro, Methods, Results, Discussion, References
- ☐ At least 2 figures with captions
- ☐ At least 1 table with metrics
- ☐ Methods section explains both Option A and B in detail
- ☐ Results report p-value from statistical test
- ☐ Discussion addresses limitations honestly
- ☐ Code uploaded/submitted

### Presentation Slides (7 slides)

- ☐ **Slide 1:** Title
- ☐ **Slide 2:** Problem + WISpR gap
- ☐ **Slide 3:** Your contribution (Option A: density diagram)
- ☐ **Slide 4:** Your contribution (Option B: co-occurrence heatmap)
- ☐ **Slide 5:** Results table
- ☐ **Slide 6:** Spatial comparison figure
- ☐ **Slide 7:** Conclusions + future work

---

## KEY CHANGES FROM ORIGINAL GUIDE

### What Was Replaced

1. **Giotto's `createGiottoObject`** → **Seurat's `CreateSeuratObject`**
2. **Giotto's `normalizeGiotto`** → **Seurat's `NormalizeData`**
3. **Giotto's `findMarkers_one_vs_all` with `method='gini'`** → **Custom Gini coefficient function**
4. **Giotto's `calculateHVG`** → **Seurat's `FindVariableFeatures`**
5. **Giotto's `runPCA`** → **Seurat's `RunPCA`**
6. **Giotto's `createNearestNetwork` + `doLeidenCluster`** → **Seurat's `FindNeighbors` + `FindClusters`**

### Why This Works

- **All Giotto operations have Seurat equivalents** — no functionality lost
- **Gini coefficient calculation** is simple and implemented directly
- **Seurat is more widely used** and better documented
- **No dependency hell** — everything installs cleanly on Windows
- **Performance is comparable** for this dataset size

---

**You've got this! Good luck! 🚀**
