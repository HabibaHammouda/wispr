"""
DEG Selection Using Gini Coefficient
=====================================

Selects differentially expressed genes (DEGs) for WISpR reference matrix
using Gini coefficient-based approach as described in Erdogan & Eroglu (2025).

The Gini coefficient measures inequality in gene expression across cell types.
High Gini = gene is cell-type-specific (good for deconvolution)
Low Gini = gene is uniformly expressed (housekeeping, not informative)

"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import os

print("="*80)
print("DEG SELECTION USING GINI COEFFICIENT")
print("="*80)
print()

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files (from preprocessing)
SCRNA_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/preprocessed_data/scrna_preprocessed.h5ad"
OUTPUT_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_analysis/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DEG selection parameters
MIN_GINI = 0.3  # Minimum Gini coefficient (genes below this are excluded)
MAX_GENES = 2000  # Maximum number of DEGs to select
MIN_DETECTION_RATE = 0.05  # Gene must be detected in at least 5% of cells per type

print(f"Configuration:")
print(f"  - scRNA data: {SCRNA_FILE}")
print(f"  - Output: {OUTPUT_DIR}")
print(f"  - Min Gini threshold: {MIN_GINI}")
print(f"  - Max DEGs: {MAX_GENES}")
print(f"  - Min detection rate: {MIN_DETECTION_RATE}")
print()

# ============================================================================
# STEP 1: Load Preprocessed scRNA-seq Data
# ============================================================================

print("="*80)
print("STEP 1: LOADING PREPROCESSED DATA")
print("="*80)
print()

print("Loading scRNA-seq reference data...")
adata = sc.read_h5ad(SCRNA_FILE)

print(f"  ✓ Loaded: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"  ✓ Cell types: {adata.obs['cell_type'].nunique()}")
print()

# ============================================================================
# STEP 2: Calculate Mean Expression Per Cell Type
# ============================================================================

print("="*80)
print("STEP 2: CALCULATING MEAN EXPRESSION PER CELL TYPE")
print("="*80)
print()

print("Computing average expression for each gene in each cell type...")

# Get cell types
cell_types = adata.obs['cell_type'].unique()
n_cell_types = len(cell_types)
n_genes = adata.shape[1]

print(f"  Cell types: {n_cell_types}")
print(f"  Genes: {n_genes:,}")
print()

# Create matrix: genes × cell_types
# This will become the reference matrix X for WISpR
mean_expression = np.zeros((n_genes, n_cell_types))
detection_rates = np.zeros((n_genes, n_cell_types))

for i, ct in enumerate(cell_types):
    # Get cells of this type
    mask = adata.obs['cell_type'] == ct
    cells_in_type = adata[mask]
    
    # Calculate mean expression (from normalized counts)
    if sparse.issparse(cells_in_type.X):
        mean_expr = np.array(cells_in_type.X.mean(axis=0)).flatten()
    else:
        mean_expr = cells_in_type.X.mean(axis=0)
    
    mean_expression[:, i] = mean_expr
    
    # Calculate detection rate (% of cells expressing gene > 0)
    if sparse.issparse(cells_in_type.X):
        detection = np.array((cells_in_type.X > 0).sum(axis=0)).flatten() / cells_in_type.shape[0]
    else:
        detection = (cells_in_type.X > 0).sum(axis=0) / cells_in_type.shape[0]
    
    detection_rates[:, i] = detection
    
    print(f"  Processed: {ct} ({mask.sum()} cells)")

print(f"\n  ✓ Created expression matrix: {n_genes:,} genes × {n_cell_types} cell types")
print()

# ============================================================================
# STEP 3: Calculate Gini Coefficient for Each Gene
# ============================================================================

print("="*80)
print("STEP 3: CALCULATING GINI COEFFICIENTS")
print("="*80)
print()

print("Computing Gini coefficient for each gene across cell types...")
print()
print("Gini coefficient interpretation:")
print("  - 0.0: Perfectly uniform expression (housekeeping gene)")
print("  - 1.0: Maximum inequality (highly cell-type-specific)")
print()

def gini_coefficient(x):
    """
    Calculate Gini coefficient for an array.
    
    Gini = (Σ Σ |xi - xj|) / (2n * Σ xi)
    
    Higher values = more inequality = more cell-type-specific
    """
    x = np.array(x).flatten()
    if np.sum(x) == 0:
        return 0.0
    
    # Sort values
    sorted_x = np.sort(x)
    n = len(x)
    index = np.arange(1, n + 1)
    
    # Gini coefficient formula
    gini = (2 * np.sum(index * sorted_x)) / (n * np.sum(sorted_x)) - (n + 1) / n
    
    return gini

# Calculate Gini for each gene
gini_coefficients = np.array([gini_coefficient(mean_expression[i, :]) 
                               for i in range(n_genes)])

print(f"Gini coefficient statistics:")
print(f"  - Mean: {gini_coefficients.mean():.3f}")
print(f"  - Median: {np.median(gini_coefficients):.3f}")
print(f"  - Min: {gini_coefficients.min():.3f}")
print(f"  - Max: {gini_coefficients.max():.3f}")
print()

# ============================================================================
# STEP 4: Filter Genes Based on Gini and Detection Rate
# ============================================================================

print("="*80)
print("STEP 4: FILTERING GENES")
print("="*80)
print()

print("Applying filters to select DEGs...")

# Filter 1: Gini coefficient threshold
high_gini_mask = gini_coefficients >= MIN_GINI
print(f"  Filter 1 (Gini ≥ {MIN_GINI}): {high_gini_mask.sum():,} genes pass")

# Filter 2: Detection rate (gene must be expressed in sufficient cells)
# At least MIN_DETECTION_RATE in at least one cell type
max_detection_per_gene = detection_rates.max(axis=1)
detection_mask = max_detection_per_gene >= MIN_DETECTION_RATE
print(f"  Filter 2 (Detection ≥ {MIN_DETECTION_RATE} in any type): {detection_mask.sum():,} genes pass")

# Combine filters
combined_mask = high_gini_mask & detection_mask
print(f"  Combined filters: {combined_mask.sum():,} genes pass")
print()

# Select top genes by Gini coefficient
candidate_genes = np.where(combined_mask)[0]
candidate_gini = gini_coefficients[candidate_genes]

# Sort by Gini (descending) and take top MAX_GENES
sorted_indices = np.argsort(candidate_gini)[::-1]
top_indices = candidate_genes[sorted_indices[:MAX_GENES]]

# Get gene names
deg_names = adata.var_names[top_indices].tolist()
deg_gini = gini_coefficients[top_indices]

print(f"Selected DEGs: {len(deg_names)}")
print(f"  - Mean Gini: {deg_gini.mean():.3f}")
print(f"  - Min Gini: {deg_gini.min():.3f}")
print(f"  - Max Gini: {deg_gini.max():.3f}")
print()

# ============================================================================
# STEP 5: Build Reference Matrix X
# ============================================================================

print("="*80)
print("STEP 5: BUILDING REFERENCE MATRIX X")
print("="*80)
print()

print("Creating reference matrix for WISpR...")

# Extract mean expression for selected DEGs
reference_matrix = mean_expression[top_indices, :]

print(f"  Reference matrix X: {reference_matrix.shape[0]} genes × {reference_matrix.shape[1]} cell types")
print(f"  Matrix range: [{reference_matrix.min():.3f}, {reference_matrix.max():.3f}]")
print()

# Create DataFrame for reference matrix
reference_df = pd.DataFrame(
    reference_matrix,
    index=deg_names,
    columns=cell_types
)

# Save reference matrix
reference_path = os.path.join(OUTPUT_DIR, "reference_matrix_X.csv")
reference_df.to_csv(reference_path)
print(f"  ✓ Saved reference matrix: {reference_path}")
print()

# Save DEG list with Gini coefficients
deg_info = pd.DataFrame({
    'gene': deg_names,
    'gini_coefficient': deg_gini,
    'max_detection_rate': max_detection_per_gene[top_indices]
})

# Add which cell type has highest expression for each gene
max_expr_celltype_idx = reference_matrix.argmax(axis=1)
deg_info['top_celltype'] = [cell_types[idx] for idx in max_expr_celltype_idx]
deg_info['top_celltype_expression'] = reference_matrix.max(axis=1)

deg_path = os.path.join(OUTPUT_DIR, "selected_degs.csv")
deg_info.to_csv(deg_path, index=False)
print(f"  ✓ Saved DEG information: {deg_path}")
print()

# ============================================================================
# STEP 6: Visualization
# ============================================================================

print("="*80)
print("STEP 6: GENERATING VISUALIZATIONS")
print("="*80)
print()

print("Creating DEG selection visualizations...")

# Figure 1: Gini coefficient distribution
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('DEG Selection Analysis', fontsize=16, fontweight='bold')

# Plot 1: Histogram of all Gini coefficients
axes[0, 0].hist(gini_coefficients, bins=100, edgecolor='black', alpha=0.7)
axes[0, 0].axvline(MIN_GINI, color='red', linestyle='--', linewidth=2, label=f'Threshold ({MIN_GINI})')
axes[0, 0].set_xlabel('Gini Coefficient')
axes[0, 0].set_ylabel('Number of Genes')
axes[0, 0].set_title('Distribution of Gini Coefficients (All Genes)')
axes[0, 0].legend()
axes[0, 0].set_yscale('log')

# Plot 2: Histogram of selected DEG Gini coefficients
axes[0, 1].hist(deg_gini, bins=50, edgecolor='black', alpha=0.7, color='green')
axes[0, 1].set_xlabel('Gini Coefficient')
axes[0, 1].set_ylabel('Number of Genes')
axes[0, 1].set_title(f'Distribution of Selected DEGs (n={len(deg_names)})')

# Plot 3: Detection rate vs Gini
axes[1, 0].scatter(max_detection_per_gene, gini_coefficients, 
                   alpha=0.1, s=1, color='gray', label='All genes')
axes[1, 0].scatter(max_detection_per_gene[top_indices], deg_gini,
                   alpha=0.5, s=10, color='red', label='Selected DEGs')
axes[1, 0].axhline(MIN_GINI, color='red', linestyle='--', alpha=0.5)
axes[1, 0].axvline(MIN_DETECTION_RATE, color='blue', linestyle='--', alpha=0.5)
axes[1, 0].set_xlabel('Max Detection Rate (across cell types)')
axes[1, 0].set_ylabel('Gini Coefficient')
axes[1, 0].set_title('Gene Selection Space')
axes[1, 0].legend()

# Plot 4: Top 20 DEGs by Gini
top_20_idx = np.argsort(deg_gini)[-20:]
top_20_genes = [deg_names[i] for i in top_20_idx]
top_20_gini = deg_gini[top_20_idx]

axes[1, 1].barh(range(20), top_20_gini, color='steelblue', edgecolor='black')
axes[1, 1].set_yticks(range(20))
axes[1, 1].set_yticklabels(top_20_genes, fontsize=8)
axes[1, 1].set_xlabel('Gini Coefficient')
axes[1, 1].set_title('Top 20 DEGs by Cell-Type Specificity')
axes[1, 1].invert_yaxis()

plt.tight_layout()
plot_path = os.path.join(OUTPUT_DIR, "deg_selection_analysis.pdf")
plt.savefig(plot_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved analysis plots: {plot_path}")
plt.close()

# Figure 2: Heatmap of reference matrix (top 50 DEGs)
print("Creating reference matrix heatmap...")

top_50_idx = np.argsort(deg_gini)[-50:]
top_50_matrix = reference_df.iloc[top_50_idx]

plt.figure(figsize=(12, 14))
sns.clustermap(
    top_50_matrix,
    cmap='viridis',
    figsize=(12, 14),
    row_cluster=True,
    col_cluster=True,
    cbar_kws={'label': 'Mean Expression (log-normalized)'},
    yticklabels=True,
    xticklabels=True
)
plt.suptitle('Top 50 DEGs: Expression Across Cell Types', y=1.02, fontsize=16, fontweight='bold')
heatmap_path = os.path.join(OUTPUT_DIR, "reference_matrix_heatmap.pdf")
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved heatmap: {heatmap_path}")
plt.close()

print()

# ============================================================================
# SUMMARY
# ============================================================================

print("="*80)
print("DEG SELECTION SUMMARY")
print("="*80)
print()

print(f"Input Data:")
print(f"  - Total genes: {n_genes:,}")
print(f"  - Cell types: {n_cell_types}")
print(f"  - Cells: {adata.shape[0]:,}")
print()

print(f"Selection Criteria:")
print(f"  - Min Gini coefficient: {MIN_GINI}")
print(f"  - Min detection rate: {MIN_DETECTION_RATE}")
print(f"  - Max DEGs: {MAX_GENES}")
print()

print(f"Results:")
print(f"  - DEGs selected: {len(deg_names)}")
print(f"  - Mean Gini: {deg_gini.mean():.3f}")
print(f"  - Gini range: [{deg_gini.min():.3f}, {deg_gini.max():.3f}]")
print()

print(f"Top 10 Most Cell-Type-Specific Genes:")
top_10_idx = np.argsort(deg_gini)[-10:]
for i, idx in enumerate(top_10_idx[::-1], 1):
    gene = deg_names[idx]
    gini = deg_gini[idx]
    top_ct = deg_info.iloc[idx]['top_celltype']
    print(f"  {i:2d}. {gene:15s} | Gini: {gini:.3f} | Top in: {top_ct}")
print()

print(f"Output Files:")
print(f"  1. {reference_path}")
print(f"     → Reference matrix X (genes × cell_types)")
print(f"  2. {deg_path}")
print(f"     → DEG list with Gini coefficients and metadata")
print(f"  3. {plot_path}")
print(f"     → DEG selection analysis plots")
print(f"  4. {heatmap_path}")
print(f"     → Reference matrix heatmap (top 50 DEGs)")