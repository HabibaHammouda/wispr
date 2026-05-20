"""
WISpR Preprocessing Pipeline - Python Version
==============================================

Loads both spatial transcriptomics (Visium) and scRNA-seq reference data,
performs QC, normalization, and gene matching.

This replaces the R/Seurat workflow with pure Python (Scanpy/Squidpy).
"""

import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import sparse
import os
import warnings
warnings.filterwarnings('ignore')

# Set random seed for reproducibility
np.random.seed(42)

# Configure scanpy
sc.settings.verbosity = 3  # Verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, facecolor='white')

print("="*80)
print("WISpR PREPROCESSING PIPELINE - PYTHON VERSION")
print("="*80)
print()

# ============================================================================
# CONFIGURATION: UPDATE THESE PATHS FOR YOUR SYSTEM
# ============================================================================

# INPUT PATHS
SPATIAL_DATA_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/filtered_feature_bc_matrix/"
SPATIAL_COORDS_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/tissue_positions_list.csv"
SPATIAL_IMAGE_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/"

SCRNA_DATA_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/extracted_hippocampus/hippocampus_scrna_10x/"
SCRNA_METADATA_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/data/extracted_hippocampus/hippocampus_cell_metadata.tsv"

# OUTPUT DIRECTORY
OUTPUT_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/preprocessed_data/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"Configuration:")
print(f"  - Spatial data: {SPATIAL_DATA_DIR}")
print(f"  - scRNA data: {SCRNA_DATA_DIR}")
print(f"  - Output: {OUTPUT_DIR}")
print()

# ============================================================================
# PART 1: Load Spatial Transcriptomics Data (10X Visium)
# ============================================================================

print("="*80)
print("PART 1: LOADING SPATIAL TRANSCRIPTOMICS DATA")
print("="*80)
print()

print("Loading Visium spatial data from 10X format...")

# Read 10X format data
adata_spatial = sc.read_10x_mtx(
    SPATIAL_DATA_DIR,
    var_names='gene_symbols',  # Use gene symbols as variable names
    cache=True
)

print(f"  ✓ Loaded matrix: {adata_spatial.shape[0]} spots × {adata_spatial.shape[1]} genes")

# Load spatial coordinates
spatial_coords = pd.read_csv(SPATIAL_COORDS_FILE, header=None)
spatial_coords.columns = ['barcode', 'in_tissue', 'array_row', 'array_col', 
                          'pxl_row_in_fullres', 'pxl_col_in_fullres']

print(f"  ✓ Loaded coordinates for {len(spatial_coords)} spots")

# Filter to in-tissue spots only
in_tissue_barcodes = spatial_coords[spatial_coords['in_tissue'] == 1]['barcode'].values
adata_spatial = adata_spatial[adata_spatial.obs_names.isin(in_tissue_barcodes)].copy()

print(f"  ✓ Filtered to {adata_spatial.shape[0]} in-tissue spots")

# Add spatial coordinates to adata object
coords_matched = spatial_coords.set_index('barcode').loc[adata_spatial.obs_names]
adata_spatial.obs['array_row'] = coords_matched['array_row'].values
adata_spatial.obs['array_col'] = coords_matched['array_col'].values
adata_spatial.obs['pxl_row'] = coords_matched['pxl_row_in_fullres'].values
adata_spatial.obs['pxl_col'] = coords_matched['pxl_col_in_fullres'].values

# Add spatial coordinates for Scanpy spatial plotting
adata_spatial.obsm['spatial'] = coords_matched[['pxl_col_in_fullres', 'pxl_row_in_fullres']].values

print(f"  ✓ Added spatial coordinates")
print()

# ============================================================================
# PART 2: Load scRNA-seq Reference Data
# ============================================================================

print("="*80)
print("PART 2: LOADING scRNA-seq REFERENCE DATA")
print("="*80)
print()

print("Loading hippocampus scRNA-seq reference...")

# Manually read the files since matrix.mtx is uncompressed
from scipy.io import mmread
import gzip

# Read matrix
matrix_path = os.path.join(SCRNA_DATA_DIR, "matrix.mtx")
X = mmread(matrix_path).T.tocsr()  # Transpose to cells × genes and convert to CSR
print(f"  ✓ Loaded expression matrix: {X.shape}")

# Read barcodes (cell IDs)
barcodes_path = os.path.join(SCRNA_DATA_DIR, "barcodes.tsv.gz")
with gzip.open(barcodes_path, 'rt') as f:
    barcodes = [line.strip() for line in f]
print(f"  ✓ Loaded {len(barcodes)} cell barcodes")

# Read features (genes)
features_path = os.path.join(SCRNA_DATA_DIR, "features.tsv.gz")
with gzip.open(features_path, 'rt') as f:
    features = [line.strip().split('\t') for line in f]
gene_ids = [f[0] for f in features]
gene_names = [f[1] if len(f) > 1 else f[0] for f in features]
print(f"  ✓ Loaded {len(gene_names)} genes")

# Create AnnData object
adata_scrna = sc.AnnData(
    X=X,
    obs=pd.DataFrame(index=barcodes),
    var=pd.DataFrame({'gene_ids': gene_ids}, index=gene_names)
)

print(f"  ✓ Created AnnData: {adata_scrna.shape[0]} cells × {adata_scrna.shape[1]} genes")

# Load cell type annotations
scrna_metadata = pd.read_csv(SCRNA_METADATA_FILE, sep='\t', index_col=0)

# Add metadata to adata
adata_scrna.obs['cell_type'] = scrna_metadata['cell_type'].values
adata_scrna.obs['broad_class'] = scrna_metadata['Class'].values
adata_scrna.obs['cluster_id'] = scrna_metadata['cluster_id'].values

print(f"  ✓ Added cell type annotations")
print(f"  ✓ Cell types: {adata_scrna.obs['cell_type'].nunique()}")
print()

# ============================================================================
# PART 3: Quality Control - Spatial Data
# ============================================================================

print("="*80)
print("PART 3: QUALITY CONTROL - SPATIAL DATA")
print("="*80)
print()

print("Calculating QC metrics for spatial data...")

# Calculate QC metrics
adata_spatial.var['mt'] = adata_spatial.var_names.str.startswith('mt-')  # Mitochondrial genes
sc.pp.calculate_qc_metrics(adata_spatial, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"  Spatial QC metrics:")
print(f"    - Total counts per spot: {adata_spatial.obs['total_counts'].median():.0f} (median)")
print(f"    - Genes per spot: {adata_spatial.obs['n_genes_by_counts'].median():.0f} (median)")
print(f"    - Mitochondrial %: {adata_spatial.obs['pct_counts_mt'].median():.2f}% (median)")
print()

# QC filtering
print("Applying QC filters to spatial data...")
print(f"  Before filtering: {adata_spatial.shape[0]} spots")

sc.pp.filter_cells(adata_spatial, min_genes=200)  # At least 200 genes
sc.pp.filter_cells(adata_spatial, max_genes=8000)  # Remove potential doublets
adata_spatial = adata_spatial[adata_spatial.obs['pct_counts_mt'] < 20, :]  # <20% mitochondrial

print(f"  After filtering: {adata_spatial.shape[0]} spots")
print(f"    - Removed {adata_spatial.shape[0]} spots with QC issues")
print()

# ============================================================================
# PART 4: Quality Control - scRNA-seq Data
# ============================================================================

print("="*80)
print("PART 4: QUALITY CONTROL - scRNA-seq DATA")
print("="*80)
print()

print("Calculating QC metrics for scRNA-seq data...")

# Calculate QC metrics
adata_scrna.var['mt'] = adata_scrna.var_names.str.startswith('mt-')
sc.pp.calculate_qc_metrics(adata_scrna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

print(f"  scRNA QC metrics:")
print(f"    - Total counts per cell: {adata_scrna.obs['total_counts'].median():.0f} (median)")
print(f"    - Genes per cell: {adata_scrna.obs['n_genes_by_counts'].median():.0f} (median)")
print(f"    - Mitochondrial %: {adata_scrna.obs['pct_counts_mt'].median():.2f}% (median)")
print()

# QC filtering (minimal since data was already filtered during extraction)
print("Applying QC filters to scRNA-seq data...")
print(f"  Before filtering: {adata_scrna.shape[0]} cells")

sc.pp.filter_cells(adata_scrna, min_genes=100)  # At least 100 genes
adata_scrna = adata_scrna[adata_scrna.obs['pct_counts_mt'] < 25, :]  # <25% mitochondrial

print(f"  After filtering: {adata_scrna.shape[0]} cells")
print()

# ============================================================================
# PART 5: Match Genes Between Datasets
# ============================================================================

print("="*80)
print("PART 5: GENE MATCHING")
print("="*80)
print()

print("Finding common genes between spatial and scRNA-seq datasets...")

# First, make gene names unique by adding suffixes to duplicates
print("Making gene names unique...")
adata_spatial.var_names_make_unique()
adata_scrna.var_names_make_unique()
print(f"  ✓ Gene names made unique")

# Get gene lists
spatial_genes = set(adata_spatial.var_names)
scrna_genes = set(adata_scrna.var_names)

# Find intersection
common_genes = sorted(list(spatial_genes.intersection(scrna_genes)))

print(f"  Spatial genes: {len(spatial_genes):,}")
print(f"  scRNA genes: {len(scrna_genes):,}")
print(f"  Common genes: {len(common_genes):,}")
print(f"  Overlap: {100 * len(common_genes) / max(len(spatial_genes), len(scrna_genes)):.1f}%")
print()

if len(common_genes) < 1000:
    print("  ⚠️  WARNING: Very few common genes! Check gene naming consistency.")
    print("  Spatial genes sample:", list(spatial_genes)[:5])
    print("  scRNA genes sample:", list(scrna_genes)[:5])
else:
    print(f"  ✓ Sufficient gene overlap for WISpR analysis")

# Subset both datasets to common genes
print("Subsetting both datasets to common genes...")
adata_spatial = adata_spatial[:, common_genes].copy()
adata_scrna = adata_scrna[:, common_genes].copy()

print(f"  ✓ Spatial data: {adata_spatial.shape[0]} spots × {adata_spatial.shape[1]} genes")
print(f"  ✓ scRNA data: {adata_scrna.shape[0]} cells × {adata_scrna.shape[1]} genes")
print()

# ============================================================================
# PART 6: Normalization - Spatial Data
# ============================================================================

print("="*80)
print("PART 6: NORMALIZATION - SPATIAL DATA")
print("="*80)
print()

print("Normalizing spatial data...")

# Store raw counts
adata_spatial.layers['counts'] = adata_spatial.X.copy()

# Normalize to 10,000 counts per spot (standard for Visium)
sc.pp.normalize_total(adata_spatial, target_sum=1e4)

# Log transform
sc.pp.log1p(adata_spatial)

print(f"  ✓ Normalized and log-transformed")
print()

# ============================================================================
# PART 7: Normalization - scRNA-seq Data
# ============================================================================

print("="*80)
print("PART 7: NORMALIZATION - scRNA-seq DATA")
print("="*80)
print()

print("Normalizing scRNA-seq data...")

# Store raw counts
adata_scrna.layers['counts'] = adata_scrna.X.copy()

# Normalize to 10,000 counts per cell
sc.pp.normalize_total(adata_scrna, target_sum=1e4)

# Log transform
sc.pp.log1p(adata_scrna)

print(f"  ✓ Normalized and log-transformed")
print()

# ============================================================================
# PART 8: Feature Selection
# ============================================================================

print("="*80)
print("PART 8: FEATURE SELECTION")
print("="*80)
print()

print("Identifying highly variable genes in spatial data...")
sc.pp.highly_variable_genes(adata_spatial, n_top_genes=2000, flavor='seurat_v3', layer='counts')
print(f"  ✓ {adata_spatial.var['highly_variable'].sum()} highly variable genes identified")
print()

print("Identifying highly variable genes in scRNA-seq data...")
sc.pp.highly_variable_genes(adata_scrna, n_top_genes=2000, flavor='seurat_v3', layer='counts')
print(f"  ✓ {adata_scrna.var['highly_variable'].sum()} highly variable genes identified")
print()

# ============================================================================
# PART 9: Save Preprocessed Data
# ============================================================================

print("="*80)
print("PART 9: SAVING PREPROCESSED DATA")
print("="*80)
print()

print("Saving preprocessed objects...")

# Save spatial data
spatial_output = os.path.join(OUTPUT_DIR, "spatial_preprocessed.h5ad")
adata_spatial.write_h5ad(spatial_output)
print(f"  ✓ Saved spatial data: {spatial_output}")

# Save scRNA data
scrna_output = os.path.join(OUTPUT_DIR, "scrna_preprocessed.h5ad")
adata_scrna.write_h5ad(scrna_output)
print(f"  ✓ Saved scRNA data: {scrna_output}")
print()

# ============================================================================
# PART 10: Generate QC Plots
# ============================================================================

print("="*80)
print("PART 10: GENERATING QC PLOTS")
print("="*80)
print()

print("Creating QC visualizations...")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
fig.suptitle('Quality Control Metrics', fontsize=16, fontweight='bold')

# Spatial QC plots
axes[0, 0].hist(adata_spatial.obs['total_counts'], bins=50, edgecolor='black')
axes[0, 0].set_xlabel('Total UMI counts')
axes[0, 0].set_ylabel('Number of spots')
axes[0, 0].set_title('Spatial: Total Counts Distribution')
axes[0, 0].axvline(adata_spatial.obs['total_counts'].median(), color='red', linestyle='--', label='Median')
axes[0, 0].legend()

axes[0, 1].hist(adata_spatial.obs['n_genes_by_counts'], bins=50, edgecolor='black')
axes[0, 1].set_xlabel('Number of genes')
axes[0, 1].set_ylabel('Number of spots')
axes[0, 1].set_title('Spatial: Genes per Spot Distribution')
axes[0, 1].axvline(adata_spatial.obs['n_genes_by_counts'].median(), color='red', linestyle='--', label='Median')
axes[0, 1].legend()

axes[0, 2].hist(adata_spatial.obs['pct_counts_mt'], bins=50, edgecolor='black')
axes[0, 2].set_xlabel('Mitochondrial %')
axes[0, 2].set_ylabel('Number of spots')
axes[0, 2].set_title('Spatial: Mitochondrial Percentage')
axes[0, 2].axvline(adata_spatial.obs['pct_counts_mt'].median(), color='red', linestyle='--', label='Median')
axes[0, 2].legend()

# scRNA QC plots
axes[1, 0].hist(adata_scrna.obs['total_counts'], bins=50, edgecolor='black')
axes[1, 0].set_xlabel('Total UMI counts')
axes[1, 0].set_ylabel('Number of cells')
axes[1, 0].set_title('scRNA: Total Counts Distribution')
axes[1, 0].axvline(adata_scrna.obs['total_counts'].median(), color='red', linestyle='--', label='Median')
axes[1, 0].legend()

axes[1, 1].hist(adata_scrna.obs['n_genes_by_counts'], bins=50, edgecolor='black')
axes[1, 1].set_xlabel('Number of genes')
axes[1, 1].set_ylabel('Number of cells')
axes[1, 1].set_title('scRNA: Genes per Cell Distribution')
axes[1, 1].axvline(adata_scrna.obs['n_genes_by_counts'].median(), color='red', linestyle='--', label='Median')
axes[1, 1].legend()

axes[1, 2].hist(adata_scrna.obs['pct_counts_mt'], bins=50, edgecolor='black')
axes[1, 2].set_xlabel('Mitochondrial %')
axes[1, 2].set_ylabel('Number of cells')
axes[1, 2].set_title('scRNA: Mitochondrial Percentage')
axes[1, 2].axvline(adata_scrna.obs['pct_counts_mt'].median(), color='red', linestyle='--', label='Median')
axes[1, 2].legend()

plt.tight_layout()
qc_plot_path = os.path.join(OUTPUT_DIR, "qc_metrics.pdf")
plt.savefig(qc_plot_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved QC plots: {qc_plot_path}")
plt.close()

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print()
print("="*80)
print("PREPROCESSING SUMMARY")
print("="*80)
print()

print("Spatial Transcriptomics Data:")
print(f"  - Spots (after QC): {adata_spatial.shape[0]:,}")
print(f"  - Genes: {adata_spatial.shape[1]:,}")
print(f"  - Mean UMI/spot: {adata_spatial.obs['total_counts'].mean():.1f}")
print(f"  - Mean genes/spot: {adata_spatial.obs['n_genes_by_counts'].mean():.1f}")
print(f"  - Highly variable genes: {adata_spatial.var['highly_variable'].sum()}")
print()

print("scRNA-seq Reference Data:")
print(f"  - Cells (after QC): {adata_scrna.shape[0]:,}")
print(f"  - Genes: {adata_scrna.shape[1]:,}")
print(f"  - Cell types: {adata_scrna.obs['cell_type'].nunique()}")
print(f"  - Mean UMI/cell: {adata_scrna.obs['total_counts'].mean():.1f}")
print(f"  - Mean genes/cell: {adata_scrna.obs['n_genes_by_counts'].mean():.1f}")
print(f"  - Highly variable genes: {adata_scrna.var['highly_variable'].sum()}")
print()

print("Cell type distribution:")
cell_type_counts = adata_scrna.obs['cell_type'].value_counts()
for ct, count in cell_type_counts.head(10).items():
    print(f"  - {ct}: {count} cells")
if len(cell_type_counts) > 10:
    print(f"  ... and {len(cell_type_counts) - 10} more cell types")
print()

print("Output Files Created:")
print(f"  1. {spatial_output}")
print(f"     → Preprocessed spatial data (H5AD format)")
print(f"  2. {scrna_output}")
print(f"     → Preprocessed scRNA reference (H5AD format)")
print(f"  3. {qc_plot_path}")
print(f"     → QC metrics visualization")
print()
print("✓ Preprocessing pipeline complete!")