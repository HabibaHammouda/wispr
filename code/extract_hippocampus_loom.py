"""
Extract Mouse Brain L1 Hippocampus scRNA-seq from Loom File
============================================================

This script processes the Linnarsson Lab Mouse Brain Atlas loom file
(l1_hippocampus.loom) to extract scRNA-seq reference data for WISpR analysis.

Follows the exact preprocessing steps from Erdogan & Eroglu (2025):
- Uses "Clusters" column for Level 1 cell type annotations
- Filters cell types with 25 ≤ cell count ≤ 250

"""

import loompy
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse, io
import os

# Set paths
LOOM_PATH = r"C:\Users\user\Desktop\My_First_BioinfoII_Research_Paper\BioProject\data\l1_hippocampus.loom"
OUTPUT_DIR = r"C:\Users\user\Desktop\My_First_BioinfoII_Research_Paper\BioProject\code\extracted_hippocampus"
os.makedirs(OUTPUT_DIR, exist_ok=True)

print("="*80)
print("EXTRACTING L1 HIPPOCAMPUS scRNA-seq FROM LOOM FILE")
print("="*80)
print()

# ============================================================================
# STEP 1: Load Loom File and Inspect Structure
# ============================================================================

print("Step 1: Loading loom file...")
print(f"  File: {LOOM_PATH}")

# Open loom file in read mode
with loompy.connect(LOOM_PATH, mode='r') as ds:
    print(f"  ✓ Loaded successfully")
    print()
    print(f"Dataset dimensions:")
    print(f"  - Genes (rows): {ds.shape[0]:,}")
    print(f"  - Cells (cols): {ds.shape[1]:,}")
    print()
    
    # ========================================================================
    # STEP 2: Explore Metadata Structure
    # ========================================================================
    
    print("Step 2: Exploring metadata structure...")
    print()
    print("Available row attributes (gene metadata):")
    for attr in ds.ra.keys():
        print(f"  - {attr}")
    print()
    
    print("Available column attributes (cell metadata):")
    for attr in ds.ca.keys():
        print(f"  - {attr}")
    print()
    
    # According to WISpR paper, we need "Class" and "ClusterName"
    # Let's check what's actually available
    
    if 'Class' in ds.ca.keys():
        print("✓ Found 'Class' column (broad cell type categories)")
        print(f"  Unique Classes: {np.unique(ds.ca['Class'])}")
        print()
    
    if 'Clusters' in ds.ca.keys():
        print("✓ Found 'Clusters' column (Level 1 cell type clusters)")
        cluster_ids = ds.ca['Clusters']
        unique_clusters = np.unique(cluster_ids)
        print(f"  Number of unique clusters: {len(unique_clusters)}")
        print(f"  Cluster range: {unique_clusters.min()}-{unique_clusters.max()}")
        print()
    
    # ========================================================================
    # STEP 3: Extract Gene Names
    # ========================================================================
    
    print("Step 3: Extracting gene expression matrix...")
    
    # Get gene names
    # Loom files typically have 'Gene' or 'GeneName' in row attributes
    if 'Gene' in ds.ra.keys():
        gene_names = ds.ra['Gene']
    elif 'GeneName' in ds.ra.keys():
        gene_names = ds.ra['GeneName']
    else:
        # Use first available gene identifier
        gene_attr = list(ds.ra.keys())[0]
        gene_names = ds.ra[gene_attr]
        print(f"  Using '{gene_attr}' as gene identifier")
    
    print(f"  Gene names extracted: {len(gene_names):,}")
    
    # ========================================================================
    # STEP 4: Extract Cell Type Metadata FIRST (before loading matrix)
    # ========================================================================
    
    print()
    print("Step 4: Processing cell type annotations (BEFORE loading matrix)...")
    
    # Get cell barcodes
    if 'CellID' in ds.ca.keys():
        cell_ids = ds.ca['CellID']
    else:
        cell_ids = [f"Cell_{i}" for i in range(ds.shape[1])]
    
    # Extract Class, Subclass, and Clusters
    class_labels = ds.ca['Class'] if 'Class' in ds.ca.keys() else None
    subclass_labels = ds.ca['Subclass'] if 'Subclass' in ds.ca.keys() else None
    cluster_labels = ds.ca['Clusters'] if 'Clusters' in ds.ca.keys() else None
    
    # MEMORY OPTIMIZATION: Also check for 'Excluded' cells and filter them out now
    passed_qc = ds.ca['PassedQC'] if 'PassedQC' in ds.ca.keys() else np.ones(len(cell_ids), dtype=bool)
    
    # Create metadata DataFrame
    metadata = pd.DataFrame({
        'cell_id': cell_ids,
        'Class': class_labels if class_labels is not None else 'Unknown',
        'Subclass': subclass_labels if subclass_labels is not None else 'Unknown',
        'Clusters': cluster_labels if cluster_labels is not None else -1,
        'PassedQC': passed_qc
    })
    
    # Remove excluded cells first
    if class_labels is not None:
        metadata = metadata[metadata['Class'] != 'Excluded'].copy()
        print(f"  Removed 'Excluded' class cells")
    
    # Use Clusters column for cell type identity
    # This is the Level 1 (L1) clustering from Linnarsson Lab
    # Create descriptive labels: Cluster_X_ClassName_Subclass
    if cluster_labels is not None:
        metadata['cell_type'] = 'Cluster_' + metadata['Clusters'].astype(str) + '_' + \
                                 metadata['Class'].astype(str)
        # Also keep the simple cluster ID for reference
        metadata['cluster_id'] = metadata['Clusters']
    else:
        raise ValueError("No Clusters column found in loom file!")
    
    print(f"  Cell types (clusters) identified: {metadata['cell_type'].nunique()}")
    print(f"  Cluster range: {metadata['cluster_id'].min()}-{metadata['cluster_id'].max()}")
    print()
    
    # ========================================================================
    # STEP 5: Apply Cell Type Filtering (25-250 cells per type) BEFORE LOADING MATRIX
    # ========================================================================
    
    print("Step 5: Filtering cell types (25 ≤ cells ≤ 250) BEFORE loading matrix...")
    print()
    
    # Count cells per type
    cell_type_counts = metadata['cell_type'].value_counts()
    print("Cell type distribution BEFORE filtering:")
    print(cell_type_counts.head(20))
    print(f"  Total cell types: {len(cell_type_counts)}")
    print()
    
    # Apply filtering criteria from Andersson et al. 2020
    valid_cell_types = cell_type_counts[(cell_type_counts >= 25) & (cell_type_counts <= 250)].index
    
    print(f"Cell types passing filter: {len(valid_cell_types)}")
    print(f"Cell types REMOVED (< 25 cells): {((cell_type_counts < 25).sum())}")
    print(f"Cell types REMOVED (> 250 cells): {((cell_type_counts > 250).sum())}")
    print()
    
    # Filter metadata to only valid cell types
    metadata_filtered = metadata[metadata['cell_type'].isin(valid_cell_types)].copy()
    
    # Get indices of cells to keep (in original matrix)
    cell_indices_to_keep = np.where(metadata['cell_type'].isin(valid_cell_types))[0]
    
    print(f"Cells to load: {len(cell_indices_to_keep):,} (from original {len(metadata):,})")
    print(f"This reduces memory usage by ~{100 * (1 - len(cell_indices_to_keep)/len(metadata)):.1f}%")
    print()
    
    # ========================================================================
    # STEP 6: Load ONLY filtered cells from matrix (MEMORY EFFICIENT)
    # ========================================================================
    
    print("Step 6: Loading expression matrix for FILTERED cells only...")
    
    # Load only the columns we need
    counts_filtered = ds[:, cell_indices_to_keep]
    
    print(f"  ✓ Matrix loaded: {counts_filtered.shape[0]:,} genes × {counts_filtered.shape[1]:,} cells")
    print(f"  Matrix sparsity: {(counts_filtered == 0).sum() / counts_filtered.size * 100:.1f}% zeros")
    print()
    
    # Show final cell type distribution
    final_counts = metadata_filtered['cell_type'].value_counts()
    print("Final cell type distribution:")
    print(final_counts)
    print()

# ============================================================================
# STEP 7: Quality Control - Remove Low-Quality Genes
# ============================================================================

print("Step 7: Gene-level quality control...")

# Remove genes with zero expression across all cells
gene_sums = counts_filtered.sum(axis=1)
genes_to_keep = gene_sums > 0

gene_names_filtered = gene_names[genes_to_keep]
counts_filtered = counts_filtered[genes_to_keep, :]

print(f"  Removed {(~genes_to_keep).sum():,} genes with zero expression")
print(f"  Retained genes: {counts_filtered.shape[0]:,}")
print()

# ============================================================================
# STEP 8: Save as Scanpy AnnData Object (Python-friendly format)
# ============================================================================

print("Step 8: Creating Scanpy AnnData object...")

# Create AnnData object
adata = sc.AnnData(
    X=sparse.csr_matrix(counts_filtered.T),  # Transpose to cells × genes for Scanpy
    obs=metadata_filtered.reset_index(drop=True),  # Cell metadata
    var=pd.DataFrame({'gene_name': gene_names_filtered})  # Gene metadata
)

# Add raw counts layer
adata.raw = adata.copy()

print(f"  AnnData object created:")
print(f"    - Shape: {adata.shape[0]} cells × {adata.shape[1]} genes")
print(f"    - Cell types: {adata.obs['cell_type'].nunique()}")
print()

# Save as h5ad
adata_path = os.path.join(OUTPUT_DIR, "hippocampus_scrna_filtered.h5ad")
adata.write_h5ad(adata_path)
print(f"  ✓ Saved AnnData: {adata_path}")
print()

# ============================================================================
# STEP 9: Export to 10X Format (for R/Seurat compatibility)
# ============================================================================

print("Step 9: Exporting to 10X MTX format for Seurat...")

mtx_dir = os.path.join(OUTPUT_DIR, "hippocampus_scrna_10x/")
os.makedirs(mtx_dir, exist_ok=True)

# Save matrix
io.mmwrite(
    os.path.join(mtx_dir, "matrix.mtx"),
    sparse.csr_matrix(counts_filtered)
)

# Save gene names (features.tsv format: gene_id, gene_name, feature_type)
features_df = pd.DataFrame({
    'gene_id': gene_names_filtered,
    'gene_name': gene_names_filtered,
    'feature_type': ['Gene Expression'] * len(gene_names_filtered)
})
features_df.to_csv(
    os.path.join(mtx_dir, "features.tsv.gz"),
    sep='\t',
    header=False,
    index=False,
    compression='gzip'
)

# Save cell barcodes
metadata_filtered[['cell_id']].to_csv(
    os.path.join(mtx_dir, "barcodes.tsv.gz"),
    header=False,
    index=False,
    compression='gzip'
)

print(f"  ✓ Saved 10X format files:")
print(f"    - {mtx_dir}matrix.mtx")
print(f"    - {mtx_dir}features.tsv.gz")
print(f"    - {mtx_dir}barcodes.tsv.gz")
print()

# ============================================================================
# STEP 10: Save Cell Type Annotations Separately
# ============================================================================

print("Step 10: Saving cell type metadata...")

metadata_path = os.path.join(OUTPUT_DIR, "hippocampus_cell_metadata.tsv")
metadata_filtered.to_csv(metadata_path, sep='\t', index=False)

print(f"  ✓ Saved metadata: {metadata_path}")
print()

# ============================================================================
# STEP 11: Generate Summary Statistics
# ============================================================================

print("="*80)
print("EXTRACTION SUMMARY")
print("="*80)
print()

print(f"Original dataset:")
print(f"  - Genes: {len(gene_names):,}")
print(f"  - Cells (before metadata filter): {len(metadata):,}")
print(f"  - Cell types (before filter): {len(cell_type_counts)}")
print()

print(f"Filtered dataset (ready for WISpR):")
print(f"  - Genes: {counts_filtered.shape[0]:,}")
print(f"  - Cells: {counts_filtered.shape[1]:,}")
print(f"  - Cell types: {len(valid_cell_types)}")
print()

print(f"Cell type filter criteria (from WISpR paper):")
print(f"  - Minimum cells per type: 25")
print(f"  - Maximum cells per type: 250")
print()

print(f"Output files created:")
print(f"  1. {adata_path}")
print(f"     → Python/Scanpy AnnData object")
print(f"  2. {mtx_dir}")
print(f"     → 10X format for R/Seurat (matrix.mtx, features.tsv.gz, barcodes.tsv.gz)")
print(f"  3. {metadata_path}")
print(f"     → Cell type annotations (TSV)")