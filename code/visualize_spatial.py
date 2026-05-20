#!/usr/bin/env python3
"""
Spatial Visualization of WISpR Results
======================================

Creates spatial maps showing where each cell type is located in the tissue.
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
import os

print("="*80)
print("SPATIAL VISUALIZATION OF WISpR RESULTS")
print("="*80)
print()

# ============================================================================
# CONFIGURATION
# ============================================================================

RESULTS_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_results/spatial_deconvolved.h5ad"
OUTPUT_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_results/"

print(f"Loading results from: {RESULTS_FILE}")
print()

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading deconvolved spatial data...")
adata = sc.read_h5ad(RESULTS_FILE)

print(f"  ✓ Loaded: {adata.shape[0]} spots")
print(f"  ✓ Spatial coordinates available: {adata.obsm['spatial'].shape}")
print()

# Get cell type columns
cell_type_cols = [col for col in adata.obs.columns if col.startswith('prop_')]
cell_types = [col.replace('prop_', '') for col in cell_type_cols]
print(f"  ✓ Cell types: {len(cell_types)}")
print()

# ============================================================================
# IDENTIFY TOP CELL TYPES BY PREVALENCE
# ============================================================================

print("Calculating cell type prevalence...")

# Calculate how many spots each cell type appears in (proportion > 1%)
prevalence = {}
mean_proportion = {}

for ct in cell_types:
    props = adata.obs[f'prop_{ct}'].values
    prevalence[ct] = (props > 0.01).sum()
    mean_proportion[ct] = props.mean()

# Sort by prevalence
sorted_celltypes = sorted(cell_types, key=lambda x: prevalence[x], reverse=True)

print(f"\nTop 10 most prevalent cell types:")
for i, ct in enumerate(sorted_celltypes[:10], 1):
    print(f"  {i:2d}. {ct:30s}: {prevalence[ct]:4d} spots ({100*prevalence[ct]/adata.shape[0]:5.1f}%)")
print()

# ============================================================================
# FIGURE 1: SPATIAL MAPS OF TOP 9 CELL TYPES
# ============================================================================

print("Creating spatial maps for top 9 cell types...")

top_9_celltypes = sorted_celltypes[:9]

fig, axes = plt.subplots(3, 3, figsize=(18, 18))
axes = axes.flatten()

for i, ct in enumerate(top_9_celltypes):
    ax = axes[i]
    
    # Get proportions
    props = adata.obs[f'prop_{ct}'].values
    
    # Get spatial coordinates
    coords = adata.obsm['spatial']
    
    # Create scatter plot
    scatter = ax.scatter(
        coords[:, 0], 
        coords[:, 1],
        c=props,
        cmap='viridis',
        s=8,
        vmin=0,
        vmax=np.percentile(props[props > 0], 95) if (props > 0).sum() > 0 else 1  # Cap at 95th percentile for better contrast
    )
    
    # Format title
    n_detected = (props > 0.01).sum()
    mean_prop = props[props > 0.01].mean() if n_detected > 0 else 0
    
    ax.set_title(
        f'{ct}\n{n_detected} spots ({100*n_detected/adata.shape[0]:.1f}%) | Mean: {mean_prop:.3f}',
        fontsize=9
    )
    ax.axis('equal')
    ax.invert_yaxis()
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Proportion', fontsize=8)
    cbar.ax.tick_params(labelsize=7)

plt.suptitle('Spatial Distribution: Top 9 Cell Types', fontsize=16, fontweight='bold')
plt.tight_layout()

output_path_top9 = os.path.join(OUTPUT_DIR, "spatial_maps_top9.pdf")
plt.savefig(output_path_top9, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_path_top9}")
plt.close()

# ============================================================================
# FIGURE 2: COMPOSITE VIEW - DOMINANT CELL TYPE PER SPOT
# ============================================================================

print("Creating dominant cell type map...")

# For each spot, find the cell type with highest proportion
dominant_celltype_idx = np.zeros(adata.shape[0], dtype=int)
dominant_proportion = np.zeros(adata.shape[0])

for i in range(adata.shape[0]):
    props = np.array([adata.obs[f'prop_{ct}'].values[i] for ct in cell_types])
    dominant_celltype_idx[i] = np.argmax(props)
    dominant_proportion[i] = props[dominant_celltype_idx[i]]

# Map to cell type names
dominant_celltype_names = [cell_types[idx] for idx in dominant_celltype_idx]

# Create color map
unique_celltypes = list(set(dominant_celltype_names))
color_map = dict(zip(unique_celltypes, sns.color_palette('tab20', len(unique_celltypes))))

# Plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Dominant cell type (categorical)
coords = adata.obsm['spatial']
for ct in unique_celltypes:
    mask = np.array(dominant_celltype_names) == ct
    ax1.scatter(
        coords[mask, 0],
        coords[mask, 1],
        c=[color_map[ct]],
        s=8,
        label=ct if prevalence[ct] > 50 else None,  # Only label common types
        alpha=0.8
    )

ax1.set_title('Dominant Cell Type per Spot', fontsize=14, fontweight='bold')
ax1.axis('equal')
ax1.invert_yaxis()
ax1.set_xticks([])
ax1.set_yticks([])
ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8, ncol=1)

# Right: Dominant proportion (confidence)
scatter = ax2.scatter(
    coords[:, 0],
    coords[:, 1],
    c=dominant_proportion,
    cmap='Reds',
    s=8,
    vmin=0,
    vmax=1
)
ax2.set_title('Confidence (Proportion of Dominant Type)', fontsize=14, fontweight='bold')
ax2.axis('equal')
ax2.invert_yaxis()
ax2.set_xticks([])
ax2.set_yticks([])
plt.colorbar(scatter, ax=ax2, label='Proportion', fraction=0.046, pad=0.04)

plt.tight_layout()
output_path_dominant = os.path.join(OUTPUT_DIR, "spatial_dominant_celltype.pdf")
plt.savefig(output_path_dominant, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_path_dominant}")
plt.close()

# ============================================================================
# FIGURE 3: SPARSITY MAP
# ============================================================================

print("Creating sparsity visualization...")

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left: Number of cell types per spot
n_celltypes = adata.obs['n_celltypes'].values
scatter1 = ax1.scatter(
    coords[:, 0],
    coords[:, 1],
    c=n_celltypes,
    cmap='plasma',
    s=8,
    vmin=n_celltypes.min(),
    vmax=n_celltypes.max()
)
ax1.set_title(f'Sparsity: Number of Cell Types per Spot\n(Mean: {n_celltypes.mean():.1f}, Median: {np.median(n_celltypes):.0f})', 
              fontsize=12, fontweight='bold')
ax1.axis('equal')
ax1.invert_yaxis()
ax1.set_xticks([])
ax1.set_yticks([])
plt.colorbar(scatter1, ax=ax1, label='# Cell Types', fraction=0.046, pad=0.04)

# Right: Reconstruction error (RMSE)
rmse = adata.obs['rmse'].values
scatter2 = ax2.scatter(
    coords[:, 0],
    coords[:, 1],
    c=rmse,
    cmap='YlOrRd',
    s=8,
    vmin=np.percentile(rmse, 5),
    vmax=np.percentile(rmse, 95)
)
ax2.set_title(f'Reconstruction Error (RMSE)\n(Mean: {rmse.mean():.3f}, Median: {np.median(rmse):.3f})', 
              fontsize=12, fontweight='bold')
ax2.axis('equal')
ax2.invert_yaxis()
ax2.set_xticks([])
ax2.set_yticks([])
plt.colorbar(scatter2, ax=ax2, label='RMSE', fraction=0.046, pad=0.04)

plt.tight_layout()
output_path_sparsity = os.path.join(OUTPUT_DIR, "spatial_sparsity_and_error.pdf")
plt.savefig(output_path_sparsity, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_path_sparsity}")
plt.close()

# ============================================================================
# FIGURE 4: CELL TYPE CO-OCCURRENCE HEATMAP
# ============================================================================

print("Creating cell type co-occurrence matrix...")

# Build co-occurrence matrix: how often do cell types appear together?
cooccurrence = np.zeros((len(cell_types), len(cell_types)))

for i in range(adata.shape[0]):
    # Get cell types present in this spot (proportion > 1%)
    present = np.array([adata.obs[f'prop_{ct}'].values[i] > 0.01 for ct in cell_types])
    
    # Update co-occurrence matrix
    for j in range(len(cell_types)):
        for k in range(len(cell_types)):
            if present[j] and present[k]:
                cooccurrence[j, k] += 1

# Normalize by prevalence
for j in range(len(cell_types)):
    if cooccurrence[j, j] > 0:
        cooccurrence[j, :] /= cooccurrence[j, j]

# Plot
plt.figure(figsize=(14, 12))
sns.heatmap(
    cooccurrence,
    xticklabels=cell_types,
    yticklabels=cell_types,
    cmap='YlOrRd',
    cbar_kws={'label': 'Co-occurrence Frequency'},
    vmin=0,
    vmax=1
)
plt.title('Cell Type Co-occurrence Matrix\n(How often do cell types appear together?)', 
          fontsize=14, fontweight='bold')
plt.xticks(rotation=45, ha='right', fontsize=8)
plt.yticks(rotation=0, fontsize=8)
plt.tight_layout()

output_path_cooccur = os.path.join(OUTPUT_DIR, "celltype_cooccurrence.pdf")
plt.savefig(output_path_cooccur, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved: {output_path_cooccur}")
plt.close()

# ============================================================================
# SUMMARY
# ============================================================================

print()
print("="*80)
print("VISUALIZATION SUMMARY")
print("="*80)
print()

print("Generated visualizations:")
print(f"  1. {output_path_top9}")
print(f"     → Spatial maps of 9 most prevalent cell types")
print(f"  2. {output_path_dominant}")
print(f"     → Dominant cell type and confidence per spot")
print(f"  3. {output_path_sparsity}")
print(f"     → Sparsity (# cell types) and reconstruction error maps")
print(f"  4. {output_path_cooccur}")
print(f"     → Cell type co-occurrence matrix")
print()

print("Key Observations:")
print(f"  - Most prevalent: {sorted_celltypes[0]} ({prevalence[sorted_celltypes[0]]} spots)")
print(f"  - Least prevalent: {sorted_celltypes[-1]} ({prevalence[sorted_celltypes[-1]]} spots)")
print(f"  - Mean sparsity: {n_celltypes.mean():.2f} cell types/spot")
print(f"  - Mean RMSE: {rmse.mean():.4f}")
print()

if n_celltypes.mean() > 10:
    print("⚠️  NOTE: Sparsity is high (>10 cell types/spot)")
    print("    This suggests hyperparameters may need adjustment for")
    print("    stronger sparsity enforcement (higher τ and λ).")
    print()
elif n_celltypes.mean() < 2:
    print("⚠️  NOTE: Sparsity is very low (<2 cell types/spot)")
    print("    This suggests hyperparameters may be too aggressive.")
    print()
else:
    print("✓ Sparsity looks biologically reasonable (2-10 cell types/spot)")
    print()

print("✓ Visualization complete!")
print()