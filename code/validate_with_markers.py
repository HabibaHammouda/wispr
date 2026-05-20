#!/usr/bin/env python3
"""
Marker Gene Validation for WISpR Results
=========================================

Validates deconvolution by checking if predicted cell types correlate
with expression of known marker genes.

For example:
- Neurons → High Slc17a7, Snap25
- Oligodendrocytes → High Mbp, Mog
- Immune cells → High Cx3cr1, Ptprc
- Vascular → High Cldn5, Pecam1
"""

import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr
import os

print("="*80)
print("MARKER GENE VALIDATION")
print("="*80)
print()

# ============================================================================
# CONFIGURATION
# ============================================================================

SPATIAL_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/preprocessed_data/spatial_preprocessed.h5ad"
RESULTS_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_results/spatial_deconvolved.h5ad"
OUTPUT_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_results/"

# Known marker genes for mouse brain cell types
MARKER_GENES = {
    'Neurons': ['Slc17a7', 'Snap25', 'Rbfox3', 'Tubb3', 'Map2'],
    'Oligodendrocytes': ['Mbp', 'Mog', 'Plp1', 'Olig1', 'Olig2'],
    'Astrocytes': ['Gfap', 'Aqp4', 'Aldh1l1', 'Slc1a2', 'Slc1a3'],
    'Immune': ['Cx3cr1', 'Ptprc', 'Aif1', 'C1qa', 'C1qb'],
    'Vascular': ['Cldn5', 'Pecam1', 'Flt1', 'Vwf', 'Tie1'],
    'Ependymal': ['Foxj1', 'Tmem212', 'Dynlrb2']
}

print(f"Configuration:")
print(f"  - Spatial data (full gene set): {SPATIAL_FILE}")
print(f"  - WISpR results: {RESULTS_FILE}")
print(f"  - Output: {OUTPUT_DIR}")
print()

# ============================================================================
# LOAD DATA
# ============================================================================

print("="*80)
print("STEP 1: LOADING DATA")
print("="*80)
print()

print("Loading spatial data with full gene set...")
adata_full = sc.read_h5ad(SPATIAL_FILE)
print(f"  ✓ Loaded: {adata_full.shape[0]} spots × {adata_full.shape[1]} genes")

print("Loading WISpR results...")
adata_results = sc.read_h5ad(RESULTS_FILE)
print(f"  ✓ Loaded: {adata_results.shape[0]} spots")
print()

# Get cell type proportions
cell_type_cols = [col for col in adata_results.obs.columns if col.startswith('prop_')]
cell_types = [col.replace('prop_', '') for col in cell_type_cols]
print(f"  ✓ Cell types: {len(cell_types)}")
print()

# ============================================================================
# CHECK MARKER GENE AVAILABILITY
# ============================================================================

print("="*80)
print("STEP 2: CHECKING MARKER GENE AVAILABILITY")
print("="*80)
print()

available_markers = {}
missing_markers = {}

for category, genes in MARKER_GENES.items():
    available = [g for g in genes if g in adata_full.var_names]
    missing = [g for g in genes if g not in adata_full.var_names]
    
    available_markers[category] = available
    missing_markers[category] = missing
    
    print(f"{category}:")
    print(f"  Available: {len(available)}/{len(genes)} - {available}")
    if missing:
        print(f"  Missing: {missing}")
    print()

# ============================================================================
# CALCULATE MARKER GENE EXPRESSION SCORES
# ============================================================================

print("="*80)
print("STEP 3: CALCULATING MARKER EXPRESSION SCORES")
print("="*80)
print()

print("Computing average marker expression for each category...")

marker_scores = {}

for category, genes in available_markers.items():
    if len(genes) > 0:
        # Get expression of marker genes
        marker_expr = adata_full[:, genes].X
        
        # Average across markers (log-normalized values)
        if hasattr(marker_expr, 'toarray'):
            marker_expr = marker_expr.toarray()
        
        score = marker_expr.mean(axis=1)
        marker_scores[category] = score
        
        print(f"  {category}: {len(genes)} markers, score range [{score.min():.3f}, {score.max():.3f}]")
    else:
        print(f"  {category}: No markers available!")

print()

# ============================================================================
# CORRELATE MARKERS WITH CELL TYPE PROPORTIONS
# ============================================================================

print("="*80)
print("STEP 4: CORRELATING MARKERS WITH PREDICTED PROPORTIONS")
print("="*80)
print()

print("Computing correlations between marker scores and WISpR predictions...")
print()

# For each broad category, find which predicted cell types correlate
correlation_results = []

for category, score in marker_scores.items():
    print(f"{category} markers:")
    
    # Correlate with all cell type proportions
    for ct in cell_types:
        props = adata_results.obs[f'prop_{ct}'].values
        
        # Pearson and Spearman correlations
        pearson_r, pearson_p = pearsonr(score, props)
        spearman_r, spearman_p = spearmanr(score, props)
        
        correlation_results.append({
            'marker_category': category,
            'cell_type': ct,
            'pearson_r': pearson_r,
            'pearson_p': pearson_p,
            'spearman_r': spearman_r,
            'spearman_p': spearman_p
        })
        
        # Print significant correlations
        if abs(pearson_r) > 0.3 and pearson_p < 0.01:
            print(f"  → {ct:30s}: r={pearson_r:+.3f} (p={pearson_p:.2e})")
    
    print()

# Save correlation results
corr_df = pd.DataFrame(correlation_results)
corr_path = os.path.join(OUTPUT_DIR, "marker_correlations.csv")
corr_df.to_csv(corr_path, index=False)
print(f"✓ Saved correlations: {corr_path}")
print()

# ============================================================================
# VISUALIZE MARKER VS PREDICTION CORRELATIONS
# ============================================================================

print("="*80)
print("STEP 5: GENERATING VALIDATION PLOTS")
print("="*80)
print()

print("Creating marker gene validation visualizations...")

# Create correlation heatmap
pivot_pearson = corr_df.pivot(index='cell_type', columns='marker_category', values='pearson_r')

plt.figure(figsize=(10, 12))
sns.heatmap(
    pivot_pearson,
    cmap='RdBu_r',
    center=0,
    vmin=-0.5,
    vmax=0.5,
    cbar_kws={'label': 'Pearson r'},
    linewidths=0.5,
    linecolor='gray'
)
plt.title('Marker Gene Validation: Correlation Heatmap\n(Marker Expression vs Predicted Proportions)', 
          fontsize=12, fontweight='bold')
plt.xlabel('Marker Gene Category', fontsize=10)
plt.ylabel('Predicted Cell Type', fontsize=10)
plt.xticks(rotation=45, ha='right')
plt.yticks(rotation=0, fontsize=8)
plt.tight_layout()

heatmap_path = os.path.join(OUTPUT_DIR, "marker_correlation_heatmap.pdf")
plt.savefig(heatmap_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved heatmap: {heatmap_path}")
plt.close()

# Create scatter plots for top correlations
print("Creating scatter plots for top correlations...")

# Get top 6 correlations
top_corrs = corr_df.nlargest(6, 'pearson_r')

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, (idx, row) in enumerate(top_corrs.iterrows()):
    if i >= 6:
        break
    
    ax = axes[i]
    category = row['marker_category']
    ct = row['cell_type']
    r = row['pearson_r']
    p = row['pearson_p']
    
    # Get data
    x = marker_scores[category]
    y = adata_results.obs[f'prop_{ct}'].values
    
    # Scatter plot
    ax.scatter(x, y, alpha=0.3, s=5, color='steelblue')
    ax.set_xlabel(f'{category} Marker Score', fontsize=9)
    ax.set_ylabel(f'{ct} Proportion', fontsize=9)
    ax.set_title(f'{category} → {ct}\nr={r:.3f}, p={p:.2e}', fontsize=10)
    
    # Add regression line
    z = np.polyfit(x, y, 1)
    p_line = np.poly1d(z)
    x_sorted = np.sort(x)
    ax.plot(x_sorted, p_line(x_sorted), "r--", linewidth=2, alpha=0.8)

plt.suptitle('Top 6 Marker-Prediction Correlations', fontsize=14, fontweight='bold')
plt.tight_layout()

scatter_path = os.path.join(OUTPUT_DIR, "top_marker_correlations.pdf")
plt.savefig(scatter_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved scatter plots: {scatter_path}")
plt.close()

# ============================================================================
# SUMMARY
# ============================================================================

print()
print("="*80)
print("VALIDATION SUMMARY")
print("="*80)
print()

print("Marker Gene Categories Tested:")
for category, genes in available_markers.items():
    print(f"  - {category}: {len(genes)} markers")
print()

print("Strongest Marker-Prediction Correlations:")
top_5 = corr_df.nlargest(5, 'pearson_r')
for i, (idx, row) in enumerate(top_5.iterrows(), 1):
    print(f"  {i}. {row['marker_category']:20s} → {row['cell_type']:30s}: r={row['pearson_r']:+.3f}")
print()

print("Output Files:")
print(f"  1. {corr_path}")
print(f"     → Full correlation table")
print(f"  2. {heatmap_path}")
print(f"     → Correlation heatmap")
print(f"  3. {scatter_path}")
print(f"     → Top 6 scatter plots")
print()

# Check if validation makes sense
high_corr_count = (corr_df['pearson_r'].abs() > 0.4).sum()
sig_corr_count = ((corr_df['pearson_r'].abs() > 0.3) & (corr_df['pearson_p'] < 0.01)).sum()

print(f"Quality Check:")
print(f"  - Strong correlations (|r| > 0.4): {high_corr_count}")
print(f"  - Significant correlations (|r| > 0.3, p < 0.01): {sig_corr_count}")
print()

if sig_corr_count > 10:
    print("✓ Good validation! Multiple cell types show expected marker correlations.")
else:
    print("⚠️  Few significant correlations detected.")
    print("    This could indicate:")
    print("    - Hyperparameters need further tuning")
    print("    - Cell type labels are too granular (cluster IDs vs broad types)")
    print("    - Spatial heterogeneity within cell type clusters")

print()
print("✓ Marker validation complete!")
print()