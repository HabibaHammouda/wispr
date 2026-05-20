#!/usr/bin/env python3
"""
HYPERPARAMETER SENSITIVITY - USING YOUR WORKING WISPR CODE
No pysindy - uses the exact same method that already worked for you!
"""

import scanpy as sc
import pandas as pd
import numpy as np
from scipy.optimize import nnls
from tqdm import tqdm
import os
import warnings
warnings.filterwarnings('ignore')

print("="*80)
print("WISPR HYPERPARAMETER SENSITIVITY - SIMPLE & WORKING")
print("="*80)
print()

# ============================================================================
# PATHS - USING YOUR EXACT FORMAT
# ============================================================================

BASE = "C:\\Users\\user\\Desktop\\My_First_BioinfoII_Research_Paper\\BioProject"
SPATIAL_FILE = f"{BASE}\\preprocessed_data\\spatial_preprocessed.h5ad"
REFERENCE_FILE = f"{BASE}\\wispr_analysis\\reference_matrix_X.csv"
OUTPUT_DIR = f"{BASE}\\sensitivity_results"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ============================================================================
# WISPR ALGORITHM (FROM YOUR WORKING CODE)
# ============================================================================

def wispr_iterative_thresholding(y, X, tau, lambda_reg, max_iter=100, tol=1e-6):
    """
    WISpR algorithm - EXACT copy from your working wispr_deconvolution.py
    """
    n_genes, n_celltypes = X.shape
    
    # Initialize with non-negative least squares
    beta, _ = nnls(X, y)
    
    for iteration in range(max_iter):
        beta_old = beta.copy()
        
        # Soft thresholding
        residual = y - X @ beta
        gradient = -2 * X.T @ residual
        beta_temp = beta - (1.0 / (2 * lambda_reg)) * gradient
        
        # Apply threshold
        beta = np.maximum(beta_temp - tau, 0)
        
        # Check convergence
        if np.linalg.norm(beta - beta_old) < tol:
            break
    
    return beta

def deconvolve_all_spots(adata, X, tau, lambda_reg):
    """
    Run WISpR on all spots with FIXED parameters
    """
    n_spots, n_genes = adata.shape
    n_celltypes = X.shape[1]
    
    proportions = np.zeros((n_spots, n_celltypes))
    rmse_values = []
    
    # Get expression matrix
    Y = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X
    
    for i in tqdm(range(n_spots), desc="Deconvolving"):
        y = Y[i, :]
        
        if np.sum(y) == 0:
            # Empty spot
            rmse_values.append(np.nan)
            continue
        
        # Run WISpR
        beta = wispr_iterative_thresholding(y, X, tau, lambda_reg)
        proportions[i, :] = beta
        
        # Calculate RMSE
        y_pred = X @ beta
        rmse = np.sqrt(np.mean((y - y_pred) ** 2))
        rmse_values.append(rmse)
    
    return proportions, np.array(rmse_values)

# ============================================================================
# LOAD DATA
# ============================================================================

print("Loading data...")
adata = sc.read_h5ad(SPATIAL_FILE)
ref_df = pd.read_csv(REFERENCE_FILE, index_col=0)

# Subset to shared genes
adata = adata[:, ref_df.index].copy()
X = ref_df.values
cell_types = ref_df.columns.tolist()

print(f"✓ Spatial: {adata.shape[0]} spots × {adata.shape[1]} genes")
print(f"✓ Reference: {X.shape[0]} genes × {X.shape[1]} cell types")
print()

# ============================================================================
# EXPERIMENTS - 4 SIMPLE TESTS
# ============================================================================

experiments = {
    'exp1_tau_loose': {'tau': 0.05, 'lambda': 10, 'desc': 'LOOSE threshold'},
    'exp2_tau_strict': {'tau': 0.30, 'lambda': 10, 'desc': 'STRICT threshold'},
    'exp3_lambda_weak': {'tau': 0.15, 'lambda': 5, 'desc': 'WEAK penalty'},
    'exp4_lambda_strong': {'tau': 0.15, 'lambda': 50, 'desc': 'STRONG penalty'}
}

print("EXPERIMENTS:")
for name, cfg in experiments.items():
    print(f"  {name}: τ={cfg['tau']}, λ={cfg['lambda']} - {cfg['desc']}")
print()
print("Estimated time: ~2 hours total")
print()

input("Press ENTER to start...")

# ============================================================================
# RUN ALL EXPERIMENTS
# ============================================================================

results_summary = []

for exp_name, config in experiments.items():
    print(f"\n{'='*80}")
    print(f"RUNNING: {exp_name}")
    print(f"  {config['desc']}")
    print(f"  τ={config['tau']}, λ={config['lambda']}")
    print(f"{'='*80}\n")
    
    # Deconvolve
    props, rmse_vals = deconvolve_all_spots(
        adata, X,
        tau=config['tau'],
        lambda_reg=config['lambda']
    )
    
    # Calculate statistics
    n_celltypes_per_spot = (props > 0).sum(axis=1)
    mean_sparsity = n_celltypes_per_spot.mean()
    median_sparsity = np.median(n_celltypes_per_spot)
    mean_rmse = np.nanmean(rmse_vals)
    median_rmse = np.nanmedian(rmse_vals)
    
    print(f"\n✓ COMPLETE!")
    print(f"  Mean sparsity: {mean_sparsity:.2f} cell types/spot")
    print(f"  Median sparsity: {median_sparsity:.0f}")
    print(f"  Mean RMSE: {mean_rmse:.4f}")
    print(f"  Median RMSE: {median_rmse:.4f}")
    
    # Save
    props_df = pd.DataFrame(props, columns=cell_types)
    props_df.to_csv(f"{OUTPUT_DIR}\\{exp_name}_proportions.csv", index=False)
    
    pd.DataFrame({
        'spot_id': range(len(rmse_vals)),
        'tau': config['tau'],
        'lambda': config['lambda'],
        'n_celltypes': n_celltypes_per_spot,
        'rmse': rmse_vals
    }).to_csv(f"{OUTPUT_DIR}\\{exp_name}_metrics.csv", index=False)
    
    print(f"  Saved: {exp_name}_proportions.csv")
    print(f"  Saved: {exp_name}_metrics.csv")
    
    # Store summary
    results_summary.append({
        'experiment': exp_name,
        'description': config['desc'],
        'tau': config['tau'],
        'lambda': config['lambda'],
        'mean_sparsity': mean_sparsity,
        'median_sparsity': median_sparsity,
        'mean_rmse': mean_rmse,
        'median_rmse': median_rmse,
        'pct_over5_celltypes': (n_celltypes_per_spot > 5).mean() * 100
    })

# ============================================================================
# FINAL SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print("SUMMARY OF ALL EXPERIMENTS")
print(f"{'='*80}\n")

summary_df = pd.DataFrame(results_summary)
print(summary_df.to_string(index=False))

summary_df.to_csv(f"{OUTPUT_DIR}\\summary_statistics.csv", index=False)
print(f"\n✓ Saved: summary_statistics.csv")

print(f"\n{'='*80}")
print("YOUR BASELINE (for comparison):")
print(f"{'='*80}")
print("  Mean sparsity: 4.22 cell types/spot")
print("  Median sparsity: 4")
print("  Mean RMSE: 0.4545")
print("  Median RMSE: 0.2801")

print(f"\n{'='*80}")
print("DONE! All results in:")
print(f"  {OUTPUT_DIR}")
print(f"{'='*80}\n")