"""
WISpR Deconvolution Algorithm
==============================

Implements the Weight-Induced Sparse Regression (WISpR) algorithm from
Erdogan & Eroglu (2025) for spatial transcriptomics deconvolution.

For each spot:
  1. GridSearch to find optimal τ (threshold) and λ (penalty) parameters
  2. Solve sparse regression: y ≈ X·β with iterative thresholding
  3. Return cell type proportions β (sparse, non-negative)

"""

import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import Lasso, Ridge
from sklearn.metrics import mean_squared_error
from scipy.optimize import nnls  # Non-negative least squares
from scipy import sparse
import warnings
warnings.filterwarnings('ignore')
import os
from tqdm import tqdm  # Progress bar

print("="*80)
print("WISpR DECONVOLUTION ALGORITHM")
print("="*80)
print()

# ============================================================================
# CONFIGURATION
# ============================================================================

# Input files
SPATIAL_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/preprocessed_data/spatial_preprocessed.h5ad"
REFERENCE_MATRIX_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_analysis/reference_matrix_X.csv"
DEG_FILE = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_analysis/selected_degs.csv"

OUTPUT_DIR = "C:/Users/user/Desktop/My_First_BioinfoII_Research_Paper/BioProject/wispr_results/"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# # WISpR hyperparameters - OLD RANGES (too low, leading to dense solutions)
# TAU_RANGE = np.linspace(0.001, 0.1, 10)  # Threshold values to try
# LAMBDA_RANGE = np.logspace(-3, 1, 10)    # Penalty values to try (log scale)
# MAX_ITER = 100  # Maximum iterations for convergence
# TOL = 1e-6      # Convergence tolerance

# WISpR hyperparameters - CORRECTED RANGES
TAU_RANGE = np.linspace(0.01, 0.5, 10)    # Higher thresholds
LAMBDA_RANGE = np.logspace(0, 3, 10)      # Higher penalties (1 to 1000)
MAX_ITER = 100
TOL = 1e-6

# For speed, you can test on subset first
TEST_MODE = False  # Set True to test on 100 spots only
N_TEST_SPOTS = 100

print(f"Configuration:")
print(f"  - Spatial data: {SPATIAL_FILE}")
print(f"  - Reference matrix: {REFERENCE_MATRIX_FILE}")
print(f"  - Output: {OUTPUT_DIR}")
print(f"  - τ (threshold) range: [{TAU_RANGE.min():.4f}, {TAU_RANGE.max():.4f}] ({len(TAU_RANGE)} values)")
print(f"  - λ (penalty) range: [{LAMBDA_RANGE.min():.4f}, {LAMBDA_RANGE.max():.4f}] ({len(LAMBDA_RANGE)} values)")
print(f"  - Test mode: {TEST_MODE}")
if TEST_MODE:
    print(f"  - Testing on first {N_TEST_SPOTS} spots only")
print()

# ============================================================================
# STEP 1: Load Data
# ============================================================================

print("="*80)
print("STEP 1: LOADING DATA")
print("="*80)
print()

print("Loading spatial data...")
adata_spatial = sc.read_h5ad(SPATIAL_FILE)
print(f"  ✓ Spatial: {adata_spatial.shape[0]} spots × {adata_spatial.shape[1]} genes")

print("Loading reference matrix X...")
reference_df = pd.read_csv(REFERENCE_MATRIX_FILE, index_col=0)
X = reference_df.values  # genes × cell_types
cell_types = reference_df.columns.tolist()
deg_names = reference_df.index.tolist()

print(f"  ✓ Reference matrix X: {X.shape[0]} genes × {X.shape[1]} cell types")
print(f"  ✓ Cell types: {len(cell_types)}")
print()

# Subset spatial data to DEGs only
print("Subsetting spatial data to DEGs...")
adata_spatial = adata_spatial[:, deg_names].copy()
print(f"  ✓ Spatial data subset: {adata_spatial.shape[0]} spots × {adata_spatial.shape[1]} genes")
print()

# Test mode: use only first N spots
if TEST_MODE:
    print(f"⚠️  TEST MODE: Using only first {N_TEST_SPOTS} spots")
    adata_spatial = adata_spatial[:N_TEST_SPOTS].copy()
    print(f"  Reduced to {adata_spatial.shape[0]} spots")
    print()

# ============================================================================
# STEP 2: WISpR Algorithm Functions
# ============================================================================

print("="*80)
print("STEP 2: DEFINING WISpR FUNCTIONS")
print("="*80)
print()

def wispr_iterative_thresholding(y, X, tau, lambda_reg, max_iter=MAX_ITER, tol=TOL):
    """
    WISpR iterative thresholding algorithm (Algorithm 1 from paper).
    
    Solves: min ||y - X·β||² + λ||β||₁
    Subject to: β ≥ 0 (non-negativity)
    
    Uses soft thresholding with threshold τ to enforce sparsity.
    
    Parameters:
    -----------
    y : array (n_genes,)
        Gene expression vector for one spot
    X : array (n_genes, n_celltypes)
        Reference matrix
    tau : float
        Threshold parameter (controls sparsity)
    lambda_reg : float
        Regularization parameter (L1 penalty)
    max_iter : int
        Maximum iterations
    tol : float
        Convergence tolerance
    
    Returns:
    --------
    beta : array (n_celltypes,)
        Cell type proportions (sparse, non-negative)
    """
    n_genes, n_celltypes = X.shape
    
    # Initialize beta with non-negative least squares
    beta, _ = nnls(X, y)
    
    # Iterative thresholding
    for iteration in range(max_iter):
        beta_old = beta.copy()
        
        # Gradient step (least squares update)
        residual = y - X @ beta
        gradient = -2 * X.T @ residual + 2 * lambda_reg * beta
        
        # Update with step size (using line search approximation)
        step_size = 0.01
        beta = beta - step_size * gradient
        
        # Soft thresholding
        beta = np.sign(beta) * np.maximum(np.abs(beta) - tau, 0)
        
        # Non-negativity constraint
        beta = np.maximum(beta, 0)
        
        # Check convergence
        if np.linalg.norm(beta - beta_old) < tol:
            break
    
    return beta


def gridsearch_hyperparameters(y, X, tau_range, lambda_range):
    """
    GridSearch to find optimal τ and λ for a single spot.
    
    Tries all combinations of τ and λ, selects the one with minimum
    reconstruction error (RMSE).
    
    Parameters:
    -----------
    y : array (n_genes,)
        Gene expression vector for one spot
    X : array (n_genes, n_celltypes)
        Reference matrix
    tau_range : array
        Threshold values to try
    lambda_range : array
        Penalty values to try
    
    Returns:
    --------
    best_beta : array (n_celltypes,)
        Cell type proportions with optimal parameters
    best_tau : float
        Optimal threshold
    best_lambda : float
        Optimal penalty
    best_error : float
        Reconstruction error (RMSE)
    """
    best_error = np.inf
    best_beta = None
    best_tau = None
    best_lambda = None
    
    for tau in tau_range:
        for lambda_reg in lambda_range:
            # Run WISpR with these parameters
            beta = wispr_iterative_thresholding(y, X, tau, lambda_reg)
            
            # Calculate reconstruction error
            y_pred = X @ beta
            error = np.sqrt(mean_squared_error(y, y_pred))
            
            # Update best parameters if error improved
            if error < best_error:
                best_error = error
                best_beta = beta
                best_tau = tau
                best_lambda = lambda_reg
    
    return best_beta, best_tau, best_lambda, best_error


print("WISpR algorithm functions defined:")
print("  1. wispr_iterative_thresholding() - Core sparse regression")
print("  2. gridsearch_hyperparameters() - Hyperparameter optimization")
print()

# ============================================================================
# STEP 3: Run WISpR on All Spots
# ============================================================================

print("="*80)
print("STEP 3: RUNNING WISpR DECONVOLUTION")
print("="*80)
print()

n_spots = adata_spatial.shape[0]
n_celltypes = X.shape[1]

print(f"Deconvolving {n_spots} spots...")
print(f"  - Grid search: {len(TAU_RANGE)} τ × {len(LAMBDA_RANGE)} λ = {len(TAU_RANGE)*len(LAMBDA_RANGE)} combinations per spot")
print(f"  - Total evaluations: {n_spots * len(TAU_RANGE) * len(LAMBDA_RANGE):,}")
print()
print("This may take 10-30 minutes depending on your system...")
print("Progress bar will show estimated time remaining.")
print()

# Initialize result arrays
cell_type_proportions = np.zeros((n_spots, n_celltypes))
optimal_tau = np.zeros(n_spots)
optimal_lambda = np.zeros(n_spots)
reconstruction_errors = np.zeros(n_spots)

# Get expression matrix as dense array
if sparse.issparse(adata_spatial.X):
    Y = adata_spatial.X.toarray()
else:
    Y = adata_spatial.X

# Run WISpR on each spot with progress bar
for i in tqdm(range(n_spots), desc="Deconvolving spots"):
    y = Y[i, :]  # Gene expression for spot i
    
    # GridSearch for optimal parameters
    beta, tau_opt, lambda_opt, error = gridsearch_hyperparameters(
        y, X, TAU_RANGE, LAMBDA_RANGE
    )
    
    # Store results
    cell_type_proportions[i, :] = beta
    optimal_tau[i] = tau_opt
    optimal_lambda[i] = lambda_opt
    reconstruction_errors[i] = error

print()
print("  ✓ Deconvolution complete!")
print()

# ============================================================================
# STEP 4: Post-Processing and Analysis
# ============================================================================

print("="*80)
print("STEP 4: POST-PROCESSING RESULTS")
print("="*80)
print()

print("Analyzing deconvolution results...")

# Calculate sparsity (number of non-zero cell types per spot)
sparsity = (cell_type_proportions > 1e-6).sum(axis=1)

print(f"Sparsity analysis:")
print(f"  - Mean cell types per spot: {sparsity.mean():.2f}")
print(f"  - Median cell types per spot: {np.median(sparsity):.0f}")
print(f"  - Range: [{sparsity.min():.0f}, {sparsity.max():.0f}]")
print()

print(f"Hyperparameter statistics:")
print(f"  - τ (threshold): mean={optimal_tau.mean():.4f}, median={np.median(optimal_tau):.4f}")
print(f"  - λ (penalty): mean={optimal_lambda.mean():.4f}, median={np.median(optimal_lambda):.4f}")
print()

print(f"Reconstruction error:")
print(f"  - Mean RMSE: {reconstruction_errors.mean():.4f}")
print(f"  - Median RMSE: {np.median(reconstruction_errors):.4f}")
print()

# Calculate cell type prevalence (how many spots have each cell type)
cell_type_prevalence = (cell_type_proportions > 1e-6).sum(axis=0)

print(f"Cell type prevalence (detected in how many spots):")
for i, ct in enumerate(cell_types):
    prevalence = cell_type_prevalence[i]
    pct = 100 * prevalence / n_spots
    mean_prop = cell_type_proportions[:, i].mean()
    print(f"  - {ct:30s}: {prevalence:4d} spots ({pct:5.1f}%) | Mean prop: {mean_prop:.4f}")
print()

# ============================================================================
# STEP 5: Save Results
# ============================================================================

print("="*80)
print("STEP 5: SAVING RESULTS")
print("="*80)
print()

print("Saving deconvolution results...")

# Add results to spatial AnnData object
for i, ct in enumerate(cell_types):
    adata_spatial.obs[f'prop_{ct}'] = cell_type_proportions[:, i]

adata_spatial.obs['n_celltypes'] = sparsity
adata_spatial.obs['optimal_tau'] = optimal_tau
adata_spatial.obs['optimal_lambda'] = optimal_lambda
adata_spatial.obs['rmse'] = reconstruction_errors

# Save annotated spatial data
spatial_results_path = os.path.join(OUTPUT_DIR, "spatial_deconvolved.h5ad")
adata_spatial.write_h5ad(spatial_results_path)
print(f"  ✓ Saved annotated spatial data: {spatial_results_path}")

# Save cell type proportions as CSV
proportions_df = pd.DataFrame(
    cell_type_proportions,
    index=adata_spatial.obs_names,
    columns=cell_types
)
proportions_path = os.path.join(OUTPUT_DIR, "cell_type_proportions.csv")
proportions_df.to_csv(proportions_path)
print(f"  ✓ Saved proportions: {proportions_path}")

# Save hyperparameters
hyperparams_df = pd.DataFrame({
    'spot': adata_spatial.obs_names,
    'tau': optimal_tau,
    'lambda': optimal_lambda,
    'rmse': reconstruction_errors,
    'n_celltypes': sparsity
})
hyperparams_path = os.path.join(OUTPUT_DIR, "hyperparameters.csv")
hyperparams_df.to_csv(hyperparams_path, index=False)
print(f"  ✓ Saved hyperparameters: {hyperparams_path}")

print()

# ============================================================================
# STEP 6: Generate Visualizations
# ============================================================================

print("="*80)
print("STEP 6: GENERATING VISUALIZATIONS")
print("="*80)
print()

print("Creating deconvolution visualizations...")

# Figure 1: Overall statistics
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
fig.suptitle('WISpR Deconvolution Results', fontsize=16, fontweight='bold')

# Plot 1: Sparsity distribution
axes[0, 0].hist(sparsity, bins=range(0, int(sparsity.max())+2), edgecolor='black', alpha=0.7)
axes[0, 0].axvline(sparsity.mean(), color='red', linestyle='--', linewidth=2, label=f'Mean: {sparsity.mean():.1f}')
axes[0, 0].set_xlabel('Number of Cell Types per Spot')
axes[0, 0].set_ylabel('Number of Spots')
axes[0, 0].set_title('Sparsity: Cell Types per Spot')
axes[0, 0].legend()

# Plot 2: Reconstruction error distribution
axes[0, 1].hist(reconstruction_errors, bins=50, edgecolor='black', alpha=0.7)
axes[0, 1].axvline(reconstruction_errors.mean(), color='red', linestyle='--', linewidth=2, 
                   label=f'Mean: {reconstruction_errors.mean():.3f}')
axes[0, 1].set_xlabel('RMSE')
axes[0, 1].set_ylabel('Number of Spots')
axes[0, 1].set_title('Reconstruction Error Distribution')
axes[0, 1].legend()

# Plot 3: Optimal τ distribution
axes[1, 0].hist(optimal_tau, bins=50, edgecolor='black', alpha=0.7, color='green')
axes[1, 0].set_xlabel('Optimal τ (Threshold)')
axes[1, 0].set_ylabel('Number of Spots')
axes[1, 0].set_title('Spot-Specific Threshold Distribution')

# Plot 4: Optimal λ distribution
axes[1, 1].hist(np.log10(optimal_lambda), bins=50, edgecolor='black', alpha=0.7, color='purple')
axes[1, 1].set_xlabel('Optimal log₁₀(λ) (Penalty)')
axes[1, 1].set_ylabel('Number of Spots')
axes[1, 1].set_title('Spot-Specific Penalty Distribution')

plt.tight_layout()
stats_plot_path = os.path.join(OUTPUT_DIR, "deconvolution_statistics.pdf")
plt.savefig(stats_plot_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved statistics plots: {stats_plot_path}")
plt.close()

# Figure 2: Cell type prevalence
fig, ax = plt.subplots(1, 1, figsize=(12, 6))
sorted_idx = np.argsort(cell_type_prevalence)[::-1]
sorted_celltypes = [cell_types[i] for i in sorted_idx]
sorted_prevalence = cell_type_prevalence[sorted_idx]

ax.barh(range(len(cell_types)), sorted_prevalence, color='steelblue', edgecolor='black')
ax.set_yticks(range(len(cell_types)))
ax.set_yticklabels(sorted_celltypes, fontsize=9)
ax.set_xlabel('Number of Spots Detected')
ax.set_title('Cell Type Prevalence Across Tissue', fontsize=14, fontweight='bold')
ax.invert_yaxis()

plt.tight_layout()
prevalence_plot_path = os.path.join(OUTPUT_DIR, "cell_type_prevalence.pdf")
plt.savefig(prevalence_plot_path, dpi=300, bbox_inches='tight')
print(f"  ✓ Saved prevalence plot: {prevalence_plot_path}")
plt.close()

print()

# ============================================================================
# SUMMARY
# ============================================================================

print("="*80)
print("WISpR DECONVOLUTION SUMMARY")
print("="*80)
print()

print(f"Input:")
print(f"  - Spots deconvolved: {n_spots:,}")
print(f"  - Genes (DEGs): {X.shape[0]:,}")
print(f"  - Cell types: {n_celltypes}")
print()

print(f"Algorithm Parameters:")
print(f"  - τ range: [{TAU_RANGE.min():.4f}, {TAU_RANGE.max():.4f}]")
print(f"  - λ range: [{LAMBDA_RANGE.min():.4f}, {LAMBDA_RANGE.max():.4f}]")
print(f"  - Grid size: {len(TAU_RANGE)} × {len(LAMBDA_RANGE)} = {len(TAU_RANGE)*len(LAMBDA_RANGE)} combinations")
print()

print(f"Results:")
print(f"  - Mean sparsity: {sparsity.mean():.2f} cell types/spot")
print(f"  - Mean RMSE: {reconstruction_errors.mean():.4f}")
print(f"  - Most prevalent cell type: {cell_types[np.argmax(cell_type_prevalence)]}")
print(f"  - Least prevalent cell type: {cell_types[np.argmin(cell_type_prevalence)]}")
print()

print(f"Output Files:")
print(f"  1. {spatial_results_path}")
print(f"     → Annotated spatial data with cell type proportions")
print(f"  2. {proportions_path}")
print(f"     → Cell type proportion matrix (spots × cell_types)")
print(f"  3. {hyperparams_path}")
print(f"     → Optimal hyperparameters per spot")
print(f"  4. {stats_plot_path}")
print(f"     → Deconvolution statistics plots")
print(f"  5. {prevalence_plot_path}")
print(f"     → Cell type prevalence across tissue")
print()
print("✓ WISpR deconvolution complete!")