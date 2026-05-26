# Important Note: this README file is created by Claude AI.

---

# Data-Driven Hyperparameter Optimization for Accurate Sparse Deconvolution of Spatial Transcriptomic Data

**Author:** Habiba Hammouda  
**Course:** Bioinformatics II Final Research Paper

---

## Links

- **GitHub Repository:** https://github.com/HabibaHammouda/wispr
- **Original WISpR Paper:** Erdogan & Eroglu (2025). *Sparse deconvolution of cell type medleys in spatial transcriptomics.* PLOS Computational Biology. https://doi.org/10.1371/journal.pcbi.1013169
- **WISpR Code (Zenodo):** https://doi.org/10.5281/zenodo.11109636

### Datasets

**Human Developing Heart (scRNA-seq + Spatial):**
- https://www.spatialresearch.org/resources-published-datasets/doi-10-1016-j-cell-2019-11-025/

**Mouse Brain Visium Data:**
- Available through 10X Genomics datasets
- L1 Hippocampus loom file from Linnarsson Lab Mouse Brain Atlas

**Note on Data Extraction:** Datasets were extracted using R scripts with specialized libraries (loompy, Seurat). The extraction process required navigating format conversions (.loom → .h5ad → 10X MTX format) and handling complex nested metadata structures. See `extract_hippocampus_loom.py` and `preprocess_wispr_seurat.R` for implementation details.

---

## Problem Definition

### The Challenge: Mixed Cellular Signals in Spatial Transcriptomics

**Spatial transcriptomics** (10X Visium) measures gene expression while preserving spatial tissue coordinates. However, each capture spot (~55 μm diameter) contains 1-10 cells whose signals are **blended together**. The individual cell types are hidden.

**Deconvolution** aims to unmix these signals:
```
y = X · β

where:
  y = spot's measured gene expression (observed)
  X = reference cell-type signatures (from scRNA-seq)
  β = cell-type proportions (unknown — what we want to estimate)
```

### Why This Problem Matters

Current deconvolution methods either:
1. Use **fixed hyperparameters** across all spots (fast but potentially inaccurate)
2. Or perform **spot-specific optimization** (slow but theoretically more accurate)

**No one had systematically tested whether the slow, careful approach is actually necessary.**

---

## The WISpR Algorithm

**WISpR** (Weight-Induced Sparse Regression) by Erdogan & Eroglu (2025) is a sparse deconvolution method designed for spatial transcriptomics.

### Key Features

1. **Biological Sparsity Enforcement**  
   Forces each spot to contain only a few cell types (realistic: a 55 μm spot physically cannot hold many cells)

2. **Spot-Specific Hyperparameter Tuning**  
   Optimizes two critical parameters **separately for every single spot**:
   - **τ (tau):** Threshold for removing weak cell-type predictions (sparsity control)
   - **λ (lambda):** L2 regularization penalty (stability control)

3. **Iterative Thresholding**  
   Inspired by SINDy (Sparse Identification of Nonlinear Dynamics), repeatedly prunes weak predictions until only well-supported cell types remain

### The GridSearch Baseline

The original WISpR paper uses **GridSearchCV** to find optimal (τ, λ) for each spot:
- **τ range:** [0.001, 0.500]
- **λ range:** {0.001, 0.01, 0.1, 1.0, 10.0}
- **5-fold cross-validation, 3 repeats**
- **Runtime:** ~32 minutes for 2,110 spots

---

## Research Question

**Is spot-specific GridSearch optimization actually necessary, or could simple fixed hyperparameters work just as well?**

### Hypothesis

If spot-specific tuning is truly necessary, then **all fixed manual settings should produce:**
- Higher reconstruction error (RMSE)
- Biologically impossible cell counts (too many or too few cell types per spot)

---

## Experimental Design

### Dataset
- **Tissue:** Mouse hippocampus (10X Visium)
- **Spots:** 2,110
- **Cell Types:** 26 (from Linnarsson Lab Mouse Brain Atlas L1 clustering)
- **Reference Cells:** 1,754 cells (filtered to 25-250 cells per type)
- **DEGs Selected:** 2,000 genes (Gini coefficient ≥ 0.30, mean Gini = 0.857)

### Comparison: Baseline vs. Manual Settings

| Method | τ (Threshold) | λ (Penalty) | Description |
|--------|---------------|-------------|-------------|
| **GridSearch Baseline** | Spot-specific | Spot-specific | Optimized per spot (gold standard) |
| **Manual 1: Loose τ** | 0.05 | 10 | Very permissive threshold |
| **Manual 2: Strict τ** | 0.30 | 10 | Very aggressive threshold |
| **Manual 3: Weak λ** | 0.15 | 5 | Low regularization |
| **Manual 4: Strong λ** | 0.15 | 50 | High regularization |

### Evaluation Metrics

1. **RMSE** (Root Mean Squared Error): Reconstruction accuracy  
2. **Sparsity:** Number of cell types per spot (biological realism)  
3. **Marker Gene Validation:** Pearson correlation with canonical markers  
4. **Statistical Tests:** Wilcoxon rank-sum, Cohen's d, Cliff's δ

---

## Project Structure

```
wispr/
├── README.md                              # This file
├── data/                                  # Raw data (not included - download separately)
│   └── l1_hippocampus.loom               # Linnarsson Mouse Brain Atlas
├── preprocessed_data/                     # Processed data files
│   ├── scrna_preprocessed.h5ad           # scRNA reference
│   └── spatial_preprocessed.h5ad         # Spatial transcriptomics
├── wispr_analysis/                        # Reference matrix and DEGs
│   ├── reference_matrix_X.csv            # Cell-type gene signatures
│   └── selected_degs.csv                 # DEG metadata (Gini coefficients)
├── wispr_results/                         # Deconvolution outputs
│   ├── spatial_deconvolved.h5ad          # Final predictions
│   ├── marker_correlations.csv           # Validation results
│   └── *.pdf                             # Visualization plots
├── sensitivity_results/                   # Manual hyperparameter experiments
│   ├── exp1_tau_loose_proportions.csv
│   ├── exp2_tau_strict_proportions.csv
│   ├── exp3_lambda_weak_proportions.csv
│   ├── exp4_lambda_strong_proportions.csv
│   └── summary_statistics.csv
└── code/                                  # Analysis scripts
    ├── extract_hippocampus_loom.py       # Data extraction from loom
    ├── preprocess_spatial_and_scrna.py   # Normalization & QC
    ├── preprocess_wispr_seurat.R         # R preprocessing pipeline
    ├── deg_selection_gini.py             # DEG selection (Gini coefficient)
    ├── wispr_deconvolution.py            # WISpR baseline implementation
    ├── hyperparameter_sensitivity.py     # Manual hyperparameter tests
    ├── statistical_validation.py         # Statistical comparisons
    ├── validate_with_markers.py          # Marker gene validation
    ├── visualize_spatial.py              # Spatial maps & plots
    └── subsample_data.R                  # Data subsampling utilities
```

### Script Ownership

- **Capital letter R scripts** (from `zenodo/` folder): Original WISpR paper code
- **Lowercase/descriptive Python/R scripts** (listed above): My implementations

---

## Analysis Pipeline

### 1. Data Extraction & Preprocessing
**Scripts:** `extract_hippocampus_loom.py`, `preprocess_spatial_and_scrna.py`, `preprocess_wispr_seurat.R`

- Extract mouse brain scRNA-seq from `.loom` file
- Filter cell types: 25 ≤ cells ≤ 250 (following Andersson et al. 2020)
- Normalize to 10,000 counts per cell/spot
- Log-transform: log₁₀(count + 1)
- Select top 2,000 variable genes

### 2. DEG Selection
**Script:** `deg_selection_gini.py`

- Calculate **Gini coefficient** for each gene across cell types
- Gini = measure of expression inequality (0 = uniform, 1 = highly cell-type-specific)
- Filter: Gini ≥ 0.30, detection rate ≥ 5%
- Select top 2,000 DEGs (mean Gini = 0.857)
- Build reference matrix **X** (genes × cell types)

### 3. WISpR Baseline Deconvolution
**Script:** `wispr_deconvolution.py`

- Implement WISpR algorithm:
  1. Initialize with non-negative least squares (NNLS)
  2. Iterative soft thresholding: β = max(β_temp - τ, 0)
  3. Convergence check
- GridSearchCV for spot-specific (τ, λ)
- Save proportions, RMSE, sparsity metrics

### 4. Manual Hyperparameter Testing
**Script:** `hyperparameter_sensitivity.py`

- Run WISpR with **4 fixed (τ, λ) settings**
- Same algorithm, no optimization per spot
- Compare RMSE, sparsity, runtime

### 5. Validation
**Scripts:** `validate_with_markers.py`, `statistical_validation.py`

- **Marker Gene Validation:** Correlate predictions with canonical markers (Mbp, Gfap, Cx3cr1, Cldn5, etc.)
- **Statistical Tests:** Wilcoxon rank-sum, effect sizes (Cohen's d, Cliff's δ)

### 6. Visualization
**Script:** `visualize_spatial.py`

- Spatial maps of cell-type distributions
- Sparsity heatmaps
- RMSE distributions
- Co-occurrence matrices

---

## Key Results

### GridSearch Baseline Performance

| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Mean Sparsity** | 4.22 cell types/spot | Biologically realistic |
| **Median Sparsity** | 4 | Consistent with Visium resolution |
| **Mean RMSE** | 0.4545 | Low reconstruction error |
| **Median RMSE** | 0.2801 | Accurate predictions |
| **Optimal τ** | ~0.50 (82.7% of spots) | Aggressive sparsity needed |
| **Optimal λ** | ~1.0 (median) | Moderate regularization |

### Biological Validation

**Marker Gene Correlations (Pearson r, all p < 10⁻²³⁰):**
- Oligodendrocytes (Mbp, Mog) → **r = 0.95**
- Immune cells (Cx3cr1, Ptprc) → **r = 0.81**
- Vascular (Cldn5, Pecam1) → **r = 0.79**

**Spatial Patterns:**
- Oligodendrocytes: Central white-matter regions (consistent with biology)
- Vascular: Scattered network-like distribution (consistent with biology)
- Immune: Diffuse surveillance pattern (consistent with biology)
- Rare neurons: Specific anatomical zones only (consistent with biology)

---

## Manual Hyperparameter Failure

### Quantitative Comparison

| Method | Mean Sparsity | RMSE Increase | Cohen's d | p-value |
|--------|---------------|---------------|-----------|---------|
| **GridSearch Baseline** | 4.22 | — | — | — |
| Manual 1: Loose τ (0.05, 10) | 12.3 | **+91.3%** | -2.64 | < 10⁻¹⁰⁰ |
| Manual 2: Strict τ (0.30, 10) | 18.7 | **+90.6%** | -2.41 | < 10⁻¹⁰⁰ |
| Manual 3: Weak λ (0.15, 5) | 26.1 | **+95.7%** | -2.76 | < 10⁻¹⁰⁰ |
| Manual 4: Strong λ (0.15, 50) | 14.5 | **+41.6%** | -0.77 | < 10⁻¹⁰⁰ |

**Effect Size Interpretation:**
- All **p < 10⁻¹⁰⁰** → statistically significant
- **Cliff's δ ≈ -1.0** → baseline wins on nearly every single spot
- **Cohen's d > 0.8** → large practical effect

### Biological Implausibility

- **12-26 cell types per spot** is physically impossible (55 μm spot holds ~5-8 cells)
- Manual settings either **overfit noise** (detecting ghost cell types) or **fail to converge** properly
- Spatial patterns become **biologically incoherent**

---

## Key Findings

### 1. Spot-Specific Optimization is Essential — Not Optional

Every fixed hyperparameter configuration failed dramatically. The "shortcut" approach of manually setting τ and λ does not work.

### 2. Tissue Heterogeneity Defeats One-Size-Fits-All

Different tissue regions (white matter, grey matter, vasculature, edges) have:
- Different cell compositions
- Different signal-to-noise ratios
- Different expression variability

**No single (τ, λ) can fit all regions.** GridSearch adapts to local context; fixed settings cannot.

### 3. Aggressive Sparsity is Genuinely Needed

82.7% of spots pushed τ to the upper grid bound (~0.50). This suggests:
- The default grid may need expansion (τ > 0.50)
- Biological sparsity constraints are tighter than initially assumed

### 4. The 32-Minute Cost is Worth It

GridSearch is slow, but it's the difference between:
- **Biologically meaningful results** (baseline)
- **Nonsensical, impossible predictions** (manual shortcuts)

### 5. Marker Gene Validation is Powerful

Independent validation with canonical markers:
- Confirms predictions are real biology (not artifacts)
- **Caught a reference mislabel:** A cluster labeled "Neurons" actually correlated r=0.955 with oligodendrocyte markers (Mbp, Mog)

---

## Dependencies

### Python
```
scanpy >= 1.9
pandas >= 1.5
numpy >= 1.23
scipy >= 1.10
matplotlib >= 3.6
seaborn >= 0.12
scikit-learn >= 1.2
loompy >= 3.0
tqdm
```

### R
```
Seurat >= 4.3
loomR
tidyverse
spatstat
Matrix
```

---

## Usage

### Quick Start

```bash
# 1. Extract data from loom file
python code/extract_hippocampus_loom.py

# 2. Preprocess scRNA and spatial data
python code/preprocess_spatial_and_scrna.py
Rscript code/preprocess_wispr_seurat.R

# 3. Select DEGs using Gini coefficient
python code/deg_selection_gini.py

# 4. Run WISpR baseline (GridSearch)
python code/wispr_deconvolution.py

# 5. Test manual hyperparameters
python code/hyperparameter_sensitivity.py

# 6. Validate with markers
python code/validate_with_markers.py

# 7. Generate visualizations
python code/visualize_spatial.py

# 8. Statistical comparison
python code/statistical_validation.py
```

### Expected Runtime
- Data extraction: ~5 minutes
- Preprocessing: ~10 minutes
- DEG selection: ~2 minutes
- **WISpR GridSearch: ~32 minutes** (2,110 spots)
- Manual tests: ~2 hours total (4 experiments × ~30 min each)
- Validation: ~5 minutes
- Visualization: ~3 minutes

---

## Limitations

1. **Single Tissue:** Only mouse hippocampus tested
2. **Matched Reference:** No cross-species or cross-tissue mismatch scenarios
3. **Grid Boundary Issue:** 82.7% of spots hit τ upper bound → wider range may improve results
4. **No Method Comparison:** Only WISpR tested; did not benchmark against RCTD, Cell2location, SPOTlight, etc.
5. **Computational Cost:** GridSearch is slow; not yet practical for very large datasets (Visium-HD with millions of spots)

---

## Future Directions

### Immediate Extensions
1. **Test on additional tissues:** Heart, tumor microenvironment, developing embryo
2. **Expand hyperparameter grid:** τ ∈ [0.001, 1.0], λ ∈ [0.001, 100]
3. **Benchmark against competitors:** RCTD, Cell2location, Stereoscope, SPOTlight

### Methodological Improvements
1. **Faster optimization:** Bayesian optimization (GP-based), surrogate models, or transfer learning from similar spots
2. **Visium-HD scaling:** Adapt to subcellular resolution (2 μm bins, millions of spots)
3. **Reference mismatch robustness:** Cross-species, cross-tissue, cross-platform scenarios
4. **Adaptive grid design:** Automatically adjust search ranges based on tissue characteristics

---

## References

1. **Erdogan & Eroglu (2025).** *Sparse deconvolution of cell type medleys in spatial transcriptomics.* PLOS Computational Biology. https://doi.org/10.1371/journal.pcbi.1013169

2. **Andersson et al. (2020).** *Single-cell and spatial transcriptomics enables probabilistic inference of cell type topography.* Communications Biology. https://doi.org/10.1038/s42003-020-01247-y

3. **10X Genomics Visium.** https://www.10xgenomics.com/products/visium

4. **Linnarsson Lab Mouse Brain Atlas.** http://mousebrain.org/

---

## Acknowledgments

- **Erdogan & Eroglu** for developing the WISpR algorithm and making code publicly available
- **Linnarsson Lab** for the Mouse Brain Atlas dataset
- **10X Genomics** for Visium spatial transcriptomics technology
- **Course instructor and TA** for guidance throughout the project

---

**Last Updated:** May 2026
