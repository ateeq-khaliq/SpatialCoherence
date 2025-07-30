# 🧬 SpatialCoherence Tutorial: Complete User Guide

[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue.svg)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/Version-1.0.0-orange.svg)](https://github.com/ateeq-khaliq/SpatialCoherence)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

**Comprehensive guide for analyzing spatial organization patterns of multicellular ecotypes in spatial transcriptomics data**

---

## 📋 Table of Contents

1. [What is SpatialCoherence?](#what-is-spatialcoherence)
2. [Installation & Setup](#installation--setup)
3. [Understanding Spatial Ecotypes](#understanding-spatial-ecotypes)
4. [Complete Parameter Guide](#complete-parameter-guide)
5. [Step-by-Step Tutorial](#step-by-step-tutorial)
6. [Real-World Example: PDAC Analysis](#real-world-example-pdac-analysis)
7. [Parameter Selection Guide](#parameter-selection-guide)
8. [Results Interpretation](#results-interpretation)
9. [Troubleshooting](#troubleshooting)
10. [Publication Guidelines](#publication-guidelines)

---

## 🔬 What is SpatialCoherence?

SpatialCoherence quantifies how **spatial ecotypes** (multicellular tissue environments) organize within spatial transcriptomics data. Unlike cell-type analysis, it focuses on **tissue-level organization patterns**.

### 🎯 **Core Concept: Discovery, Not Validation**

**❌ What SpatialCoherence Does NOT Do:**
- Assume which ecotypes "should" be organized
- Use predetermined biological knowledge for classification
- Validate pre-existing hypotheses about organization

** What SpatialCoherence DOES Do:**
- **Discovers** which spatial ecotypes are organized vs disorganized
- **Measures** spatial coherence objectively using neighborhood statistics
- **Classifies** based on user-defined, data-driven criteria
- **Enables** biological interpretation after discovery

### 🧪 **Scientific Workflow:**
1. **Tool measures** spatial coherence scores (0-1 scale)
2. **User defines** classification threshold based on their data
3. **Tool classifies** ecotypes as Organized/Disorganized
4. **User interprets** biological meaning of discoveries

---

## 📦 Installation & Setup

### Prerequisites

```r
# Check R version (≥ 4.1.0 required)
R.version.string

# Install required packages
install.packages(c("devtools", "ggplot2", "reshape2", "dplyr"))

# Optional packages for enhanced features
if (!require(BiocManager)) install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
install.packages(c("circlize", "yaml"))
```

### Install SpatialCoherence

```r
# Install from GitHub
devtools::install_github("ateeq-khaliq/SpatialCoherence")

# Load and verify
library(SpatialCoherence)
get_package_version()
#> [1] "1.0.0"
```

### Verify Installation

```r
# Test basic functionality
get_package_info()
#> $package
#> [1] "SpatialCoherence"
#> 
#> $version
#> [1] "1.0.0"
#> 
#> $author
#> [1] "Ateeq Khaliq"

# Check available functions
ls("package:SpatialCoherence")
```

---

## 📊 Computational & Statistical Methods

This section details the mathematical foundations, statistical tests, and computational algorithms underlying SpatialCoherence analysis.

---

### 🧮 **Core Spatial Coherence Algorithm**

#### **1. Neighborhood Definition**

For each spot *i*, the spatial neighborhood *N(i)* is defined as:

```
N(i) = {j : d(i,j) ≤ τ + ε, j ≠ i}
```

Where:
- `d(i,j)` = Euclidean distance between spots i and j
- `τ` = typical distance between adjacent spots (auto-calculated)
- `ε` = tolerance factor (30% of typical distance)

**Distance Calculation:**
```
d(i,j) = √[(x_i - x_j)² + (y_i - y_j)²]
```

For Visium data, typical distance is estimated as:
```
τ = median({d(i,j) : i,j ∈ random_sample, d(i,j) > 0})
```

#### **2. Observed Spatial Score**

For ecotype *k*, the observed spatial score is:

```
S_obs(k) = (1/|V_k|) × Σ_{i ∈ V_k} (1/|N(i)|) × Σ_{j ∈ N(i)} I(c_j = k)
```

Where:
- `V_k` = set of spots belonging to ecotype k
- `|V_k|` = number of spots in ecotype k  
- `N(i)` = neighborhood of spot i
- `I(c_j = k)` = indicator function (1 if spot j is ecotype k, 0 otherwise)
- `c_j` = ecotype label of spot j

**Biological interpretation:** Average proportion of same-ecotype neighbors per spot.

#### **3. Theoretical Maximum Score**

The theoretical maximum (complete spatial clustering) is calculated using hexagonal packing geometry:

```
S_max(k) = max(0, (6|V_k| - 12√a - 6) / |V_k|)
```

Where:
```
a = √(4|V_k| / (6√3))
```

This represents the maximum possible spatial coherence assuming optimal hexagonal arrangement.

#### **4. Expected Random Score**

The expected score under random spatial arrangement is estimated through permutation testing:

```
S_random(k) = (1/R) × Σ_{r=1}^{R} S_obs^{(r)}(k)
```

Where:
- `R` = number of random permutations (default: 100)
- `S_obs^{(r)}(k)` = observed score for ecotype k in permutation r

#### **5. Final Coherence Score**

The normalized spatial coherence score is:

```
Coherence(k) = (S_obs(k) - S_random(k)) / (S_max(k) - S_random(k))
```

**Properties:**
- Range: [0, 1] where 0 = random, 1 = perfectly organized
- Normalized by theoretical maximum for the given number of spots
- Accounts for random spatial arrangement baseline

---

### 📈 **Statistical Tests and Corrections**

#### **1. Permutation Testing**

**Null Hypothesis:** Spatial arrangement of ecotype k is random

**Test Statistic:** 
```
T = S_obs(k)
```

**P-value Calculation:**
```
p-value = (1 + Σ_{r=1}^{R} I(S_obs^{(r)}(k) ≥ S_obs(k))) / (R + 1)
```

Where `S_obs^{(r)}(k)` is the observed score in the r-th random permutation.

**Statistical Power:** With R=100 permutations, minimum detectable p-value = 0.01

#### **2. Treatment Effect Analysis**

**Student's t-test** for comparing coherence between treatment groups:

**Test Statistic:**
```
t = (x̄₁ - x̄₂) / (s_pooled × √(1/n₁ + 1/n₂))
```

Where:
- `x̄₁, x̄₂` = mean coherence scores in groups 1 and 2
- `s_pooled` = pooled standard deviation
- `n₁, n₂` = sample sizes

**Pooled Standard Deviation:**
```
s_pooled = √[((n₁-1)s₁² + (n₂-1)s₂²) / (n₁+n₂-2)]
```

**Degrees of Freedom:**
```
df = n₁ + n₂ - 2
```

#### **3. Effect Size Calculation (Cohen's d)**

**Formula:**
```
d = (x̄₁ - x̄₂) / s_pooled
```

**Interpretation Guidelines:**
- `|d| < 0.2`: Negligible effect
- `0.2 ≤ |d| < 0.5`: Small effect  
- `0.5 ≤ |d| < 0.8`: Medium effect
- `|d| ≥ 0.8`: Large effect

**Confidence Interval for Cohen's d:**
```
CI = d ± t_(α/2,df) × SE(d)
```

Where:
```
SE(d) = √[(n₁+n₂)/(n₁×n₂) + d²/(2(n₁+n₂))]
```

#### **4. Multiple Testing Correction**

When testing multiple ecotypes, p-values are adjusted using:

**Benjamini-Hochberg (FDR) Correction:**
```
p_adj(i) = min(1, p(i) × m / i)
```

Where:
- `p(i)` = i-th smallest p-value
- `m` = total number of tests
- Applied sequentially from largest to smallest p-value

**Bonferroni Correction (Conservative):**
```
p_adj = min(1, p × m)
```

---

### 🔢 **Spatial Statistics & Autocorrelation**

#### **1. Moran's I (Global Spatial Autocorrelation)**

**Formula:**
```
I = (n/W) × [Σᵢ Σⱼ wᵢⱼ(xᵢ - x̄)(xⱼ - x̄)] / [Σᵢ(xᵢ - x̄)²]
```

Where:
- `n` = number of spatial units
- `W = Σᵢ Σⱼ wᵢⱼ` = sum of all spatial weights
- `wᵢⱼ` = spatial weight between units i and j
- `xᵢ` = value at location i
- `x̄` = mean value

**Interpretation:**
- `I > 0`: Positive spatial autocorrelation (clustering)
- `I = 0`: No spatial autocorrelation (random)
- `I < 0`: Negative spatial autocorrelation (dispersion)

**Expected Value under Null:**
```
E[I] = -1/(n-1)
```

**Variance:**
```
Var[I] = [n²-3n+3]S₁ - nS₂ + 3W² / [(n-1)(n-2)(n-3)W²]
```

#### **2. Local Indicators of Spatial Association (LISA)**

**Local Moran's I:**
```
Iᵢ = (xᵢ - x̄) × Σⱼ wᵢⱼ(xⱼ - x̄) / σ²
```

Where `σ²` is the variance of x.

**Significance Testing:**
Under null hypothesis of no local spatial autocorrelation:
```
Z(Iᵢ) = (Iᵢ - E[Iᵢ]) / √Var[Iᵢ]
```

#### **3. Spatial Weight Matrices**

**Binary Contiguity (Used in SpatialCoherence):**
```
wᵢⱼ = {1 if spots i and j are neighbors
       {0 otherwise
```

**Row-Standardized Weights:**
```
w'ᵢⱼ = wᵢⱼ / Σⱼ wᵢⱼ
```

**Distance-Based Weights (Alternative):**
```
wᵢⱼ = {1/d(i,j) if d(i,j) ≤ threshold
       {0        otherwise
```

---

### 🧪 **Quality Control Metrics**

#### **1. Neighborhood Completeness**

**Metric:** Proportion of spots with complete neighborhoods
```
Completeness = |{i : |N(i)| = k}| / |V|
```

Where:
- `k` = expected number of neighbors (typically 6 for Visium)
- `V` = all spots

#### **2. Spatial Coverage**

**Metric:** Effective spatial resolution
```
Coverage = Area(ConvexHull(spots)) / (|V| × typical_area_per_spot)
```

#### **3. Ecotype Representation**

**Minimum Sample Size Check:**
```
Adequacy(k) = {1 if |V_k| ≥ threshold
              {0 otherwise
```

**Spatial Distribution Check:**
```
Spread(k) = |{samples : V_k ∩ sample ≠ ∅}| / |samples|
```

---

### 🔬 **Advanced Statistical Methods**

#### **1. Compositional Analysis**

**Relative Abundance Difference:**
```
RelAbund(k) = P(k|Organized) - P(k|Disorganized)
```

Where:
```
P(k|Organized) = |V_k ∩ V_organized| / |V_organized|
P(k|Disorganized) = |V_k ∩ V_disorganized| / |V_disorganized|
```

**Statistical Test:** Chi-square test for independence
```
χ² = Σₖ [(Oₖ - Eₖ)² / Eₖ]
```

Where:
- `Oₖ` = observed count of ecotype k in organized regions
- `Eₖ` = expected count under independence

#### **2. Threshold Sensitivity Analysis**

**Stability Metric:**
```
Stability(t₁, t₂) = |Classifications(t₁) ∩ Classifications(t₂)| / |V|
```

**Optimal Threshold Selection:**
Maximize separation between organized and disorganized groups:
```
t* = argmax_t [μ_organized(t) - μ_disorganized(t)] / σ_pooled(t)
```

#### **3. Sample-Level Analysis**

**Sample Coherence Score:**
```
Sample_Coherence(s) = Σₖ wₖ × Coherence(k,s)
```

Where `wₖ` is the proportion of ecotype k in sample s.

**Between-Sample Variability:**
```
ICC = (MSB - MSW) / (MSB + (k-1)MSW)
```

Where:
- `MSB` = mean square between samples
- `MSW` = mean square within samples  
- `k` = average number of spots per sample

---

### 💻 **Computational Complexity**

#### **Time Complexity**

**Basic Algorithm:**
- Neighborhood construction: O(n²) worst case, O(n log n) with spatial indexing
- Coherence calculation: O(n × k) where k = average neighborhood size
- Permutation testing: O(R × n × k) where R = number of permutations

**Total:** O(R × n²) worst case, O(R × n log n) optimized

**Space Complexity:**
- Spatial coordinates: O(n)
- Neighborhood matrices: O(n × k)
- Permutation storage: O(R × n) if cached

**Total:** O(R × n) typical case

#### **Optimization Strategies**

**1. Spatial Indexing:**
Use k-d trees or spatial hashing to reduce neighbor search from O(n) to O(log n)

**2. Early Termination:**
For permutation testing, stop when sufficient precision achieved:
```
Stop when: SE(p-value) < tolerance
```

**3. Parallel Processing:**
Permutations are embarrassingly parallel:
```
Speedup ≈ min(R, number_of_cores)
```

---

### 📐 **Geometric Considerations**

#### **1. Boundary Effects**

**Edge Correction Factor:**
```
Correction(i) = |N_theoretical(i)| / |N_observed(i)|
```

Applied to coherence calculation:
```
S_corrected(k) = Σᵢ Correction(i) × LocalScore(i,k)
```

#### **2. Anisotropy Detection**

**Directional Analysis:**
Calculate coherence in different directions to detect anisotropic patterns.

**Ripley's K-function (Isotropic Test):**
```
K(r) = λ⁻¹ × E[number of points within distance r]
```

Where λ is the intensity of the point process.

#### **3. Scale Dependencies**

**Multi-scale Analysis:**
Test coherence at different neighborhood sizes:
```
Coherence(k,r) = f(neighborhood_radius = r)
```

**Scale Selection:**
Choose r that maximizes signal-to-noise ratio:
```
r* = argmax_r [Signal(r) / Noise(r)]
```

---

### 🔍 **Model Diagnostics**

#### **1. Residual Analysis**

**Spatial Residuals:**
```
e(i) = Observed(i) - Expected(i)
```

**Moran's I on Residuals:**
Test for remaining spatial autocorrelation in residuals.

#### **2. Goodness of Fit**

**R-squared for Spatial Models:**
```
R² = 1 - SS_residual / SS_total
```

**Pseudo R-squared (for non-linear models):**
```
Pseudo-R² = (LogLik_full - LogLik_null) / LogLik_null
```

#### **3. Cross-Validation**

**Spatial Cross-Validation:**
- Split data spatially (not randomly) to avoid spatial leakage
- Use buffer zones between training and testing regions

**Leave-One-Out:**
```
CV_score = (1/n) × Σᵢ L(yᵢ, ŷ₋ᵢ)
```

Where `ŷ₋ᵢ` is prediction for point i using all other points.

---

### 📚 **Statistical Assumptions**

#### **1. Spatial Independence (Null Model)**

**Assumption:** Under null hypothesis, ecotype labels are randomly distributed across spatial locations.

**Validation:** Permutation testing preserves marginal distributions while breaking spatial structure.

#### **2. Stationarity**

**Assumption:** Spatial process is stationary (statistical properties don't vary across space).

**Testing:** Use moving window analysis to detect non-stationarity.

#### **3. Isotropy**

**Assumption:** Spatial relationships are the same in all directions.

**Testing:** Compare coherence in different directional sectors.

#### **4. Normal Approximation**

**When Valid:** For large sample sizes (n > 30), coherence scores approximately normal.

**Central Limit Theorem Application:**
```
(Coherence - μ) / (σ/√n) ~ N(0,1)
```

---

### 🎯 **Implementation Notes**

#### **Numerical Stability**

**1. Division by Zero Protection:**
```
Coherence = {0                           if S_max = S_random
            {(S_obs - S_random) / ε      if |S_max - S_random| < ε  
            {(S_obs - S_random) / (S_max - S_random)  otherwise
```

**2. Overflow Prevention:**
Use log-space calculations for very large or very small probabilities.

**3. Precision Control:**
All calculations use double precision (64-bit) floating point arithmetic.

#### **Algorithmic Choices**

**1. Permutation Strategy:**
- **Complete permutation:** Shuffle all ecotype labels
- **Preserves:** Marginal ecotype frequencies
- **Breaks:** Spatial structure only

**2. Tie Handling:**
When multiple neighbors at same distance, include all (may exceed k neighbors).

**3. Missing Data:**
Spots with insufficient neighborhood information are excluded from analysis.

---

This comprehensive statistical foundation ensures that SpatialCoherence provides robust, well-validated measures of spatial organization with proper statistical inference and effect size quantification.

---

## 🔬 Understanding Spatial Ecotypes

### What are Spatial Ecotypes?

**Spatial ecotypes** are multicellular tissue environments characterized by:

- **Coordinated gene expression** across multiple cell types
- **Spatially coherent organization** within tissue sections
- **Functional specialization** (e.g., hypoxic regions, immune niches)
- **Consistent composition** across samples

### Examples in Cancer Research:

| Spatial Ecotype | Description | Expected Organization |
|------------------|-------------|----------------------|
| **SE01** | Hypoxic Environment | High (coordinated stress response) |
| **SE02** | Immune Infiltration | Low (dynamic, dispersed) |
| **SE03** | Stromal Barrier | High (structural organization) |
| **SE04** | Proliferative Zone | High (coordinated growth) |
| **SE05** | Vascular Rich | Medium (organized but linear) |
| **SE06** | Metabolic Active | Medium (clustered activity) |
| **SE07** | Ductal Organized | High (epithelial structure) |
| **SE08** | Acinar Remnant | Medium (degraded organization) |
| **SE09** | Inflammatory Zone | Low (dynamic immune response) |
| **SE10** | Mixed Transitional | Low (heterogeneous mixture) |

**Key Point:** These are **hypotheses** - let the tool discover the actual organization patterns!

---

## ⚙️ Complete Parameter Guide

### 🎯 **Classification Parameters**

#### **`coherence_threshold`** (Default: 0.47)
**What it controls:** The cutoff for classifying ecotypes as "Organized" vs "Disorganized"

**How to choose:**
```r
# Method 1: Use default (well-validated across studies)
coherence_threshold = 0.47

# Method 2: Data-driven selection
distribution <- analyze_coherence_distribution(your_data, "ecotype_col", "sample_col")
print(distribution$suggested_thresholds)
# Choose based on your data distribution

# Method 3: Biological reasoning
coherence_threshold = 0.55  # Strict (fewer organized)
coherence_threshold = 0.40  # Lenient (more organized)
```

**Guidelines by tissue type:**
- **Normal tissue:** 0.55-0.65 (well-organized)
- **Primary tumors:** 0.45-0.55 (mixed patterns)
- **Metastatic sites:** 0.35-0.45 (more chaotic)
- **Treated tissue:** 0.40-0.50 (therapy-altered)

#### **`min_spots_per_ecotype`** (Default: 50)
**What it controls:** Minimum spots required for reliable analysis

**Quality guidelines:**
- **50+ spots:** Basic analysis (minimum)
- **100+ spots:** Good reliability
- **250+ spots:** Excellent statistical power
- **500+ spots:** Robust for all analyses

```r
# Check your data
ecotype_counts <- table(your_seurat$ecotype_column)
print(ecotype_counts)

# Adjust based on data quality
min_spots_per_ecotype = 30   # Exploratory (low quality ok)
min_spots_per_ecotype = 100  # Publication (good quality)
min_spots_per_ecotype = 200  # High precision (strict quality)
```

#### **`min_samples_per_ecotype`** (Default: 3)
**What it controls:** Minimum samples needed for statistical testing

**Statistical guidelines:**
- **3 samples:** Minimum for basic statistics
- **5 samples:** Good for t-tests and comparisons
- **10+ samples:** Robust statistical testing

### 🔬 **Spatial Analysis Parameters**

#### **`n_neighbors`** (Default: 6)
**What it controls:** Number of spatial neighbors considered for each spot

**Technical considerations:**
```r
# Visium platform (hexagonal grid)
n_neighbors = 6  # Standard (captures immediate neighbors)

# Other platforms or custom analysis
n_neighbors = 4  # Conservative (closest neighbors only)
n_neighbors = 8  # Expanded (includes diagonal neighbors)
n_neighbors = 12 # Extended (larger neighborhood)
```

**How to choose:**
- **Dense platforms:** Higher values (8-12)
- **Sparse platforms:** Lower values (4-6)
- **Noisy data:** Lower values (4-6)
- **Clean data:** Higher values (6-8)

### 📊 **Statistical Parameters**

#### **`n_permutations`** (Default: 100)
**What it controls:** Number of random permutations for statistical significance

**Precision vs Speed:**
```r
# Quick exploration
n_permutations = 25   # Fast, less precise

# Standard analysis
n_permutations = 100  # Good balance

# Publication quality
n_permutations = 200  # High precision

# Maximum precision
n_permutations = 500  # Slow but very robust
```

#### **`random_seed`** (Default: NULL)
**What it controls:** Reproducibility of random processes

```r
# For reproducible results
random_seed = 123

# For exploration (different each time)
random_seed = NULL
```

### 🔬 **Analysis Mode Parameters**

#### **`enhanced_analysis`** (Default: FALSE)
**What it controls:** Whether to use advanced spatial algorithms

```r
# Basic analysis (faster, good for exploration)
enhanced_analysis = FALSE

# Advanced analysis (slower, publication quality)
enhanced_analysis = TRUE
```

**When to use enhanced:**
-  Publication-quality results needed
-  Treatment effect analysis
-  Compositional analysis required
-  High-precision spatial patterns

**When to use basic:**
-  Initial data exploration
-  Parameter optimization
-  Large datasets (speed concerns)
-  Proof-of-concept analysis

---

## 📚 Step-by-Step Tutorial

### Step 1: Data Preparation

```r
# Load required libraries
library(SpatialCoherence)
library(Seurat)

# Ensure your Seurat object has required components
# 1. Spatial ecotype assignments
# 2. Sample identifiers  
# 3. Spatial coordinates (in @images)
# 4. Optional: treatment information

# Validate your data
validation <- validate_spatial_data(
  seurat_object = your_seurat,
  ecotype_column = "your_ecotype_column",
  sample_column = "your_sample_column"
)

if (!validation$valid) {
  stop("Data validation failed: ", validation$messages)
} else {
  cat(" Data validation passed\n")
  if (length(validation$messages) > 0) {
    cat("Warnings:\n")
    for(msg in validation$messages) cat(" ", msg, "\n")
  }
}
```

### Step 2: Understand Your Data Distribution

```r
# Analyze your data to guide parameter selection
distribution_analysis <- analyze_coherence_distribution(
  seurat_object = your_seurat,
  ecotype_column = "your_ecotype_column", 
  sample_column = "your_sample_column"
)

# View statistical summary
cat("=== YOUR DATA DISTRIBUTION ===\n")
print(distribution_analysis$summary)

# View suggested thresholds
cat("\n=== SUGGESTED THRESHOLDS ===\n")
print(distribution_analysis$suggested_thresholds)

# Visualize distribution
hist(distribution_analysis$scores, 
     main = "Coherence Score Distribution",
     xlab = "Spatial Coherence Score", 
     ylab = "Frequency",
     breaks = 20)
abline(v = 0.47, col = "red", lwd = 2, lty = 2)  # Default threshold
```

### Step 3: Threshold Sensitivity Analysis

```r
# Test multiple thresholds to understand sensitivity
threshold_analysis <- perform_threshold_analysis(
  seurat_object = your_seurat,
  ecotype_column = "your_ecotype_column",
  sample_column = "your_sample_column",
  thresholds = c(0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60)
)

# View results
cat("=== THRESHOLD SENSITIVITY ANALYSIS ===\n")
print(threshold_analysis$threshold_analysis)

# Choose your threshold based on:
# 1. Biological interpretation
# 2. Desired classification balance
# 3. Data distribution characteristics
chosen_threshold <- 0.47  # Or your preferred value
```

### Step 4: Run Basic Discovery Analysis

```r
# Start with basic analysis for exploration
results_basic <- run_spatial_analysis(
  seurat_object = your_seurat,
  ecotype_column = "your_ecotype_column",
  sample_column = "your_sample_column",
  output_dir = "Basic_Discovery_Analysis",
  
  # Your chosen parameters
  coherence_threshold = chosen_threshold,
  min_spots_per_ecotype = 50,
  random_seed = 123,  # For reproducibility
  
  # Basic settings
  enhanced_analysis = FALSE,
  save_csvs = TRUE,
  create_plots = TRUE,
  verbose = TRUE
)

# Examine discoveries
cat("=== SPATIAL ECOTYPE DISCOVERIES ===\n")
discovery_summary <- data.frame(
  Ecotype = names(results_basic$mean_coherence),
  Coherence_Score = round(results_basic$mean_coherence, 3),
  Classification = results_basic$organization_results[names(results_basic$mean_coherence)]
)
discovery_summary <- discovery_summary[order(-discovery_summary$Coherence_Score), ]
print(discovery_summary)
```

### Step 5: Enhanced Analysis (If Basic Looks Good)

```r
# Run enhanced analysis for publication-quality results
results_enhanced <- run_spatial_analysis(
  seurat_object = your_seurat,
  ecotype_column = "your_ecotype_column",
  sample_column = "your_sample_column",
  treatment_column = "your_treatment_column",  # If available
  output_dir = "Enhanced_Publication_Analysis",
  
  # Optimized parameters
  coherence_threshold = chosen_threshold,
  min_spots_per_ecotype = 100,  # Higher quality
  n_neighbors = 6,
  n_permutations = 200,  # Higher precision
  random_seed = 123,
  
  # Enhanced features
  enhanced_analysis = TRUE,
  analyze_treatments = TRUE,  # If treatment column provided
  save_csvs = TRUE,
  create_plots = TRUE,
  verbose = TRUE
)
```

### Step 6: Treatment Analysis (Optional)

```r
# If you have treatment data, analyze effects
if (!is.null(results_enhanced$treatment_results)) {
  cat("=== TREATMENT EFFECTS DISCOVERED ===\n")
  
  treatment_effects <- results_enhanced$treatment_results
  significant_effects <- treatment_effects[treatment_effects$p_value < 0.05, ]
  
  if (nrow(significant_effects) > 0) {
    cat("Significant treatment effects found:\n")
    for (i in 1:nrow(significant_effects)) {
      ecotype <- significant_effects$ecotype[i]
      effect_size <- significant_effects$effect_size[i]
      p_val <- significant_effects$p_value[i]
      
      direction <- ifelse(effect_size > 0, "INCREASED", "DECREASED")
      magnitude <- ifelse(abs(effect_size) > 0.5, "LARGE", 
                         ifelse(abs(effect_size) > 0.2, "MEDIUM", "SMALL"))
      
      cat(sprintf("- %s: %s organization (%s effect, p = %.4f)\n", 
                  ecotype, direction, magnitude, p_val))
    }
  } else {
    cat("No significant treatment effects detected.\n")
  }
}
```

---

## 🩺 Real-World Example: PDAC Analysis

### Complete PDAC Analysis Workflow

```r
# ============================================================================
# COMPREHENSIVE PDAC SPATIAL ECOTYPE ANALYSIS
# ============================================================================

library(SpatialCoherence)

# Load your PDAC data
# pdac <- readRDS("path/to/your/pdac_data.rds")
# Ensure it has: CompositionCluster_SE, library_id, nac_treatment columns

# Step 1: Validate PDAC data
cat("=== VALIDATING PDAC DATA ===\n")
validation <- validate_spatial_data(
  pdac, 
  ecotype_column = "CompositionCluster_SE",
  sample_column = "library_id"
)

if (!validation$valid) stop(validation$messages)
cat(" PDAC data validated\n")
cat("Samples:", length(unique(pdac$library_id)), "\n")
cat("Spots:", ncol(pdac), "\n") 
cat("Spatial ecotypes:", length(unique(pdac$CompositionCluster_SE)), "\n\n")

# Step 2: Understand PDAC data distribution
cat("=== ANALYZING PDAC DATA DISTRIBUTION ===\n")
pdac_distribution <- analyze_coherence_distribution(
  pdac, "CompositionCluster_SE", "library_id"
)

print(pdac_distribution$summary)
print(pdac_distribution$suggested_thresholds)

# Step 3: PDAC threshold sensitivity
cat("\n=== PDAC THRESHOLD SENSITIVITY ===\n")
pdac_thresholds <- perform_threshold_analysis(
  pdac, "CompositionCluster_SE", "library_id",
  thresholds = c(0.35, 0.40, 0.47, 0.55, 0.60)
)
print(pdac_thresholds$threshold_analysis)

# Step 4: Choose threshold for PDAC
# For this example, using data-driven median approach
pdac_threshold <- pdac_distribution$suggested_thresholds$Threshold[
  pdac_distribution$suggested_thresholds$Approach == "Balanced (Mean)"
]
cat("\nChosen threshold for PDAC:", pdac_threshold, "\n\n")

# Step 5: PDAC spatial ecotype annotations
pdac_annotations <- c(
  "SE01" = "Hypoxic_Environment",
  "SE02" = "Immune_Infiltration",
  "SE03" = "Stromal_Barrier", 
  "SE04" = "Proliferative_Zone",
  "SE05" = "Vascular_Rich",
  "SE06" = "Metabolic_Active",
  "SE07" = "Ductal_Organized",
  "SE08" = "Acinar_Remnant",
  "SE09" = "Inflammatory_Zone",
  "SE10" = "Mixed_Transitional"
)

# Step 6: Comprehensive PDAC analysis
cat("=== RUNNING COMPREHENSIVE PDAC ANALYSIS ===\n")
pdac_results <- run_spatial_analysis(
  seurat_object = pdac,
  ecotype_column = "CompositionCluster_SE",
  sample_column = "library_id",
  treatment_column = "nac_treatment",
  output_dir = "PDAC_SpatialCoherence_Analysis",
  
  # PDAC-optimized parameters
  coherence_threshold = pdac_threshold,
  min_spots_per_ecotype = 100,
  n_neighbors = 6,
  n_permutations = 200,
  random_seed = 42,  # For PDAC reproducibility
  
  # Enhanced PDAC analysis
  enhanced_analysis = TRUE,
  analyze_treatments = TRUE,
  save_csvs = TRUE,
  create_plots = TRUE,
  
  # PDAC annotations
  ecotype_annotations = pdac_annotations,
  
  verbose = TRUE
)

# Step 7: Interpret PDAC discoveries
cat("\n=== PDAC SPATIAL ECOTYPE DISCOVERIES ===\n")
pdac_summary <- data.frame(
  Spatial_Ecotype = names(pdac_results$mean_coherence),
  Functional_Annotation = pdac_annotations[names(pdac_results$mean_coherence)],
  Coherence_Score = round(pdac_results$mean_coherence, 3),
  Classification = pdac_results$organization_results[names(pdac_results$mean_coherence)]
)
pdac_summary <- pdac_summary[order(-pdac_summary$Coherence_Score), ]
print(pdac_summary)

# Step 8: PDAC treatment effects
if (!is.null(pdac_results$treatment_results)) {
  cat("\n=== NAC TREATMENT EFFECTS IN PDAC ===\n")
  nac_effects <- pdac_results$treatment_results
  nac_significant <- nac_effects[nac_effects$p_value < 0.05, ]
  
  if (nrow(nac_significant) > 0) {
    cat("NAC significantly affected these spatial ecotypes:\n")
    for (i in 1:nrow(nac_significant)) {
      ecotype <- nac_significant$ecotype[i]
      annotation <- pdac_annotations[ecotype]
      effect_size <- nac_significant$effect_size[i]
      p_val <- nac_significant$p_value[i]
      
      direction <- ifelse(effect_size > 0, "increased", "decreased")
      cat(sprintf("- %s (%s): %s organization (effect = %.3f, p = %.4f)\n",
                  ecotype, annotation, direction, abs(effect_size), p_val))
    }
  }
}

# Step 9: PDAC biological interpretation
cat("\n=== PDAC BIOLOGICAL INSIGHTS ===\n")
organized_ecotypes <- names(pdac_results$organization_results)[
  pdac_results$organization_results == "Organized"
]
disorganized_ecotypes <- names(pdac_results$organization_results)[
  pdac_results$organization_results == "Disorganized"
]

cat("Organized PDAC environments:\n")
for (ecotype in organized_ecotypes) {
  cat(sprintf("- %s (%s): %.3f coherence\n", 
              ecotype, pdac_annotations[ecotype], 
              pdac_results$mean_coherence[ecotype]))
}

cat("\nDisorganized PDAC environments:\n")
for (ecotype in disorganized_ecotypes) {
  cat(sprintf("- %s (%s): %.3f coherence\n", 
              ecotype, pdac_annotations[ecotype], 
              pdac_results$mean_coherence[ecotype]))
}

cat("\n=== PDAC ANALYSIS COMPLETE ===\n")
cat("Results saved to: PDAC_SpatialCoherence_Analysis/\n")
cat("Key findings ready for publication!\n")
```

---

## 🎯 Parameter Selection Guide

### 🔬 **Study Type-Specific Parameters**

#### **Exploratory Studies**
```r
run_spatial_analysis(
  coherence_threshold = 0.40,      # More permissive
  min_spots_per_ecotype = 30,      # Lower requirements
  n_permutations = 50,            # Faster analysis
  enhanced_analysis = FALSE,       # Basic discovery
  random_seed = NULL              # Allow variation
)
```

#### **Publication Studies**
```r
run_spatial_analysis(
  coherence_threshold = 0.47,      # Well-validated default
  min_spots_per_ecotype = 100,     # Quality control
  n_permutations = 200,           # High precision
  enhanced_analysis = TRUE,        # Advanced algorithms
  random_seed = 123,              # Reproducible
  analyze_treatments = TRUE        # If applicable
)
```

#### **Clinical Validation**
```r
run_spatial_analysis(
  coherence_threshold = 0.55,      # Stringent criteria
  min_spots_per_ecotype = 200,     # High quality
  n_permutations = 500,           # Maximum precision
  enhanced_analysis = TRUE,
  random_seed = 42,               # Fixed seed
  min_samples_per_ecotype = 10     # Statistical power
)
```

### 📊 **Data Type-Specific Parameters**

#### **High-Quality Visium Data**
```r
run_spatial_analysis(
  coherence_threshold = 0.47,
  min_spots_per_ecotype = 100,
  n_neighbors = 6,               # Standard hexagonal
  n_permutations = 200,
  enhanced_analysis = TRUE
)
```

#### **Lower-Quality or Sparse Data**
```r
run_spatial_analysis(
  coherence_threshold = 0.40,      # More lenient
  min_spots_per_ecotype = 50,      # Lower requirements
  n_neighbors = 4,                # Conservative neighborhood
  n_permutations = 100,           # Standard precision
  enhanced_analysis = FALSE        # Less demanding
)
```

#### **Large Datasets (>100k spots)**
```r
run_spatial_analysis(
  coherence_threshold = 0.47,
  min_spots_per_ecotype = 200,     # Higher quality bar
  n_neighbors = 6,
  n_permutations = 100,           # Balance speed/precision
  enhanced_analysis = TRUE,
  parallel = TRUE                 # If implemented
)
```

### 🎯 **Threshold Selection Decision Tree**

```
Do you have prior knowledge about expected organization?
├─ YES: Use biological reasoning
│   ├─ Expect high organization → threshold = 0.55-0.60
│   ├─ Expect mixed patterns → threshold = 0.45-0.50  
│   └─ Expect low organization → threshold = 0.35-0.40
└─ NO: Use data-driven approach
    ├─ Run analyze_coherence_distribution()
    ├─ Examine suggested thresholds
    ├─ Consider distribution shape
    └─ Choose based on desired balance

Want to compare with literature?
├─ Use established threshold (0.47) for comparability
└─ Report sensitivity analysis with multiple thresholds

Publishing results?
├─ Include threshold sensitivity analysis
├─ Justify threshold choice in methods
└─ Report classification statistics
```

---

## 📊 Results Interpretation

### 🎯 **Understanding Coherence Scores**

#### **Score Ranges and Biological Meaning**

| Score Range | Classification | Biological Interpretation | Examples |
|-------------|---------------|---------------------------|----------|
| **0.8-1.0** | Highly Organized | Clear spatial structure, coordinated architecture | Epithelial ducts, organized stroma |
| **0.6-0.8** | Organized | Moderate spatial coherence, structured regions | Proliferative zones, hypoxic cores |
| **0.47-0.6** | Moderately Organized | Some spatial structure, partially organized | Metabolic clusters, vascular regions |
| **0.3-0.47** | Disorganized | Mixed patterns, limited spatial structure | Immune infiltration, transitional zones |
| **0.0-0.3** | Highly Disorganized | Random distribution, no spatial organization | Inflammatory responses, invasion fronts |

#### **Key Metrics to Report**

```r
# Calculate summary statistics for your publication
results_summary <- function(results) {
  coherence_scores <- results$mean_coherence
  organization <- results$organization_results
  
  # Basic statistics
  n_total <- length(coherence_scores)
  n_organized <- sum(organization == "Organized")
  n_disorganized <- sum(organization == "Disorganized")
  
  # Score statistics
  mean_coherence <- mean(coherence_scores)
  median_coherence <- median(coherence_scores)
  score_range <- range(coherence_scores)
  
  # Organized vs disorganized comparison
  organized_scores <- coherence_scores[organization == "Organized"]
  disorganized_scores <- coherence_scores[organization == "Disorganized"]
  
  mean_organized <- mean(organized_scores)
  mean_disorganized <- mean(disorganized_scores)
  
  # Return summary
  list(
    total_ecotypes = n_total,
    organized_count = n_organized,
    organized_percentage = round(n_organized/n_total*100, 1),
    disorganized_count = n_disorganized,
    overall_mean_coherence = round(mean_coherence, 3),
    overall_median_coherence = round(median_coherence, 3),
    coherence_range = round(score_range, 3),
    mean_organized_coherence = round(mean_organized, 3),
    mean_disorganized_coherence = round(mean_disorganized, 3),
    organization_separation = round(mean_organized - mean_disorganized, 3)
  )
}

# Use with your results
summary_stats <- results_summary(your_results)
print(summary_stats)
```

### 📈 **Treatment Effect Interpretation**

#### **Effect Size Guidelines (Cohen's d)**

| Effect Size | Interpretation | Biological Significance |
|-------------|---------------|------------------------|
| **d ≥ 0.8** | Large effect | Major reorganization, strong treatment response |
| **0.5 ≤ d < 0.8** | Medium effect | Moderate reorganization, clear treatment impact |
| **0.2 ≤ d < 0.5** | Small effect | Subtle reorganization, mild treatment effect |
| **d < 0.2** | Negligible | No meaningful spatial reorganization |

#### **Treatment Effect Analysis**

```r
interpret_treatment_effects <- function(treatment_results, ecotype_annotations = NULL) {
  if (is.null(treatment_results)) {
    cat("No treatment effects to interpret.\n")
    return(NULL)
  }
  
  # Add effect magnitude classification
  treatment_results$effect_magnitude <- sapply(abs(treatment_results$effect_size), function(d) {
    if (d >= 0.8) return("Large")
    if (d >= 0.5) return("Medium")  
    if (d >= 0.2) return("Small")
    return("Negligible")
  })
  
  # Add direction
  treatment_results$direction <- ifelse(treatment_results$effect_size > 0, 
                                       "Increased Organization", 
                                       "Decreased Organization")
  
  # Sort by effect size magnitude
  treatment_results <- treatment_results[order(abs(treatment_results$effect_size), decreasing = TRUE), ]
  
  # Print interpretation
  cat("=== TREATMENT EFFECT INTERPRETATION ===\n")
  
  significant_effects <- treatment_results[treatment_results$p_value < 0.05, ]
  if (nrow(significant_effects) > 0) {
    cat("Significant treatment effects (p < 0.05):\n")
    for (i in 1:nrow(significant_effects)) {
      ecotype <- significant_effects$ecotype[i]
      annotation <- if (!is.null(ecotype_annotations)) ecotype_annotations[ecotype] else "Unknown"
      effect_size <- significant_effects$effect_size[i]
      magnitude <- significant_effects$effect_magnitude[i]
      direction <- significant_effects$direction[i]
      p_value <- significant_effects$p_value[i]
      
      cat(sprintf("- %s (%s): %s (%s effect, d=%.3f, p=%.4f)\n",
                  ecotype, annotation, direction, magnitude, abs(effect_size), p_value))
    }
  } else {
    cat("No significant treatment effects detected.\n")
  }
  
  # Effect magnitude summary
  cat(sprintf("\nEffect magnitude distribution:\n"))
  magnitude_table <- table(treatment_results$effect_magnitude)
  for (magnitude in names(magnitude_table)) {
    cat(sprintf("- %s effects: %d ecotypes\n", magnitude, magnitude_table[magnitude]))
  }
  
  return(treatment_results)
}

# Use with your results
if (!is.null(your_results$treatment_results)) {
  interpreted_effects <- interpret_treatment_effects(
    your_results$treatment_results,
    your_ecotype_annotations  # If you have them
  )
}
```

### 🔬 **Sample-Level Organization Patterns**

#### **Sample Classification**

```r
classify_samples <- function(results) {
  # Calculate mean coherence per sample
  if ("coherence_matrix" %in% names(results)) {
    if (is.data.frame(results$coherence_matrix)) {
      # Enhanced analysis format
      coherence_matrix <- results$coherence_matrix[, !names(results$coherence_matrix) %in% "metaprogram"]
      coherence_matrix <- as.matrix(coherence_matrix)
      sample_coherence <- colMeans(coherence_matrix, na.rm = TRUE)
    } else {
      # Basic analysis format
      sample_coherence <- colMeans(results$coherence_matrix, na.rm = TRUE)
    }
  } else {
    # Calculate from detailed results
    sample_coherence <- tapply(results$detailed_results$coherence_score, 
                              results$detailed_results$sample, mean, na.rm = TRUE)
  }
  
  # Classify samples
  threshold_high <- quantile(sample_coherence, 0.75, na.rm = TRUE)
  threshold_low <- quantile(sample_coherence, 0.25, na.rm = TRUE)
  
  sample_classification <- sapply(sample_coherence, function(score) {
    if (is.na(score)) return("Unknown")
    if (score >= threshold_high) return("Highly Organized")
    if (score <= threshold_low) return("Poorly Organized") 
    return("Moderately Organized")
  })
  
  # Summary
  cat("=== SAMPLE ORGANIZATION CLASSIFICATION ===\n")
  classification_table <- table(sample_classification)
  for (class in names(classification_table)) {
    cat(sprintf("- %s: %d samples (%.1f%%)\n", 
                class, classification_table[class],
                classification_table[class]/length(sample_classification)*100))
  }
  
  return(list(
    sample_coherence = sample_coherence,
    sample_classification = sample_classification,
    thresholds = c(low = threshold_low, high = threshold_high)
  ))
}

# Use with your results
sample_analysis <- classify_samples(your_results)
```

---

## 🚨 Troubleshooting

### ❌ **Common Errors and Solutions**

#### **Error: "threshold must be between 0 and 1"**
```r
# Problem: Invalid threshold value
coherence_threshold = 1.5  # ❌ WRONG

# Solution: Use valid range
coherence_threshold = 0.47  #  CORRECT (0-1 range)
```

#### **Error: "Ecotype column not found in metadata"**
```r
# Problem: Column name mismatch
ecotype_column = "wrong_column_name"  # ❌ WRONG

# Solution: Check your data
colnames(your_seurat@meta.data)  # See available columns
ecotype_column = "CompositionCluster_SE"  #  CORRECT
```

#### **Warning: "min_spots_per_ecotype < 10 may lead to unreliable results"**
```r
# Problem: Too few spots per ecotype
min_spots_per_ecotype = 5  # ❌ TOO LOW

# Solution: Check your data and adjust
ecotype_counts <- table(your_seurat@meta.data$ecotype_column)
print(ecotype_counts)  # See actual counts
min_spots_per_ecotype = 50  #  REASONABLE
```

#### **Issue: "All ecotypes classified as disorganized"**
```r
# Problem: Threshold too high for your data
coherence_threshold = 0.8  # ❌ TOO STRICT

# Solution: Analyze your data distribution
distribution <- analyze_coherence_distribution(your_seurat, "ecotype_col", "sample_col")
print(distribution$summary)
print(distribution$suggested_thresholds)

# Use data-appropriate threshold
coherence_threshold = 0.40  #  ADJUSTED FOR DATA
```

#### **Issue: "No treatment effects detected"**
```r
# Check treatment group balance
table(your_seurat@meta.data$treatment_column)
# Ensure: 1) Sufficient samples per group (≥5)
#         2) Reasonable treatment contrast
#         3) Correct column name
```

### 🔧 **Performance Issues**

#### **Analysis Too Slow**
```r
# Speed up analysis
run_spatial_analysis(
  n_permutations = 50,        # Reduce from 200
  enhanced_analysis = FALSE,   # Use basic mode
  min_spots_per_ecotype = 100  # Filter small ecotypes
)
```

#### **Memory Issues with Large Datasets**
```r
# Optimize for large data
run_spatial_analysis(
  enhanced_analysis = FALSE,   # Less memory intensive
  save_intermediate_results = FALSE,
  create_plots = FALSE,       # Skip plotting initially
  verbose = FALSE             # Reduce output
)
```

### 📊 **Data Quality Issues**

#### **Sparse or Low-Quality Data**
```r
# Adjust parameters for low-quality data
run_spatial_analysis(
  coherence_threshold = 0.35,  # More lenient
  min_spots_per_ecotype = 30,  # Lower requirements
  n_neighbors = 4,            # Conservative neighborhood
  n_permutations = 50         # Faster analysis
)
```

#### **Missing Spatial Coordinates**
```r
# The tool will use simulated coordinates with warning
# For better results, ensure your Seurat object has proper spatial data:
# your_seurat@images should contain coordinate information
```

---

## 📄 Publication Guidelines

### 📖 **Methods Section Template**

```markdown
### Spatial Ecotype Organization Analysis

Spatial ecotype organization was quantified using SpatialCoherence v1.0.0 
(https://github.com/ateeq-khaliq/SpatialCoherence). The analysis measures 
spatial coherence of multicellular ecotypes using neighborhood-based algorithms 
with permutation testing for statistical validation.

**Parameters:** Spatial coherence was calculated using a 6-neighbor approach 
with [X] random permutations. Ecotypes were classified as "organized" 
(coherence ≥ [threshold]) or "disorganized" (coherence < [threshold]) based on 
[justification for threshold choice]. Only ecotypes with ≥[X] spots across 
≥[X] samples were included in the analysis.

**Threshold Selection:** The coherence threshold ([value]) was selected based on 
[data-driven analysis/biological reasoning/literature precedent]. Threshold 
sensitivity analysis was performed using values from [X] to [X] to ensure 
robustness of classification.

**Statistical Testing:** [If treatment analysis] Treatment effects on spatial 
organization were assessed using t-tests comparing coherence scores between 
treatment groups, with Cohen's d calculated for effect size estimation. 
P-values < 0.05 were considered statistically significant.

**Reproducibility:** All analyses used random seed [X] for reproducible results.
```

### 📊 **Results Section Template**

```markdown
### Spatial Organization Patterns

We analyzed [X] spatial ecotypes across [X] samples containing [X] total spots. 
Of the [X] ecotypes analyzed, [X] ([X]%) exhibited organized spatial patterns 
(coherence ≥ [threshold]), while [X] ([X]%) showed disorganized patterns.

Coherence scores ranged from [X] to [X] (mean = [X] ± [X] SD). Organized 
ecotypes showed significantly higher coherence scores compared to disorganized 
ecotypes (mean organized: [X] ± [X] vs. disorganized: [X] ± [X], p < 0.001).

**Most Organized Ecotypes:** [List top 3 with scores and functional annotations]

**Most Disorganized Ecotypes:** [List bottom 3 with scores and functional annotations]

**[If treatment analysis]** Treatment significantly altered spatial organization 
in [X] ecotypes (p < 0.05). [Treatment name] [increased/decreased] organization 
in [specific ecotypes with effect sizes].
```

### 📈 **Figure Legends**

#### **Main Coherence Figure**
```markdown
**Figure X. Spatial ecotype organization analysis.** 
(A) Spatial coherence scores for each ecotype, ranked from most to least 
organized. Bars colored by organization classification (green = organized, 
red = disorganized) using threshold = [X]. 
(B) Sample-level organization showing mean coherence across all ecotypes per sample. 
Points represent individual ecotypes within each sample; black bars show sample 
means ± SD. 
(C) [If applicable] Treatment effects on spatial organization. Effect sizes 
(Cohen's d) for each ecotype; positive values indicate increased organization 
with treatment. Solid bars indicate significant effects (p < 0.05).
```

#### **Supplementary Threshold Analysis**
```markdown
**Supplementary Figure X. Threshold sensitivity analysis.** 
Percentage of ecotypes classified as organized using different coherence 
thresholds. The selected threshold ([X], red dashed line) balances biological 
interpretability with data-driven considerations. Gray shaded region shows 
commonly used threshold range (0.40-0.55).
```

### 📋 **Data Availability Statement**

```markdown
**Data and Code Availability:** 
Spatial coherence analysis was performed using SpatialCoherence v1.0.0, 
available at https://github.com/ateeq-khaliq/SpatialCoherence. Analysis 
parameters and complete results are provided in Supplementary Table X. 
Code for reproducing all analyses is available at [your repository/supplement].
```

### 🏆 **Key Points for High-Impact Journals**

1. **Emphasize Discovery Approach:** "We employed an unbiased, data-driven approach to discover spatial organization patterns without predetermined assumptions about ecotype behavior."

2. **Highlight Methodological Rigor:** "Threshold selection was validated through sensitivity analysis and data distribution assessment."

3. **Demonstrate Reproducibility:** "All analyses used fixed random seeds and documented parameters for complete reproducibility."

4. **Show Biological Relevance:** Connect spatial organization patterns to known biology and functional consequences.

5. **Include Comprehensive Validation:** Show results are robust across multiple threshold values and parameter settings.

---

## 🎯 Best Practices Summary

###  **Do's**

1. **Always validate your data first**
   ```r
   validate_spatial_data(your_seurat, ecotype_column, sample_column)
   ```

2. **Understand your data distribution before choosing thresholds**
   ```r
   analyze_coherence_distribution(your_seurat, ecotype_column, sample_column)
   ```

3. **Perform threshold sensitivity analysis**
   ```r
   perform_threshold_analysis(your_seurat, ecotype_column, sample_column)
   ```

4. **Use reproducible settings for publication**
   ```r
   random_seed = 123  # Fixed seed
   n_permutations = 200  # Adequate precision
   ```

5. **Document your parameter choices**
   - Save analysis_parameters.csv
   - Justify threshold selection
   - Report sensitivity analysis

6. **Start with basic analysis, then enhance**
   ```r
   # First: enhanced_analysis = FALSE (exploration)
   # Then: enhanced_analysis = TRUE (publication)
   ```

### ❌ **Don'ts**

1. **Don't assume which ecotypes should be organized**
   - Let the tool discover patterns objectively

2. **Don't use arbitrary thresholds without justification**
   - Base threshold choice on data or biology

3. **Don't ignore quality control parameters**
   - Check min_spots_per_ecotype requirements
   - Ensure adequate sample sizes

4. **Don't skip parameter documentation**
   - Always save and report analysis settings

5. **Don't use only one threshold**
   - Include sensitivity analysis in publications

6. **Don't ignore warnings**
   - Address data quality issues appropriately

---

## 🔬 Advanced Applications

### 🧬 **Multi-Condition Comparisons**

```r
# Compare spatial organization across multiple conditions
conditions <- c("Control", "Treatment_A", "Treatment_B")
condition_results <- list()

for (condition in conditions) {
  condition_data <- subset(your_seurat, condition_column == condition)
  
  condition_results[[condition]] <- run_spatial_analysis(
    seurat_object = condition_data,
    ecotype_column = "ecotype_column",
    sample_column = "sample_column", 
    output_dir = paste0("Analysis_", condition),
    coherence_threshold = 0.47,  # Same threshold for comparison
    random_seed = 123,           # Same seed for reproducibility
    verbose = FALSE
  )
}

# Compare results across conditions
compare_conditions <- function(results_list) {
  # Extract mean coherence for each condition
  coherence_comparison <- data.frame()
  
  for (condition in names(results_list)) {
    condition_coherence <- results_list[[condition]]$mean_coherence
    temp_df <- data.frame(
      Condition = condition,
      Ecotype = names(condition_coherence),
      Coherence = condition_coherence
    )
    coherence_comparison <- rbind(coherence_comparison, temp_df)
  }
  
  return(coherence_comparison)
}

multi_condition_data <- compare_conditions(condition_results)
```

### 📊 **Longitudinal Analysis**

```r
# Analyze spatial organization changes over time
timepoints <- c("Baseline", "Week_4", "Week_8", "Week_12")
longitudinal_results <- list()

for (timepoint in timepoints) {
  timepoint_data <- subset(your_seurat, timepoint_column == timepoint)
  
  longitudinal_results[[timepoint]] <- run_spatial_analysis(
    seurat_object = timepoint_data,
    ecotype_column = "ecotype_column",
    sample_column = "sample_column",
    output_dir = paste0("Longitudinal_", timepoint),
    coherence_threshold = 0.47,
    random_seed = 123,
    verbose = FALSE
  )
}

# Track organization changes over time
track_organization_changes <- function(longitudinal_results) {
  # Create time series of coherence scores
  time_series_data <- data.frame()
  
  for (timepoint in names(longitudinal_results)) {
    timepoint_coherence <- longitudinal_results[[timepoint]]$mean_coherence
    temp_df <- data.frame(
      Timepoint = timepoint,
      Ecotype = names(timepoint_coherence),
      Coherence = timepoint_coherence
    )
    time_series_data <- rbind(time_series_data, temp_df)
  }
  
  return(time_series_data)
}

time_series_analysis <- track_organization_changes(longitudinal_results)
```

### 🔍 **Biomarker Discovery**

```r
# Identify spatial organization patterns associated with outcomes
discover_spatial_biomarkers <- function(results, clinical_data) {
  # Extract sample-level coherence scores
  sample_coherence <- classify_samples(results)$sample_coherence
  
  # Merge with clinical outcomes
  biomarker_data <- data.frame(
    Sample = names(sample_coherence),
    Mean_Coherence = sample_coherence
  )
  
  # Add clinical data
  biomarker_data <- merge(biomarker_data, clinical_data, by = "Sample")
  
  # Test association with outcomes
  if ("Survival_Status" %in% colnames(biomarker_data)) {
    # Survival analysis
    high_coherence <- biomarker_data$Mean_Coherence > median(biomarker_data$Mean_Coherence)
    biomarker_data$Coherence_Group <- ifelse(high_coherence, "High", "Low")
    
    # Test survival difference
    survival_test <- survival::survdiff(
      survival::Surv(Survival_Time, Survival_Status) ~ Coherence_Group, 
      data = biomarker_data
    )
    
    cat("Spatial organization biomarker analysis:\n")
    cat("High vs Low coherence survival difference: p =", 
        round(1 - pchisq(survival_test$chisq, 1), 4), "\n")
  }
  
  return(biomarker_data)
}

# Use with your clinical data
# biomarker_analysis <- discover_spatial_biomarkers(your_results, your_clinical_data)
```

---

## 📚 Additional Resources

### 🔗 **Links and References**

- **GitHub Repository:** https://github.com/ateeq-khaliq/SpatialCoherence
- **Issues and Support:** https://github.com/ateeq-khaliq/SpatialCoherence/issues
- **Documentation:** Complete function documentation in package help files

### 📖 **Recommended Reading**

1. **Spatial Transcriptomics Methods:**
   - Ståhl et al., Science 2016 (Visium methodology)
   - Moses & Pachter, Nature Reviews Genetics 2022 (Spatial transcriptomics review)

2. **Spatial Analysis in Cancer:**
   - Greenwald et al., Cell 2024 (Glioblastoma spatial organization)
   - Luca et al., Cell 2021 (EcoTyper multicellular programs)

3. **Statistical Methods:**
   - Cohen, Psychological Bulletin 1988 (Effect size interpretation)
   - Moran, Biometrika 1950 (Spatial autocorrelation principles)

### 💬 **Getting Help**

1. **Check this tutorial first** for parameter guidance
2. **Review troubleshooting section** for common issues  
3. **Search GitHub issues** for similar problems
4. **Create new GitHub issue** with:
   - Complete error message
   - Minimal reproducible example
   - Session info (`sessionInfo()`)
   - Parameter values used

### 🤝 **Contributing**

Contributions welcome! Please:
1. Fork the repository
2. Create feature branch
3. Add tests for new functionality
4. Submit pull request with clear description

---

## 📄 Citation

If you use SpatialCoherence in your research, please cite:

```bibtex
@software{khaliq2024spatialcoherence,
  title = {SpatialCoherence: Spatial Ecotype Organization Analysis for Spatial Transcriptomics},
  author = {Ateeq Khaliq},
  year = {2026},
  version = {1.0.0},
  url = {https://github.com/ateeq-khaliq/SpatialCoherence},
  note = {R package for discovering spatial organization patterns in multicellular ecotypes}
}
```

---

## 🎉 Conclusion

SpatialCoherence provides a powerful, objective framework for discovering spatial organization patterns in spatial transcriptomics data. By following this tutorial, you can:

 **Discover** which spatial ecotypes are organized vs disorganized in your data  
 **Analyze** treatment effects on spatial tissue architecture  
 **Generate** publication-quality results with proper statistical validation  
 **Interpret** findings in biological context with confidence  

**Happy analyzing! 🧬✨**

---

*SpatialCoherence v1.0.0 - Making spatial transcriptomics analysis accessible to everyone.*

**Author:** Ateeq Khaliq | **Email:** akhaliq@iu.edu | **Institution:** Indiana University
