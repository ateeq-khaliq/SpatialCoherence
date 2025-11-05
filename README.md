# SpatialCoherence

## Comprehensive Spatial Organization Analysis for Spatial Transcriptomics

[![R-CMD-check](https://github.com/ateeq-khaliq/SpatialCoherence/workflows/R-CMD-check/badge.svg)](https://github.com/ateeq-khaliq/SpatialCoherence/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**SpatialCoherence** is a comprehensive R package for analyzing spatial organization patterns in tissues using **spatial transcriptomics data**. Quantify tissue architecture, classify organized vs. disorganized regions, and analyze treatment effects on spatial patterns.

> **Developed by:** [Ateeq Khaliq](https://scholar.google.com/citations?user=uciT_dkAAAAJ&hl=en) | Indiana University  
> **📧 Contact:** [akhaliq@iu.edu](mailto:akhaliq@iu.edu) | **🆔 ORCID:** [0000-0001-5200-081X](https://orcid.org/0000-0001-5200-081X)

## 🚀 **Quick Start**

### **Installation**
```r
# Install from GitHub
devtools::install_github("ateeq-khaliq/SpatialCoherence")
library(SpatialCoherence)
```

### **Basic Usage**
```r
# Run comprehensive spatial analysis
results <- run_spatial_analysis(
  seurat_object = your_spatial_data,
  ecotype_column = "seurat_clusters",
  sample_column = "orig.ident", 
  save_csvs = TRUE,
  output_dir = "my_results"
)

# View results
print(results$mean_coherence)
print(results$organization_results)

# Access saved CSV files in output directory
```

## 🔬 **Key Features**

- ** Sample-wise Analysis**: Per-sample, per-ecotype coherence calculation
- ** CSV Export**: Automatic saving of all results to CSV files
- ** Treatment Analysis**: Compare spatial organization across conditions
- ** Publication Plots**: Professional visualization suite
- ** Multiple Platforms**: Visium, Slide-seq, MERFISH, seqFISH+ support

## 📊 **Output Files**

The package automatically generates:
- `detailed_coherence_results.csv` - Per sample, per ecotype results
- `mean_coherence_by_ecotype.csv` - Average coherence by ecotype  
- `coherence_matrix_samples_x_ecotypes.csv` - Full sample × ecotype matrix
- `treatment_effects.csv` - Treatment effect analysis (if applicable)

## 📋 **Platform Compatibility**

| Platform | Spot Resolution | Status | Features |
|----------|----------------|--------|----------|
| **10X Visium** | 55μm |  Tested | Full analysis suite |
| **Slide-seq** | 10μm |  Compatible | High-resolution analysis |
| **MERFISH** | Single-cell |  Compatible | Advanced neighbor detection |
| **seqFISH+** | Single-cell |  Compatible | Multi-scale analysis |
| **Xenium** | Subcellular |  Compatible | Ultra-high resolution |

## 🎯 **Citation**

If you use SpatialCoherence in your research, please cite:

```
Khaliq, A. (2024). SpatialCoherence: Spatial Organization Analysis for 
Spatial Transcriptomics. R package version 1.0.0. 
https://github.com/ateeq-khaliq/SpatialCoherence
```

## 📚 **References**

Based on methodology from:
Tirosh, I., et al. (2024). Integrative spatial analysis reveals a multi-layered organization of glioblastoma.

## 🔧 **Support**

- **GitHub Issues**: [Report bugs or request features](https://github.com/ateeq-khaliq/SpatialCoherence/issues)
- **Email**: akhaliq@iu.edu
- **Documentation**: See function help files (`?run_spatial_analysis`)

---

**SpatialCoherence v1.0.0** - Professional spatial transcriptomics analysis made simple.
