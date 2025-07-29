# SpatialCoherence 1.0.0

## Initial Release (2024-12-19)

This is the initial release of SpatialCoherence, an R package for analyzing spatial organization patterns in spatial transcriptomics data.

###  New Features

* **Core Analysis Framework**
  - Main analysis function `run_spatial_analysis()` with comprehensive workflow
  - Spatial coherence calculation adapted from Tirosh et al. methodology
  - Organization classification (organized vs disorganized tissue patterns)
  - Treatment effect analysis capabilities
  - Configurable parameters via YAML files

* **Data Validation and Quality Control**
  - `validate_spatial_data()` function for comprehensive input validation
  - Automatic detection of data quality issues
  - Helpful error messages and warnings
  - Support for multiple spatial transcriptomics platforms

* **Utility Functions**
  - `calculate_spatial_coherence()` for single ecotype analysis
  - `classify_organization()` with multiple classification methods
  - `load_config()` for configuration file management
  - `plot_organization_spectrum()` for Tirosh-style visualizations

* **Configuration System**
  - Comprehensive YAML configuration template
  - User-friendly parameter customization
  - Multiple analysis scenarios supported
  - Example configurations for common use cases

* **Platform Support**
  - 10X Genomics Visium
  - 10X Genomics Xenium  
  - Slide-seq
  - MERFISH
  - seqFISH
  - Any spatial transcriptomics data in Seurat format

* **Testing Infrastructure**
  - Comprehensive test suite using testthat
  - Automated testing via GitHub Actions
  - Cross-platform compatibility (Windows, macOS, Linux)
  - Input validation testing

### 📊 Analysis Capabilities

* **Spatial Coherence Analysis**
  - Calculate spatial coherence scores for each cell type/ecotype
  - Null model comparison using random permutations
  - Customizable neighbor definitions and distance thresholds
  - Abundance filtering to focus on relevant cell types

* **Organization Classification**
  - Threshold-based classification (default: 0.47)
  - Tertile-based classification for balanced groups
  - Quantile-based classification options
  - Tirosh-style organization spectrum plots

* **Treatment Analysis**
  - Compare spatial organization across treatment conditions
  - Statistical testing for treatment effects
  - Visualization of treatment-organization interactions
  - Support for multiple treatment groups

### 🔧 Technical Features

* **Robust Input Handling**
  - Comprehensive data validation
  - Informative error messages
  - Graceful handling of edge cases
  - Memory-efficient processing

* **Flexible Configuration**
  - YAML-based parameter files
  - Column name mapping for different datasets
  - Threshold customization
  - Output directory management

### 👨‍💻 Author

**Ateeq Khaliq**
- Email: akhaliq@iu.edu
- Institution: Indiana University
- ORCID: 0000-0001-5200-081X
- GitHub: @ateeq-khaliq

### 📄 License

MIT License - see LICENSE file for details

### 🚀 Getting Started

```r
# Install from GitHub
# devtools::install_github("ateeq-khaliq/SpatialCoherence")

# Load package
library(SpatialCoherence)

# Run basic analysis
# results <- run_spatial_analysis(your_spatial_data)

# View results
# print(results$coherence_results)
