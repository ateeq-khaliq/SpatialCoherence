# 🧬 SpatialCoherence

[![GitHub](https://img.shields.io/badge/GitHub-ateeq--khaliq%2FSpatialCoherence-blue?style=flat-square&logo=github)](https://github.com/ateeq-khaliq/SpatialCoherence)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue?style=flat-square&logo=r)](https://www.r-project.org/)
[![Version](https://img.shields.io/badge/Version-1.0.0-brightgreen?style=flat-square)](https://github.com/ateeq-khaliq/SpatialCoherence/releases)

**Spatial Organization Analysis for Spatial Transcriptomics Data**

SpatialCoherence is an R package for analyzing spatial organization patterns in spatial transcriptomics data. It quantifies how spatially coherent different cell clusters are within tissue sections and classifies regions as organized or disorganized based on their spatial patterns.

---

## 🚀 Key Features

- **Spatial Coherence Calculation**: Quantify spatial organization of cell clusters
- **Organization Classification**: Automatically classify tissue regions as organized vs. disorganized
- **Flexible Parameters**: Customizable thresholds and analysis parameters
- **Visualization Tools**: Create publication-ready plots and figures
- **Treatment Analysis**: Compare spatial organization before/after treatments
- **Cancer Research**: Specialized workflows for tumor architecture analysis

---

## 📦 Installation

### Prerequisites

Ensure you have R ≥ 4.1.0 installed:

```r
# Check R version
R.version.string
```

### Install Required Dependencies

```r
# Install required packages
if (!require(devtools)) install.packages("devtools")
if (!require(ggplot2)) install.packages("ggplot2")
if (!require(Seurat)) install.packages("Seurat")
```

### Install SpatialCoherence

```r
# Install from GitHub
devtools::install_github("ateeq-khaliq/SpatialCoherence")

# Load the package
library(SpatialCoherence)

# Verify installation
get_package_version()
#> [1] "1.0.0"
```

---

## ⚡ Quick Start

### Basic Usage

```r
library(SpatialCoherence)

# Load your Seurat object with spatial data
# seurat_obj <- readRDS("your_spatial_data.rds")

# Calculate spatial coherence for all clusters
cluster_column <- "seurat_clusters"
# coherence_scores <- calculate_spatial_coherence(seurat_obj, cluster_column)

# For demonstration, using example scores:
coherence_scores <- c(
  "0" = 0.72,  # Highly organized cluster
  "1" = 0.28,  # Dispersed cluster  
  "2" = 0.55,  # Moderately organized
  "3" = 0.41   # Moderately dispersed
)

# Classify organization patterns
organization <- classify_organization(coherence_scores)
print(organization)
#>         0          1          2          3 
#> "Organized" "Disorganized" "Organized" "Disorganized"

# Create visualization
plot <- plot_organization_spectrum(coherence_scores, organization)
print(plot)
```

### Expected Output

The analysis will classify your clusters into two categories:
- **Organized**: Clusters that form spatially coherent regions (coherence > 0.47)
- **Disorganized**: Clusters that are dispersed throughout the tissue (coherence ≤ 0.47)

---

## 📊 Complete Workflow

### 1. Data Preparation

Your Seurat object must contain:
- Expression data in `@assays$RNA@counts`
- Cluster assignments in `@meta.data` (e.g., "seurat_clusters")
- Spatial coordinates in `@images`

```r
# Validate your data
validation <- validate_spatial_data(seurat_obj)
print(validation)
```

### 2. Configuration Setup

```r
# Load default configuration
config <- get_default_config()
print(config)

# Key parameters:
# - coherence_threshold: 0.47 (default threshold)
# - n_neighbors: 6 (spatial neighborhood size)
# - abundance_threshold: 0.05 (minimum cluster size)
```

### 3. Spatial Analysis

```r
# Calculate coherence for all clusters
cluster_column <- "seurat_clusters"
coherence_scores <- calculate_spatial_coherence(seurat_obj, cluster_column)

# Classify organization
organization <- classify_organization(coherence_scores)

# Print detailed results
cat("=== SPATIAL ORGANIZATION RESULTS ===\n")
for(cluster in names(coherence_scores)) {
  score <- coherence_scores[cluster]
  org <- organization[cluster]
  status <- ifelse(org == "Organized", "🟢", "🔴")
  cat(sprintf("%s Cluster %s: %s (%.2f)\n", status, cluster, org, score))
}
```

### 4. Advanced Analysis

```r
# Test threshold sensitivity
thresholds <- c(0.3, 0.47, 0.6)
for(thresh in thresholds) {
  result <- classify_organization(coherence_scores, threshold = thresh)
  organized_count <- sum(result == "Organized")
  cat(sprintf("Threshold %.2f: %d/%d organized\n", 
              thresh, organized_count, length(result)))
}

# Summary statistics
cat("Mean coherence:", round(mean(coherence_scores), 3), "\n")
cat("Range:", round(min(coherence_scores), 3), "-", 
    round(max(coherence_scores), 3), "\n")
organized_pct <- sum(organization == "Organized") / length(organization) * 100
cat("Organized clusters:", round(organized_pct, 1), "%\n")
```

---

## 🔬 Research Applications

### Cancer Research Example

```r
# Analyze tumor spatial architecture
cancer_coherence <- c(
  "Tumor_core" = 0.84,     # Dense tumor mass
  "T_cells" = 0.23,        # Infiltrating immune cells
  "Stroma" = 0.71,         # Organized stromal regions
  "NK_cells" = 0.18,       # Dispersed NK cells
  "Vasculature" = 0.67,    # Organized blood vessels
  "B_cells" = 0.29,        # Dispersed B cells
  "Macrophages" = 0.52     # Moderately organized
)

# Analyze organization
cancer_org <- classify_organization(cancer_coherence)

# Calculate tissue-specific metrics
tumor_clusters <- c("Tumor_core")
immune_clusters <- c("T_cells", "NK_cells", "B_cells")
stromal_clusters <- c("Stroma", "Vasculature", "Macrophages")

tumor_coherence <- mean(cancer_coherence[tumor_clusters])
immune_coherence <- mean(cancer_coherence[immune_clusters])
stromal_coherence <- mean(cancer_coherence[stromal_clusters])

cat("Tumor coherence:", round(tumor_coherence, 2), "\n")
cat("Immune coherence:", round(immune_coherence, 2), "\n")  
cat("Stromal coherence:", round(stromal_coherence, 2), "\n")
```

### Treatment Response Analysis

```r
# Compare before/after treatment
pre_treatment <- c(
  "Tumor" = 0.87,      # Highly organized tumor
  "CD8_T" = 0.19,      # Low T cell infiltration
  "Stroma" = 0.69      # Dense stromal barrier
)

post_treatment <- c(
  "Tumor" = 0.42,      # Disrupted tumor structure
  "CD8_T" = 0.68,      # Increased T cell infiltration  
  "Stroma" = 0.33      # Disrupted stroma
)

# Calculate treatment response
for(cell_type in names(pre_treatment)) {
  change <- post_treatment[cell_type] - pre_treatment[cell_type]
  cat(sprintf("%s: %+.2f change\n", cell_type, change))
}
```

---

## ⚙️ Parameter Guidelines

### Coherence Threshold Selection

| Threshold | Use Case | Description |
|-----------|----------|-------------|
| 0.30-0.40 | Mixed tissues | Use for highly heterogeneous samples |
| 0.45-0.50 | Balanced (default) | Recommended for most applications |
| 0.55-0.65 | Organized tissues | Use for well-structured samples |
| 0.70+ | Very strict | Only clear spatial patterns |

### Sample Size Requirements

| Spots per Cluster | Reliability | Recommendation |
|------------------|-------------|----------------|
| 50 spots | Low | Combine clusters |
| 100 spots | Moderate | Proceed with caution |
| 250 spots | Good | Good for analysis |
| 500 spots | Excellent | Ideal for analysis |
| 1000+ spots | Excellent | Robust analysis |

**Important Notes:**
- Minimum 50 spots per cluster for basic analysis
- 100+ spots recommended for reliable results
- 250+ spots ideal for statistical testing
- Consider biological significance over statistical significance

---

## 📈 Visualization

### Basic Plot

```r
# Create organization spectrum plot
plot <- plot_organization_spectrum(coherence_scores, organization)
print(plot)
```

### Enhanced Visualization

```r
# Customize plot appearance
enhanced_plot <- plot_organization_spectrum(coherence_scores, organization) +
  labs(
    title = "Spatial Organization Analysis",
    subtitle = "Clusters ranked by spatial coherence",
    caption = "Green = Organized | Red = Disorganized"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Save plot
ggsave("spatial_organization.png", enhanced_plot, 
       width = 10, height = 6, dpi = 300)
```

---

## 🚨 Troubleshooting

### Common Issues

**Error: "Input must be a Seurat object"**
```r
# Check object class
class(your_object)

# Ensure proper Seurat object
if(!inherits(your_object, "Seurat")) {
  stop("Please ensure your data is in Seurat format")
}
```

**Warning: "All clusters classified as disorganized"**
```r
# Lower threshold
result <- classify_organization(scores, threshold = 0.3)

# Check score distribution  
summary(scores)
hist(scores, main = "Coherence Score Distribution")
```

**Issue: "Non-numeric coherence scores"**
```r
# Check and convert scores
str(coherence_scores)
coherence_scores <- as.numeric(coherence_scores)
names(coherence_scores) <- cluster_names
```

---

## 📚 API Reference

### Main Functions

#### `calculate_spatial_coherence(seurat_obj, cluster_column)`
Calculate spatial coherence scores for all clusters.

**Parameters:**
- `seurat_obj`: Seurat object with spatial data
- `cluster_column`: Name of cluster column in metadata

**Returns:** Named vector of coherence scores

#### `classify_organization(coherence_scores, threshold = 0.47)`
Classify clusters as organized or disorganized.

**Parameters:**
- `coherence_scores`: Named vector of coherence scores
- `threshold`: Classification threshold (default: 0.47)

**Returns:** Named vector of organization classifications

#### `plot_organization_spectrum(coherence_scores, organization)`
Create visualization of spatial organization results.

**Parameters:**
- `coherence_scores`: Named vector of coherence scores
- `organization`: Named vector of organization classifications

**Returns:** ggplot object

### Utility Functions

- `get_default_config()`: Get default analysis parameters
- `validate_spatial_data()`: Validate Seurat object for analysis
- `get_package_version()`: Get current package version

---

## 🏥 Clinical Applications

### Biomarker Discovery
- Identify spatial organization patterns associated with treatment response
- Characterize tumor architecture and immune infiltration
- Assess tissue remodeling during disease progression

### Drug Development
- Evaluate therapeutic effects on tissue organization
- Screen compounds for spatial reorganization effects
- Optimize combination therapy strategies

### Diagnostic Applications
- Classify tissue samples based on spatial patterns
- Predict treatment outcomes from baseline organization
- Monitor disease progression through spatial changes

---

## 📖 Examples and Tutorials

### Basic Tutorial
```r
# Load example data
data("example_spatial_data")

# Run complete analysis
results <- run_spatial_coherence_analysis(example_spatial_data)

# View results
summary(results)
plot(results)
```

### Advanced Workflows
```r
# Custom parameter optimization
optimal_params <- optimize_parameters(spatial_data, 
                                     param_grid = list(
                                       threshold = seq(0.3, 0.7, 0.1),
                                       n_neighbors = c(4, 6, 8)
                                     ))

# Batch processing multiple samples
batch_results <- process_multiple_samples(sample_list, 
                                         config = optimal_params)
```

---

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

### Development Setup
```bash
git clone https://github.com/ateeq-khaliq/SpatialCoherence.git
cd SpatialCoherence
```

### Running Tests
```r
devtools::test()
```

### Code Style
We follow the [tidyverse style guide](https://style.tidyverse.org/).

---

## 📄 Citation

If you use SpatialCoherence in your research, please cite:

```bibtex
@software{khaliq2024spatialcoherence,
  title = {SpatialCoherence: Spatial Organization Analysis for Spatial Transcriptomics},
  author = {Ateeq Khaliq},
  year = {2024},
  version = {1.0.0},
  url = {https://github.com/ateeq-khaliq/SpatialCoherence},
  doi = {10.5281/zenodo.XXXXXXX}
}
```

**Original Methodology:**
> Greenwald, A.C., Galili Darnell, N., Hoefflin, R., et al. (2024). Integrative spatial analysis reveals a multi-layered organization of glioblastoma. *Cell*, 187(10), 2485-2501.e26. https://doi.org/10.1016/j.cell.2024.04.020

---

## 📞 Support

### Getting Help

- 📧 **Email**: [akhaliq@iu.edu](mailto:akhaliq@iu.edu)
- 🐛 **Bug Reports**: [GitHub Issues](https://github.com/ateeq-khaliq/SpatialCoherence/issues)
- 💬 **Discussions**: [GitHub Discussions](https://github.com/ateeq-khaliq/SpatialCoherence/discussions)

### Before Asking for Help

1. ✅ Check the documentation and examples
2. ✅ Review the troubleshooting guide  
3. ✅ Try the example code to ensure basic functionality
4. ✅ Search existing issues on GitHub

### When Reporting Issues

Please include:
- Your R version and package version
- Complete error message
- Minimal reproducible example
- Description of expected vs. actual behavior

---

## 🔗 Related Resources

### Documentation
- [Package Documentation](https://ateeq-khaliq.github.io/SpatialCoherence/)
- [Tutorials and Vignettes](https://ateeq-khaliq.github.io/SpatialCoherence/articles/)
- [API Reference](https://ateeq-khaliq.github.io/SpatialCoherence/reference/)

### Related Tools
- [Seurat](https://satijalab.org/seurat/) - Single-cell analysis framework
- [SpatialDE](https://github.com/Teichlab/SpatialDE) - Spatial gene expression
- [RCTD](https://github.com/dmcable/spacexr) - Cell type deconvolution
- [SpatialLIBD](http://spatial.libd.org/spatialLIBD/) - Spatial transcriptomics datasets

### Example Datasets
- [10X Visium Datasets](https://www.10xgenomics.com/resources/datasets)
- [STdata Package](https://github.com/jbergenstrahle/STdata)
- [SpatialDB](https://www.spatialomics.org/SpatialDB/)

---

## 🏆 Acknowledgments

**Developed by [Ateeq Khaliq](https://scholar.google.com/citations?user=uciT_dkAAAAJ&hl=en)**  
Indiana University • ORCID: [0000-0001-5200-081X](https://orcid.org/0000-0001-5200-081X)

*Special thanks to the spatial transcriptomics community and all package contributors.*

---

## 📊 Project Status

![GitHub last commit](https://img.shields.io/github/last-commit/ateeq-khaliq/SpatialCoherence)
![GitHub issues](https://img.shields.io/github/issues/ateeq-khaliq/SpatialCoherence)
![GitHub pull requests](https://img.shields.io/github/issues-pr/ateeq-khaliq/SpatialCoherence)

**Current Version**: 1.0.0  
**Status**: Active Development  
**Last Updated**: July 2025

---

### ⭐ **Found this helpful? Give us a star on GitHub!** ⭐

[![Star on GitHub](https://img.shields.io/github/stars/ateeq-khaliq/SpatialCoherence?style=social)](https://github.com/ateeq-khaliq/SpatialCoherence)

**Happy analyzing! 🧬✨**

---

*SpatialCoherence - Making spatial transcriptomics analysis accessible to everyone.*
