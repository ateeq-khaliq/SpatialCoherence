# 🧬 SpatialCoherence Tutorial & Vignettes

<div align="center">

![SpatialCoherence Logo](https://img.shields.io/badge/SpatialCoherence-v1.0.0-brightgreen?style=for-the-badge&logo=r&logoColor=white)

**Complete Guide to Spatial Organization Analysis in Spatial Transcriptomics**

[![GitHub](https://img.shields.io/badge/GitHub-ateeq--khaliq%2FSpatialCoherence-blue?style=flat-square&logo=github)](https://github.com/ateeq-khaliq/SpatialCoherence)
[![License](https://img.shields.io/badge/License-MIT-yellow?style=flat-square)](https://opensource.org/licenses/MIT)
[![R Version](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue?style=flat-square&logo=r)](https://www.r-project.org/)

*Master spatial coherence analysis with step-by-step tutorials, real-world examples, and best practices*

</div>

---

## 📚 **Tutorial Contents**

<table>
<tr>
<td align="center" width="33%">

### 🚀 **Getting Started**
![Beginner](https://img.shields.io/badge/Level-Beginner-green?style=flat-square)

Perfect for new users
- Installation & Setup
- Basic Concepts
- First Analysis

</td>
<td align="center" width="33%">

### 🔬 **Core Analysis**
![Intermediate](https://img.shields.io/badge/Level-Intermediate-orange?style=flat-square)

Master the fundamentals
- Spatial Coherence Calculation
- Organization Classification
- Parameter Selection

</td>
<td align="center" width="33%">

### 🧠 **Advanced Applications**
![Advanced](https://img.shields.io/badge/Level-Advanced-red?style=flat-square)

Real-world research
- Cancer Studies
- Treatment Analysis
- Custom Workflows

</td>
</tr>
</table>

---

## 🎯 **Quick Navigation**

<div align="center">

| 📖 Tutorial | ⏱️ Duration | 🎓 Level | 📋 Prerequisites |
|-------------|------------|----------|------------------|
| [**🔧 Installation Guide**](#installation) | 5 min | Beginner | R ≥ 4.1.0 |
| [**⚡ Quick Start**](#quick-start) | 10 min | Beginner | Spatial data |
| [**📊 Complete Workflow**](#complete-workflow) | 30 min | Intermediate | Seurat object |
| [**🔬 Cancer Research Example**](#cancer-example) | 20 min | Intermediate | Basic R |
| [**🧪 Treatment Analysis**](#treatment-analysis) | 25 min | Advanced | Statistics knowledge |
| [**⚙️ Parameter Optimization**](#parameter-guide) | 15 min | Advanced | Domain expertise |

</div>

---

## 🌟 **What You'll Learn**

<details>
<summary><b>🎓 Core Concepts</b></summary>

- **Spatial Coherence**: How to measure if your clusters form spatial patterns
- **Organization Classification**: Distinguish organized vs. disorganized tissue regions  
- **Parameter Selection**: Choose optimal settings for your research
- **Result Interpretation**: Understand what your results mean biologically

</details>

<details>
<summary><b>🛠️ Practical Skills</b></summary>

- Calculate spatial coherence scores from your clusters
- Classify tissue organization patterns
- Create publication-ready visualizations
- Optimize analysis parameters for your data
- Troubleshoot common issues

</details>

<details>
<summary><b>🔬 Research Applications</b></summary>

- **Cancer Research**: Analyze tumor architecture and treatment effects
- **Developmental Biology**: Track tissue organization during development
- **Immunology**: Study immune cell spatial distributions
- **Drug Discovery**: Evaluate therapeutic effects on tissue structure

</details>

---

## <a id="installation"></a>🔧 **Installation Guide**

<div align="center">
<img src="https://img.shields.io/badge/Step%201-Install%20Package-blue?style=for-the-badge" />
</div>

### **Prerequisites**
```r
# Check R version (must be ≥ 4.1.0)
R.version.string

# Install required packages
if (!require(devtools)) install.packages("devtools")
if (!require(ggplot2)) install.packages("ggplot2")
```

### **Install SpatialCoherence**
```r
# Install from GitHub
devtools::install_github("ateeq-khaliq/SpatialCoherence")

# Load the package
library(SpatialCoherence)

# Verify installation
get_package_version()
```

<div align="center">
✅ <b>Success!</b> You should see: <code>[1] "1.0.0"</code>
</div>

---

## <a id="quick-start"></a>⚡ **Quick Start (10 Minutes)**

<div align="center">
<img src="https://img.shields.io/badge/Quick%20Start-10%20Minutes-green?style=for-the-badge" />
</div>

### **Step 1: Get Your Data Ready**
```r
# You need a Seurat object with:
# 1. Spatial coordinates
# 2. Cluster assignments (e.g., seurat_clusters)

# Check your clusters
cluster_column <- "seurat_clusters"  # Your cluster column name
# unique_clusters <- unique(your_seurat_object@meta.data[[cluster_column]])
# print(unique_clusters)
```

### **Step 2: Calculate Spatial Coherence**
```r
# Calculate how spatially organized each cluster is
# Replace with your actual Seurat object and cluster column

# For tutorial, we'll use example scores:
coherence_scores <- c(
  "0" = 0.72,  # Highly organized cluster
  "1" = 0.28,  # Dispersed cluster
  "2" = 0.55,  # Moderately organized
  "3" = 0.41   # Moderately dispersed
)

print(coherence_scores)
```

### **Step 3: Classify Organization**
```r
# Classify clusters as organized vs. disorganized
organization <- classify_organization(coherence_scores)
print(organization)
```

**Expected Output:**
```
         0          1          2          3 
"Organized" "Disorganized" "Organized" "Disorganized" 
```

### **Step 4: Visualize Results**
```r
# Create organization spectrum plot
plot <- plot_organization_spectrum(coherence_scores, organization)
print(plot)
```

<div align="center">
🎉 <b>Congratulations!</b> You've completed your first spatial coherence analysis!
</div>

---

## <a id="complete-workflow"></a>📊 **Complete Workflow (30 Minutes)**

<div align="center">
<img src="https://img.shields.io/badge/Complete%20Workflow-30%20Minutes-orange?style=for-the-badge" />
</div>

### **Phase 1: Data Preparation**

<details>
<summary><b>📋 Click to expand: Data Requirements</b></summary>

Your Seurat object must have:
```r
# 1. Expression data
# your_seurat_object@assays$RNA@counts

# 2. Metadata with clusters
# your_seurat_object@meta.data$seurat_clusters

# 3. Spatial coordinates  
# your_seurat_object@images

# Validate your data
validation <- validate_spatial_data(your_seurat_object)
print(validation)
```

</details>

### **Phase 2: Configuration Setup**

```r
# Load default configuration
config <- get_default_config()
print(config)

# Key parameters:
# - coherence_threshold: 0.47 (adjust based on your needs)
# - n_neighbors: 6 (spatial neighborhood size)
# - abundance_threshold: 0.05 (minimum cluster size)
```

### **Phase 3: Spatial Analysis**

```r
# Calculate coherence for all clusters
cluster_column <- "seurat_clusters"

# Method 1: Real calculation (with your data)
# coherence_scores <- calculate_spatial_coherence(your_seurat_object, cluster_column)

# Method 2: Example for tutorial
coherence_scores <- c(
  "0" = 0.78, "1" = 0.32, "2" = 0.65, "3" = 0.29, 
  "4" = 0.71, "5" = 0.44, "6" = 0.58, "7" = 0.36
)

# Classify organization
organization <- classify_organization(coherence_scores)

# Print results
cat("=== SPATIAL ORGANIZATION RESULTS ===\n")
for(cluster in names(coherence_scores)) {
  score <- coherence_scores[cluster]
  org <- organization[cluster]
  status <- ifelse(org == "Organized", "🟢", "🔴")
  cat(sprintf("%s Cluster %s: %s (%.2f)\n", status, cluster, org, score))
}
```

### **Phase 4: Advanced Analysis**

```r
# Test different thresholds
thresholds <- c(0.3, 0.47, 0.6)
cat("\n=== THRESHOLD SENSITIVITY ===\n")
for(thresh in thresholds) {
  result <- classify_organization(coherence_scores, threshold = thresh)
  organized_count <- sum(result == "Organized")
  cat(sprintf("Threshold %.2f: %d/%d organized\n", thresh, organized_count, length(result)))
}

# Summary statistics
cat("\n=== SUMMARY STATISTICS ===\n")
cat("Mean coherence:", round(mean(coherence_scores), 3), "\n")
cat("Range:", round(min(coherence_scores), 3), "-", round(max(coherence_scores), 3), "\n")
organized_pct <- sum(organization == "Organized") / length(organization) * 100
cat("Organized clusters:", round(organized_pct, 1), "%\n")
```

### **Phase 5: Visualization**

```r
# Create comprehensive visualization
library(ggplot2)

# Basic plot
plot1 <- plot_organization_spectrum(coherence_scores, organization)

# Enhanced plot
plot2 <- plot1 + 
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

print(plot2)

# Save plot
# ggsave("spatial_organization.png", plot2, width = 10, height = 6, dpi = 300)
```

---

## <a id="cancer-example"></a>🔬 **Cancer Research Example (20 Minutes)**

<div align="center">
<img src="https://img.shields.io/badge/Cancer%20Research-20%20Minutes-red?style=for-the-badge" />
</div>

### **Scenario: Tumor Spatial Architecture Analysis**

```r
# Simulate typical cancer tissue spatial organization
cancer_coherence <- c(
  "0" = 0.84,  # Dense tumor core - highly organized
  "1" = 0.23,  # Infiltrating T cells - dispersed
  "2" = 0.71,  # Stromal regions - organized
  "3" = 0.18,  # NK cells - highly dispersed  
  "4" = 0.67,  # Vasculature - organized
  "5" = 0.29,  # B cells - moderately dispersed
  "6" = 0.45,  # Mixed immune - borderline
  "7" = 0.52   # Macrophages - moderately organized
)

# Analyze organization
cancer_org <- classify_organization(cancer_coherence)

# Create detailed report
cat("🔬 CANCER TISSUE SPATIAL ANALYSIS\n")
cat("==================================\n")

# Group by organization status
organized <- names(cancer_org)[cancer_org == "Organized"]
disorganized <- names(cancer_org)[cancer_org == "Disorganized"]

cat("\n🟢 ORGANIZED CLUSTERS (Form distinct regions):\n")
for(cluster in organized) {
  score <- cancer_coherence[cluster]
  cat(sprintf("   Cluster %s: %.2f", cluster, score))
  if(score > 0.7) cat(" (Highly organized)")
  cat("\n")
}

cat("\n🔴 DISORGANIZED CLUSTERS (Dispersed throughout tissue):\n")
for(cluster in disorganized) {
  score <- cancer_coherence[cluster]
  cat(sprintf("   Cluster %s: %.2f", cluster, score))
  if(score < 0.3) cat(" (Highly dispersed)")
  cat("\n")
}

# Calculate cancer-specific metrics
tumor_clusters <- c("0")  # Assuming cluster 0 is tumor
immune_clusters <- c("1", "3", "5", "6")  # Immune clusters
stromal_clusters <- c("2", "4", "7")  # Stromal/support clusters

tumor_coherence <- mean(cancer_coherence[tumor_clusters])
immune_coherence <- mean(cancer_coherence[immune_clusters])
stromal_coherence <- mean(cancer_coherence[stromal_clusters])

cat("\n📊 TISSUE COMPARTMENT ANALYSIS:\n")
cat(sprintf("   Tumor coherence: %.2f %s\n", tumor_coherence, 
            ifelse(tumor_coherence > 0.6, "(Solid tumor mass)", "(Infiltrated tumor)")))
cat(sprintf("   Immune coherence: %.2f %s\n", immune_coherence,
            ifelse(immune_coherence < 0.4, "(Good infiltration)", "(Poor infiltration)")))
cat(sprintf("   Stromal coherence: %.2f %s\n", stromal_coherence,
            ifelse(stromal_coherence > 0.5, "(Organized stroma)", "(Disrupted stroma)")))
```

### **Clinical Interpretation**
```r
# Provide clinical context
cat("\n🏥 CLINICAL IMPLICATIONS:\n")

if(tumor_coherence > 0.7) {
  cat("• High tumor coherence suggests compact tumor mass\n")
  cat("• May indicate poor drug penetration\n")
  cat("• Consider combination therapies\n")
}

if(immune_coherence < 0.3) {
  cat("• Low immune coherence indicates good T cell infiltration\n") 
  cat("• Favorable for immunotherapy response\n")
  cat("• Monitor for immune activation markers\n")
}

if(stromal_coherence > 0.6) {
  cat("• Organized stroma may create barriers\n")
  cat("• Consider anti-fibrotic agents\n")
  cat("• May affect drug delivery\n")
}
```

---

## <a id="treatment-analysis"></a>🧪 **Treatment Analysis (25 Minutes)**

<div align="center">
<img src="https://img.shields.io/badge/Treatment%20Analysis-25%20Minutes-purple?style=for-the-badge" />
</div>

### **Scenario: Before vs. After Immunotherapy**

```r
# Pre-treatment spatial organization
pre_treatment <- c(
  "Tumor_core" = 0.87,     # Highly organized tumor
  "Tumor_edge" = 0.74,     # Organized tumor periphery
  "CD8_T_cells" = 0.19,    # Low T cell infiltration
  "CD4_T_cells" = 0.22,    # Low helper T cell presence
  "Tregs" = 0.45,          # Moderate regulatory T cells
  "Macrophages" = 0.58,    # Organized macrophages
  "B_cells" = 0.31,        # Dispersed B cells
  "NK_cells" = 0.15,       # Very dispersed NK cells
  "Stroma" = 0.69,         # Dense stromal barrier
  "Vasculature" = 0.72     # Organized blood vessels
)

# Post-treatment spatial organization  
post_treatment <- c(
  "Tumor_core" = 0.42,     # Disrupted tumor structure
  "Tumor_edge" = 0.38,     # Further disruption
  "CD8_T_cells" = 0.68,    # Increased T cell infiltration
  "CD4_T_cells" = 0.61,    # Activated helper T cells
  "Tregs" = 0.28,          # Reduced Tregs
  "Macrophages" = 0.34,    # Polarized macrophages
  "B_cells" = 0.71,        # Organized B cell response
  "NK_cells" = 0.55,       # Activated NK cells
  "Stroma" = 0.33,         # Disrupted stroma
  "Vasculature" = 0.48     # Remodeled vasculature
)

# Analyze treatment response
pre_org <- classify_organization(pre_treatment)
post_org <- classify_organization(post_treatment)

cat("💊 IMMUNOTHERAPY RESPONSE ANALYSIS\n")
cat("==================================\n")

cat("\n📅 PRE-TREATMENT:\n")
for(cell_type in names(pre_treatment)) {
  score <- pre_treatment[cell_type]
  org <- pre_org[cell_type]
  status <- ifelse(org == "Organized", "🟢", "🔴")
  cat(sprintf("%s %-12s: %s (%.2f)\n", status, cell_type, org, score))
}

cat("\n📅 POST-TREATMENT:\n")
for(cell_type in names(post_treatment)) {
  score <- post_treatment[cell_type]
  org <- post_org[cell_type]
  status <- ifelse(org == "Organized", "🟢", "🔴")
  cat(sprintf("%s %-12s: %s (%.2f)\n", status, cell_type, org, score))
}

# Calculate treatment response metrics
cat("\n📈 TREATMENT RESPONSE METRICS:\n")
for(cell_type in names(pre_treatment)) {
  pre_score <- pre_treatment[cell_type]
  post_score <- post_treatment[cell_type]
  change <- post_score - pre_score
  change_pct <- (change / pre_score) * 100
  
  direction <- ifelse(change > 0, "↗️", "↘️")
  magnitude <- ifelse(abs(change) > 0.2, "MAJOR", ifelse(abs(change) > 0.1, "moderate", "minor"))
  
  cat(sprintf("%s %-12s: %+.2f (%+.0f%%) [%s change]\n", 
              direction, cell_type, change, change_pct, magnitude))
}

# Response classification
tumor_disruption <- mean(c(pre_treatment["Tumor_core"] - post_treatment["Tumor_core"],
                          pre_treatment["Tumor_edge"] - post_treatment["Tumor_edge"]))
immune_activation <- mean(c(post_treatment["CD8_T_cells"] - pre_treatment["CD8_T_cells"],
                           post_treatment["NK_cells"] - pre_treatment["NK_cells"]))

cat("\n🎯 OVERALL RESPONSE:\n")
if(tumor_disruption > 0.3 && immune_activation > 0.3) {
  cat("✅ EXCELLENT RESPONSE: Significant tumor disruption + immune activation\n")
} else if(tumor_disruption > 0.2 || immune_activation > 0.2) {
  cat("✅ GOOD RESPONSE: Clear treatment effects observed\n")
} else {
  cat("⚠️  LIMITED RESPONSE: Consider treatment modification\n")
}
```

### **Predictive Biomarkers**
```r
# Identify potential predictive markers
cat("\n🔬 POTENTIAL BIOMARKERS:\n")

# High baseline immune infiltration predicts response
baseline_immune <- mean(pre_treatment[c("CD8_T_cells", "NK_cells")])
if(baseline_immune > 0.2) {
  cat("• High baseline immune infiltration (favorable)\n")
} else {
  cat("• Low baseline immune infiltration (may need priming)\n")
}

# Stromal barrier assessment
if(pre_treatment["Stroma"] > 0.6) {
  cat("• Dense stromal barrier detected\n")
  cat("• Consider anti-fibrotic combination\n")
}

# Vascular accessibility
if(pre_treatment["Vasculature"] > 0.7) {
  cat("• Organized vasculature suggests good drug access\n")
} else {
  cat("• Poor vascular organization may limit drug delivery\n")
}
```

---

## <a id="parameter-guide"></a>⚙️ **Parameter Optimization Guide (15 Minutes)**

<div align="center">
<img src="https://img.shields.io/badge/Parameter%20Guide-15%20Minutes-gold?style=for-the-badge" />
</div>

### **Coherence Threshold Selection**

```r
# Test threshold sensitivity with your data
test_data <- c(
  "A" = 0.65, "B" = 0.42, "C" = 0.78, "D" = 0.31, 
  "E" = 0.54, "F" = 0.29, "G" = 0.71, "H" = 0.38
)

cat("🎯 THRESHOLD OPTIMIZATION\n")
cat("========================\n")

thresholds <- seq(0.2, 0.8, by = 0.05)
results_table <- data.frame(
  Threshold = thresholds,
  Organized = numeric(length(thresholds)),
  Disorganized = numeric(length(thresholds)),
  Ratio = numeric(length(thresholds))
)

for(i in 1:length(thresholds)) {
  thresh <- thresholds[i]
  result <- classify_organization(test_data, threshold = thresh)
  organized_count <- sum(result == "Organized")
  disorganized_count <- sum(result == "Disorganized")
  
  results_table$Organized[i] <- organized_count
  results_table$Disorganized[i] <- disorganized_count
  results_table$Ratio[i] <- round(organized_count / disorganized_count, 2)
}

# Print optimization table
cat("Threshold | Organized | Disorganized | Ratio\n")
cat("----------|-----------|--------------|------\n")
for(i in 1:nrow(results_table)) {
  row <- results_table[i,]
  cat(sprintf("   %.2f   |     %d     |      %d       | %.2f\n", 
              row$Threshold, row$Organized, row$Disorganized, row$Ratio))
}

# Recommendations
cat("\n💡 RECOMMENDATIONS:\n")
cat("• Threshold 0.30-0.40: Use for highly mixed tissues\n")
cat("• Threshold 0.45-0.50: Balanced (recommended default)\n") 
cat("• Threshold 0.55-0.65: Use for well-organized tissues\n")
cat("• Threshold 0.70+: Very strict (only clear patterns)\n")
```

### **Sample Size Considerations**

```r
cat("\n📊 SAMPLE SIZE GUIDELINES\n")
cat("=========================\n")

# Simulate different sample sizes
sample_sizes <- c(50, 100, 250, 500, 1000)
cat("Spots per Cluster | Reliability | Recommendation\n")
cat("------------------|-------------|---------------\n")

for(size in sample_sizes) {
  reliability <- ifelse(size < 100, "Low", 
                       ifelse(size < 250, "Moderate",
                              ifelse(size < 500, "Good", "Excellent")))
  
  recommendation <- ifelse(size < 100, "Combine clusters",
                          ifelse(size < 250, "Proceed with caution",
                                 ifelse(size < 500, "Good for analysis", "Ideal")))
  
  cat(sprintf("       %4d       |    %-8s |  %s\n", size, reliability, recommendation))
}

cat("\n⚠️  IMPORTANT NOTES:\n")
cat("• Minimum 50 spots per cluster for basic analysis\n")
cat("• 100+ spots recommended for reliable results\n") 
cat("• 250+ spots ideal for statistical testing\n")
cat("• Consider biological significance over statistical significance\n")
```

---

## 🚨 **Troubleshooting Guide**

<details>
<summary><b>❌ Common Error: "Input must be a Seurat object"</b></summary>

**Problem:** Function doesn't recognize your object as Seurat format

**Solution:**
```r
# Check object class
class(your_object)

# Convert if needed
if(!inherits(your_object, "Seurat")) {
  # Your object needs to be a proper Seurat object
  stop("Please ensure your data is in Seurat format")
}
```

</details>

<details>
<summary><b>⚠️ Warning: "All clusters classified as disorganized"</b></summary>

**Problem:** Threshold too high or data quality issues

**Solutions:**
```r
# Solution 1: Lower threshold
result <- classify_organization(scores, threshold = 0.3)

# Solution 2: Check score distribution
summary(scores)
hist(scores, main = "Coherence Score Distribution")

# Solution 3: Verify biological expectations
# Are your clusters expected to be spatially organized?
```

</details>

<details>
<summary><b>🔧 Issue: "Non-numeric coherence scores"</b></summary>

**Problem:** Scores are not in numeric format

**Solution:**
```r
# Check score format
str(coherence_scores)

# Convert to numeric if needed
coherence_scores <- as.numeric(coherence_scores)
names(coherence_scores) <- cluster_names
```

</details>

---

## 📞 **Getting Help & Support**

<div align="center">

### 🤝 **Community Support**

[![GitHub Issues](https://img.shields.io/badge/GitHub-Issues-red?style=for-the-badge&logo=github)](https://github.com/ateeq-khaliq/SpatialCoherence/issues)
[![Email](https://img.shields.io/badge/Email-akhaliq%40iu.edu-blue?style=for-the-badge&logo=gmail)](mailto:akhaliq@iu.edu)

</div>

**Before asking for help, please:**
1. ✅ Check this tutorial for your specific use case
2. ✅ Review the troubleshooting guide above  
3. ✅ Try the example code to ensure basic functionality
4. ✅ Include a reproducible example in your question

**When reporting issues, include:**
- Your R version and package version
- Complete error message
- Minimal reproducible example
- Description of expected vs. actual behavior

---

## 📚 **Additional Resources**

### **📖 Recommended Reading**
- [Original Tirosh et al. Paper (Nature, 2016)](https://www.nature.com/articles/nature20123)
- [Spatial Transcriptomics Analysis Methods](https://www.nature.com/articles/s41592-022-01409-2)
- [Seurat Spatial Analysis Tutorial](https://satijalab.org/seurat/articles/spatial_vignette.html)

### **🔗 Related Tools**
- [Seurat](https://satijalab.org/seurat/) - Spatial transcriptomics analysis
- [SpatialDE](https://github.com/Teichlab/SpatialDE) - Spatial gene expression
- [RCTD](https://github.com/dmcable/spacexr) - Cell type deconvolution

### **📊 Example Datasets**
- [10X Visium Brain Data](https://www.10xgenomics.com/resources/datasets)
- [STdata Package](https://github.com/jbergenstrahle/STdata)
- [SpatialLIBD](http://spatial.libd.org/spatialLIBD/)

---

## 📄 **Citation**

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
> Tirosh, I., et al. (2016). Single-cell RNA-seq supports a developmental hierarchy in human oligodendroglioma. *Nature*, 539(7628), 309-313.

---

## 🏆 **Acknowledgments**

<div align="center">

**Developed by [Ateeq Khaliq](https://scholar.google.com/citations?user=uciT_dkAAAAJ&hl=en)**  
Indiana University • ORCID: [0000-0001-5200-081X](https://orcid.org/0000-0001-5200-081X)

*Special thanks to the spatial transcriptomics community and all package contributors*

</div>

---

<div align="center">

### ⭐ **Found this helpful? Give us a star on GitHub!** ⭐

[![Star on GitHub](https://img.shields.io/github/stars/ateeq-khaliq/SpatialCoherence?style=social)](https://github.com/ateeq-khaliq/SpatialCoherence)

**Happy analyzing! 🧬✨**

</div>
