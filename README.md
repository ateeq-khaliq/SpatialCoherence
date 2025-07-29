# SpatialCoherence

## Spatial Organization Analysis for Spatial Transcriptomics

[![R-CMD-check](https://github.com/ateeq-khaliq/SpatialCoherence/workflows/R-CMD-check/badge.svg)](https://github.com/ateeq-khaliq/SpatialCoherence/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/release/ateeq-khaliq/SpatialCoherence.svg)](https://github.com/ateeq-khaliq/SpatialCoherence/releases)
[![R version](https://img.shields.io/badge/R-%E2%89%A5%204.1.0-blue.svg)](https://www.r-project.org/)

**SpatialCoherence** is a comprehensive R package for analyzing spatial organization patterns in tissues using **spatial transcriptomics data**. Quantify tissue architecture, classify organized vs. disorganized regions, and analyze treatment effects on spatial patterns.

> **👨‍💻 Developed by:** [Ateeq Khaliq](https://scholar.google.com/citations?user=uciT_dkAAAAJ&hl=en) | Indiana University  
> **📧 Contact:** [akhaliq@iu.edu](mailto:akhaliq@iu.edu) | **🆔 ORCID:** [0000-0001-5200-081X](https://orcid.org/0000-0001-5200-081X)

---

## 🔬 **What is Spatial Coherence?**

Spatial coherence measures how spatially clustered different cell types/Clusters are within a tissue sample:
- **High coherence** (0.7-1.0) = Organized, structured tissue architecture
- **Medium coherence** (0.4-0.7) = Moderately organized regions  
- **Low coherence** (0.0-0.4) = Disorganized, disrupted tissue patterns

##  **Key Features**

| Feature | Description | Status |
|---------|-------------|--------|
| **Spatial Coherence Analysis** | Quantify how spatially organized cell types are | v1.0.0 |
| **Organization Classification** | Classify regions as organized vs. disorganized | v1.0.0 |
| **Treatment Effect Analysis** | Analyze therapy effects on spatial patterns | Coming soon |
| **Cell Type Composition** | Cell type distributions across organization levels | Coming soon |
| **Publication Plots** | Professional publication ready visualizations | v1.0.0 |
| **User Configuration** | Fully customizable via YAML files | v1.0.0 |

##  **Research Applications**

<details>
<summary><strong> Cancer Research</strong></summary>

- Analyze tumor architecture and organization
- Study treatment effects on tissue structure
- Identify therapy-resistant organized regions
- Compare spatial patterns across cancer types
- Quantify immune infiltration patterns
</details>

<details>
<summary><strong> Developmental Biology</strong></summary>

- Track tissue organization during development
- Analyze spatial gene expression patterns
- Study morphogenetic processes
- Quantify tissue patterning
</details>

<details>
<summary><strong> Immunology</strong></summary>

- Analyze immune cell spatial distributions
- Study immune infiltration patterns
- Examine spatial immune-tumor interactions
- Quantify immune organization
</details>

<details>
<summary><strong> General Biology</strong></summary>

- Quantify tissue architecture across conditions
- Study disease-related spatial changes
- Analyze organ-specific spatial patterns
- Compare spatial organization between samples
</details>

## 📋 **Supported Platforms**

| Platform | Type | Status | Notes |
|----------|------|--------|-------|
| **10X Visium** | Standard spatial transcriptomics | ✅ Tested | Most common platform |
| **Slide-seq** | High-resolution spatial | ✅ Compatible | High-density spots |
| **MERFISH** | Multiplexed imaging | ✅ Compatible | Single-cell resolution |
| **seqFISH+** | Sequential FISH | ✅ Compatible | Single-cell resolution |
| **Xenium** | 10X high-resolution | ✅ Compatible | Subcellular resolution |
| **CosMx** | NanoString spatial | ✅ Compatible | High-plex imaging |
| **Custom** | Any XY coordinate data | ✅ Flexible | User-defined formats |

##  **Quick Start (3 Steps)**

### **Step 1: Installation**
```r
# Install from GitHub
devtools::install_github("ateeq-khaliq/SpatialCoherence")
library(SpatialCoherence)
