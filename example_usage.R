# ============================================================================
# SpatialCoherence - Example Usage Script
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script demonstrates how to use SpatialCoherence for spatial organization
# analysis of spatial transcriptomics data.

library(SpatialCoherence)

# ============================================================================
# EXAMPLE 1: BASIC ANALYSIS WITH DEFAULT PARAMETERS
# ============================================================================

cat("=== Example 1: Basic Analysis ===\n")

# Load your spatial transcriptomics data (replace with your file path)
# spatial_data <- readRDS("path/to/your/spatial_data.rds")

# For demonstration, we'll show the function call:
# results <- run_spatial_analysis(spatial_data)

cat("To run basic analysis:\n")
cat("results <- run_spatial_analysis(your_seurat_object)\n\n")

# ============================================================================
# EXAMPLE 2: ANALYSIS WITH CUSTOM CONFIGURATION
# ============================================================================

cat("=== Example 2: Custom Configuration ===\n")

# First, customize the config/user_parameters.yaml file for your data
# Then run analysis with your configuration:

cat("To run with custom config:\n")
cat("results <- run_spatial_analysis(\n")
cat("  seurat_object = your_data,\n")
cat("  config_path = 'config/user_parameters.yaml',\n")
cat("  output_dir = 'my_results'\n")
cat(")\n\n")

# ============================================================================
# EXAMPLE 3: QUICK SETUP FOR COMMON DATA TYPES
# ============================================================================

cat("=== Example 3: Quick Setup for Different Data Types ===\n")

# For 10X Visium data:
cat("For 10X Visium data:\n")
cat("- Set ecotype_column to 'seurat_clusters'\n")
cat("- Set sample_column to 'orig.ident'\n\n")

# For cancer treatment studies:
cat("For cancer treatment studies:\n")
cat("- Set ecotype_column to your cell type column\n")
cat("- Set treatment_column to your treatment column\n")
cat("- Define treated_values and untreated_values in config\n\n")

# ============================================================================
# EXAMPLE 4: ACCESSING AND INTERPRETING RESULTS
# ============================================================================

cat("=== Example 4: Working with Results ===\n")

cat("After running analysis, access results like this:\n\n")

cat("# View spatial coherence scores\n")
cat("print(results$coherence_results$coherence_scores)\n\n")

cat("# View organization classification\n") 
cat("print(results$organization_results$ecotype_organization)\n\n")

cat("# Check analysis metadata\n")
cat("print(results$metadata)\n\n")

cat("# View configuration used\n")
cat("print(results$config)\n\n")

# ============================================================================
# EXAMPLE 5: CUSTOMIZING ANALYSIS PARAMETERS
# ============================================================================

cat("=== Example 5: Customizing Parameters ===\n")

cat("Key parameters to customize in config/user_parameters.yaml:\n\n")

cat("1. Data structure:\n")
cat("   - ecotype_column: Your cell type column name\n")
cat("   - sample_column: Your sample identifier column\n") 
cat("   - treatment_column: Your treatment column (if applicable)\n\n")

cat("2. Analysis parameters:\n")
cat("   - coherence_threshold: Threshold for organized vs disorganized (default: 0.47)\n")
cat("   - n_neighbors: Number of spatial neighbors to consider (default: 6)\n")
cat("   - abundance_threshold: Minimum cell type abundance (default: 0.05)\n\n")

cat("3. Output settings:\n")
cat("   - output_dir: Where to save results\n")
cat("   - save_individual_plots: Whether to save individual plots\n\n")

# ============================================================================
# EXAMPLE 6: TROUBLESHOOTING COMMON ISSUES
# ============================================================================

cat("=== Example 6: Troubleshooting ===\n")

cat("Common issues and solutions:\n\n")

cat("1. 'Column not found' error:\n")
cat("   - Check your column names with: colnames(your_data@meta.data)\n")
cat("   - Update ecotype_column and sample_column in config accordingly\n\n")

cat("2. 'No spatial images found' error:\n")
cat("   - Ensure your Seurat object has spatial coordinates\n")
cat("   - Check with: length(your_data@images)\n\n")

cat("3. 'Configuration file not found' error:\n")
cat("   - Make sure config/user_parameters.yaml exists\n")
cat("   - Or use run_spatial_analysis(your_data) for defaults\n\n")

# ============================================================================
# EXAMPLE 7: BATCH PROCESSING MULTIPLE SAMPLES
# ============================================================================

cat("=== Example 7: Batch Processing ===\n")

cat("For multiple datasets:\n\n")

cat("datasets <- list(\n")
cat("  'sample1' = readRDS('sample1.rds'),\n")
cat("  'sample2' = readRDS('sample2.rds')\n")
cat(")\n\n")

cat("results_list <- list()\n")
cat("for(name in names(datasets)) {\n")
cat("  cat('Processing:', name, '\\n')\n")
cat("  results_list[[name]] <- run_spatial_analysis(\n")
cat("    seurat_object = datasets[[name]],\n")
cat("    output_dir = paste0('results_', name)\n")
cat("  )\n")
cat("}\n\n")

# ============================================================================
# GETTING HELP
# ============================================================================

cat("=== Getting Help ===\n")

cat("For help and support:\n")
cat("- Email: akhaliq@iu.edu\n")
cat("- GitHub Issues: https://github.com/ateeq-khaliq/SpatialCoherence/issues\n")
cat("- Documentation: Check the README.md file\n\n")

cat("=== Package Information ===\n")
cat("SpatialCoherence v1.0.0\n")
cat("Author: Ateeq Khaliq\n")
cat("Institution: Indiana University\n")
cat("ORCID: 0000-0001-5200-081X\n")
cat("License: MIT\n\n")

cat("Thank you for using SpatialCoherence!\n")
cat("This is the initial release - more features coming in future updates.\n")
