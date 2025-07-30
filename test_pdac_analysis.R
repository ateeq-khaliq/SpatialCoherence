
# ============================================================================
# TEST SPATIALCOHERENCE ON PDAC DATA
# ============================================================================

library(SpatialCoherence)

# Test with your PDAC data
cat("=== TESTING SPATIALCOHERENCE ON PDAC DATA ===\n")

# Example 1: Basic analysis with default threshold
cat("\n1. Basic analysis with default threshold (0.47):\n")
results_default <- run_spatial_analysis(
  seurat_object = pdac,
  ecotype_column = "CompositionCluster_SE",
  sample_column = "library_id",
  output_dir = "PDAC_Test_Default",
  enhanced_analysis = TRUE,
  verbose = TRUE
)

# Example 2: Custom threshold analysis
cat("\n2. Analysis with custom threshold (0.55):\n")
results_custom <- run_spatial_analysis(
  seurat_object = pdac,
  ecotype_column = "CompositionCluster_SE", 
  sample_column = "library_id",
  output_dir = "PDAC_Test_Custom",
  coherence_threshold = 0.55,  # More stringent
  enhanced_analysis = TRUE,
  verbose = TRUE
)

# Example 3: Compare thresholds
cat("\n3. Threshold sensitivity analysis:\n")
threshold_analysis <- perform_threshold_analysis(
  pdac, "CompositionCluster_SE", "library_id",
  thresholds = c(0.35, 0.40, 0.45, 0.50, 0.55, 0.60)
)
print(threshold_analysis$threshold_analysis)

# Example 4: Treatment analysis (if nac_treatment column exists)
if("nac_treatment" %in% colnames(pdac@meta.data)) {
  cat("\n4. Treatment effect analysis:\n")
  results_treatment <- run_spatial_analysis(
    seurat_object = pdac,
    ecotype_column = "CompositionCluster_SE",
    sample_column = "library_id", 
    treatment_column = "nac_treatment",
    output_dir = "PDAC_Test_Treatment",
    enhanced_analysis = TRUE,
    analyze_treatments = TRUE,
    verbose = TRUE
  )
  
  if(!is.null(results_treatment$treatment_results)) {
    cat("Treatment effects found for", nrow(results_treatment$treatment_results), "ecotypes\n")
    print(results_treatment$treatment_results)
  }
}

# Example 5: Distribution analysis
cat("\n5. Coherence distribution analysis:\n")
distribution <- analyze_coherence_distribution(pdac, "CompositionCluster_SE", "library_id")
cat("Suggested thresholds based on your data:\n")
print(distribution$suggested_thresholds)

# Compare results
cat("\n=== COMPARISON OF RESULTS ===\n")
comparison <- compare_thresholds(
  list(results_default, results_custom),
  c(0.47, 0.55)
)
print(comparison)

cat("\n=== TEST COMPLETE ===\n")
cat("Check the following directories for results:\n")
cat("- PDAC_Test_Default/\n")
cat("- PDAC_Test_Custom/\n")
if("nac_treatment" %in% colnames(pdac@meta.data)) {
  cat("- PDAC_Test_Treatment/\n")
}

