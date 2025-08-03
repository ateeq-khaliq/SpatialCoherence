#!/usr/bin/env Rscript

# ============================================================================
# SpatialCoherence Data Loader - Reads CSV Files from Tool Output
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script loads SpatialCoherence tool output CSV files and prepares them
# for plotting functions. Use this BEFORE running the plotting scripts.

# Required libraries
required_packages <- c("utils", "dplyr")

# Function to check and install packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      install.packages(pkg, quiet = TRUE)
      library(pkg, character.only = TRUE, quietly = TRUE)
    }
  }
}

check_and_install(required_packages)

# ============================================================================
# MAIN DATA LOADING FUNCTION
# ============================================================================

#' Load SpatialCoherence CSV Output Files
#' 
#' @param input_dir Path to directory containing SpatialCoherence CSV output files
#' @param verbose Print loading status messages
#' @return List containing all loaded data in format expected by plotting functions
#' 
#' @examples
#' # Load data from tool output directory
#' results <- load_spatialcoherence_data("pdac_spatial_coherence_test")
#' 
#' # Then use with plotting functions
#' all_plots <- create_all_spatialcoherence_plots(results, output_dir = "plots")
#' 
load_spatialcoherence_data <- function(input_dir, verbose = TRUE) {
  
  if (verbose) {
    cat("============================================================================\n")
    cat("Loading SpatialCoherence Data from CSV Files\n")
    cat("============================================================================\n\n")
    cat("Input directory:", input_dir, "\n\n")
  }
  
  # Check if directory exists
  if (!dir.exists(input_dir)) {
    stop("Input directory does not exist: ", input_dir)
  }
  
  # Initialize results list
  results <- list()
  
  # ========================================================================
  # REQUIRED FILES
  # ========================================================================
  
  # 1. Mean Coherence by Ecotype (REQUIRED)
  mean_coherence_file <- file.path(input_dir, "mean_coherence_by_ecotype.csv")
  if (file.exists(mean_coherence_file)) {
    if (verbose) cat("Loading mean_coherence_by_ecotype.csv...\n")
    
    mean_coherence_df <- read.csv(mean_coherence_file, stringsAsFactors = FALSE)
    
    # Extract mean coherence as named vector
    results$mean_coherence <- setNames(mean_coherence_df$mean_coherence, 
                                      mean_coherence_df$ecotype)
    
    # Extract organization results as named vector
    if ("organization" %in% colnames(mean_coherence_df)) {
      results$organization_results <- setNames(mean_coherence_df$organization, 
                                              mean_coherence_df$ecotype)
    }
    
    # Extract ecotype annotations if available
    if ("functional_annotation" %in% colnames(mean_coherence_df)) {
      annotations <- mean_coherence_df$functional_annotation
      names(annotations) <- mean_coherence_df$ecotype
      # Remove NA annotations
      annotations <- annotations[!is.na(annotations)]
      if (length(annotations) > 0) {
        results$ecotype_annotations <- annotations
      }
    }
    
    if (verbose) cat("  ✓ Loaded mean coherence for", length(results$mean_coherence), "ecotypes\n")
    
  } else {
    stop("Required file not found: mean_coherence_by_ecotype.csv")
  }
  
  # 2. Detailed Coherence Results (REQUIRED for most plots)
  detailed_file <- file.path(input_dir, "detailed_coherence_results.csv")
  if (file.exists(detailed_file)) {
    if (verbose) cat("Loading detailed_coherence_results.csv...\n")
    
    results$detailed_results <- read.csv(detailed_file, stringsAsFactors = FALSE)
    
    # Standardize column names
    if ("metaprogram" %in% colnames(results$detailed_results)) {
      # Already correct format
    } else if ("ecotype" %in% colnames(results$detailed_results)) {
      colnames(results$detailed_results)[colnames(results$detailed_results) == "ecotype"] <- "metaprogram"
    }
    
    if ("coherence_score" %in% colnames(results$detailed_results)) {
      # Already correct format
    } else if ("spatial_score" %in% colnames(results$detailed_results)) {
      colnames(results$detailed_results)[colnames(results$detailed_results) == "spatial_score"] <- "coherence_score"
    }
    
    if (verbose) cat("  ✓ Loaded detailed results:", nrow(results$detailed_results), "observations\n")
    
  } else {
    if (verbose) cat("  ! detailed_coherence_results.csv not found (some plots may not work)\n")
  }
  
  # ========================================================================
  # OPTIONAL FILES
  # ========================================================================
  
  # 3. Enhanced Coherence Matrix
  matrix_file <- file.path(input_dir, "enhanced_coherence_matrix.csv")
  if (file.exists(matrix_file)) {
    if (verbose) cat("Loading enhanced_coherence_matrix.csv...\n")
    
    coherence_matrix_df <- read.csv(matrix_file, stringsAsFactors = FALSE)
    
    # Check if this is the enhanced format with metaprogram column
    if ("metaprogram" %in% colnames(coherence_matrix_df)) {
      results$coherence_matrix <- coherence_matrix_df
      if (verbose) cat("  ✓ Loaded enhanced coherence matrix\n")
    } else {
      # Handle case where first column contains row names but no header
      # Make first column row names and remove it
      if (ncol(coherence_matrix_df) > 1) {
        # Set row names using first column, but handle duplicates
        row_names <- make.unique(as.character(coherence_matrix_df[, 1]))
        matrix_data <- coherence_matrix_df[, -1]
        rownames(matrix_data) <- row_names
        results$coherence_matrix <- as.matrix(matrix_data)
        if (verbose) cat("  ✓ Loaded coherence matrix (handled duplicate names)\n")
      }
    }
    
  } else {
    # Try alternative matrix file name
    alt_matrix_file <- file.path(input_dir, "coherence_matrix_samples_x_ecotypes.csv")
    if (file.exists(alt_matrix_file)) {
      if (verbose) cat("Loading coherence_matrix_samples_x_ecotypes.csv...\n")
      
      coherence_matrix_df <- read.csv(alt_matrix_file, stringsAsFactors = FALSE)
      
      # Handle potential duplicate row names
      if (ncol(coherence_matrix_df) > 1) {
        row_names <- make.unique(as.character(coherence_matrix_df[, 1]))
        matrix_data <- coherence_matrix_df[, -1]
        rownames(matrix_data) <- row_names
        results$coherence_matrix <- as.matrix(matrix_data)
        if (verbose) cat("  ✓ Loaded basic coherence matrix (handled duplicate names)\n")
      }
    } else {
      if (verbose) cat("  ! No coherence matrix found (heatmaps may not work)\n")
    }
  }
  
  # 4. Treatment Effects
  treatment_file <- file.path(input_dir, "treatment_effects.csv")
  if (file.exists(treatment_file)) {
    if (verbose) cat("Loading treatment_effects.csv...\n")
    
    results$treatment_results <- read.csv(treatment_file, stringsAsFactors = FALSE)
    
    if (verbose) cat("  ✓ Loaded treatment effects for", nrow(results$treatment_results), "ecotypes\n")
    
  } else {
    if (verbose) cat("  ! treatment_effects.csv not found (treatment plots will be skipped)\n")
  }
  
  # 5. Compositional Differences
  comp_file <- file.path(input_dir, "compositional_differences_data.csv")
  if (file.exists(comp_file)) {
    if (verbose) cat("Loading compositional_differences_data.csv...\n")
    
    comp_data <- read.csv(comp_file, stringsAsFactors = FALSE)
    results$compositional_results <- list(differences_data = comp_data)
    
    if (verbose) cat("  ✓ Loaded compositional differences\n")
    
  } else {
    if (verbose) cat("  ! compositional_differences_data.csv not found (compositional plots will be skipped)\n")
  }
  
  # 6. Ecotype Annotations (standalone file)
  annot_file <- file.path(input_dir, "ecotype_annotations.csv")
  if (file.exists(annot_file)) {
    if (verbose) cat("Loading ecotype_annotations.csv...\n")
    
    annot_df <- read.csv(annot_file, stringsAsFactors = FALSE)
    if ("ecotype" %in% colnames(annot_df) && "functional_annotation" %in% colnames(annot_df)) {
      annotations <- setNames(annot_df$functional_annotation, annot_df$ecotype)
      # Remove NA annotations
      annotations <- annotations[!is.na(annotations)]
      if (length(annotations) > 0) {
        results$ecotype_annotations <- annotations
      }
    }
    
    if (verbose) cat("  ✓ Loaded ecotype annotations\n")
    
  } else {
    if (verbose) cat("  ! ecotype_annotations.csv not found (functional plots will use generic names)\n")
  }
  
  # 7. Analysis Parameters
  params_file <- file.path(input_dir, "analysis_parameters.csv")
  if (file.exists(params_file)) {
    if (verbose) cat("Loading analysis_parameters.csv...\n")
    
    params_df <- read.csv(params_file, stringsAsFactors = FALSE)
    
    # Convert to named list
    analysis_params <- setNames(as.list(params_df$Value), params_df$Parameter)
    
    # Convert numeric parameters
    numeric_params <- c("coherence_threshold", "min_spots_per_ecotype", "n_neighbors", "n_permutations")
    for (param in numeric_params) {
      if (param %in% names(analysis_params)) {
        analysis_params[[param]] <- as.numeric(analysis_params[[param]])
      }
    }
    
    # Convert logical parameters
    logical_params <- c("enhanced_analysis")
    for (param in logical_params) {
      if (param %in% names(analysis_params)) {
        analysis_params[[param]] <- as.logical(analysis_params[[param]])
      }
    }
    
    results$analysis_parameters <- analysis_params
    
    if (verbose) cat("  ✓ Loaded analysis parameters\n")
    
  } else {
    if (verbose) cat("  ! analysis_parameters.csv not found (using default parameters)\n")
    # Set default parameters
    results$analysis_parameters <- list(
      coherence_threshold = 0.47,
      enhanced_analysis = FALSE
    )
  }
  
  # ========================================================================
  # DATA VALIDATION AND PREPARATION
  # ========================================================================
  
  if (verbose) cat("\nValidating loaded data...\n")
  
  # Ensure organization results exist
  if (is.null(results$organization_results) && !is.null(results$mean_coherence)) {
    threshold <- ifelse(!is.null(results$analysis_parameters$coherence_threshold), 
                       results$analysis_parameters$coherence_threshold, 0.47)
    results$organization_results <- ifelse(results$mean_coherence >= threshold, 
                                         "Organized", "Disorganized")
    names(results$organization_results) <- names(results$mean_coherence)
    if (verbose) cat("  ✓ Generated organization results using threshold =", threshold, "\n")
  }
  
  # Add metadata
  results$metadata <- list(
    n_ecotypes = length(results$mean_coherence),
    n_samples = if (!is.null(results$detailed_results)) length(unique(results$detailed_results$sample)) else NA,
    has_treatment_data = !is.null(results$treatment_results),
    has_compositional_data = !is.null(results$compositional_results),
    has_annotations = !is.null(results$ecotype_annotations),
    data_loaded_from = input_dir,
    loading_date = Sys.time()
  )
  
  # ========================================================================
  # SUMMARY
  # ========================================================================
  
  if (verbose) {
    cat("\n============================================================================\n")
    cat("Data Loading Complete\n")
    cat("============================================================================\n")
    cat("Summary:\n")
    cat("  • Ecotypes:", results$metadata$n_ecotypes, "\n")
    if (!is.na(results$metadata$n_samples)) {
      cat("  • Samples:", results$metadata$n_samples, "\n")
    }
    cat("  • Organization results:", ifelse(!is.null(results$organization_results), "✓", "✗"), "\n")
    cat("  • Detailed results:", ifelse(!is.null(results$detailed_results), "✓", "✗"), "\n")
    cat("  • Coherence matrix:", ifelse(!is.null(results$coherence_matrix), "✓", "✗"), "\n")
    cat("  • Treatment data:", ifelse(results$metadata$has_treatment_data, "✓", "✗"), "\n")
    cat("  • Compositional data:", ifelse(results$metadata$has_compositional_data, "✓", "✗"), "\n")
    cat("  • Functional annotations:", ifelse(results$metadata$has_annotations, "✓", "✗"), "\n")
    
    cat("\nReady for plotting! Use this results object with:\n")
    cat("  create_all_spatialcoherence_plots(results, output_dir = 'my_plots')\n\n")
  }
  
  return(results)
}

# ============================================================================
# CONVENIENCE FUNCTION FOR QUICK LOADING AND PLOTTING
# ============================================================================

#' Load Data and Create All Plots in One Step
#' 
#' @param input_dir Path to directory containing SpatialCoherence CSV files
#' @param output_dir Path where plots should be saved
#' @param treatment_column Name of treatment column in detailed_results (if any)
#' @param verbose Print status messages
#' 
quick_spatialcoherence_plots <- function(input_dir, 
                                        output_dir = "spatial_coherence_plots",
                                        treatment_column = "treatment",
                                        verbose = TRUE) {
  
  # Load the data
  results <- load_spatialcoherence_data(input_dir, verbose = verbose)
  
  # Check if plotting functions are loaded
  if (!exists("create_all_spatialcoherence_plots")) {
    cat("Loading plotting functions...\n")
    # Try to source the plotting scripts
    plot_scripts <- c("01_basic_coherence_plots.R", 
                     "02_treatment_analysis_plots.R",
                     "03_segmentation_analysis_plots.R", 
                     "04_heatmap_visualization_plots.R")
    
    for (script in plot_scripts) {
      if (file.exists(script)) {
        source(script)
      } else {
        cat("Warning: Could not find", script, "\n")
      }
    }
  }
  
  # Create plots if function is available
  if (exists("create_all_spatialcoherence_plots")) {
    all_plots <- create_all_spatialcoherence_plots(
      results = results,
      output_dir = output_dir,
      save_plots = TRUE,
      treatment_column = treatment_column
    )
    
    return(list(results = results, plots = all_plots))
  } else {
    cat("Plotting functions not available. Returning loaded data only.\n")
    return(results)
  }
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# Example usage (commented out):
# 
# # Method 1: Load data first, then create plots
# results <- load_spatialcoherence_data("pdac_spatial_coherence_test")
# 
# # Source plotting scripts
# source("01_basic_coherence_plots.R")
# source("02_treatment_analysis_plots.R") 
# source("03_segmentation_analysis_plots.R")
# source("04_heatmap_visualization_plots.R")
# 
# # Create all plots
# all_plots <- create_all_spatialcoherence_plots(
#   results = results,
#   output_dir = "my_visualization_output"
# )
# 
# # Method 2: Quick one-step approach
# output <- quick_spatialcoherence_plots(
#   input_dir = "pdac_spatial_coherence_test",
#   output_dir = "my_plots"
# )

cat("SpatialCoherence Data Loader loaded successfully.\n")
cat("Use load_spatialcoherence_data('your_csv_directory') to load your data.\n")
cat("Use quick_spatialcoherence_plots('input_dir', 'output_dir') for one-step plotting.\n")