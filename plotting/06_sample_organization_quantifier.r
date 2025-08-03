#!/usr/bin/env Rscript

# ============================================================================
# Sample Organization Quantifier - Rank Samples by Spatial Coherence
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script quantifies and ranks samples from most organized to most 
# disorganized based on spatial coherence scores, with comprehensive CSV outputs

# Required libraries
required_packages <- c("dplyr", "readr", "ggplot2", "scales")

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
# MAIN QUANTIFICATION FUNCTION
# ============================================================================

quantify_sample_organization <- function(results, 
                                        output_dir = "sample_organization_analysis",
                                        coherence_threshold = NULL,
                                        save_csvs = TRUE,
                                        generate_plot = TRUE,
                                        verbose = TRUE) {
  
  if (verbose) {
    cat("============================================================================\n")
    cat("SAMPLE ORGANIZATION QUANTIFICATION ANALYSIS\n")
    cat("============================================================================\n\n")
  }
  
  # Extract coherence threshold
  if (is.null(coherence_threshold)) {
    if (!is.null(results$analysis_parameters$coherence_threshold)) {
      coherence_threshold <- results$analysis_parameters$coherence_threshold
    } else {
      coherence_threshold <- 0.47  # Default PDAC threshold
    }
  }
  
  if (verbose) {
    cat("Using coherence threshold:", coherence_threshold, "\n\n")
  }
  
  # Create output directory
  if (save_csvs && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # ========================================================================
  # METHOD 1: FROM COHERENCE MATRIX (MOST COMPREHENSIVE)
  # ========================================================================
  
  sample_rankings <- NULL
  
  if (!is.null(results$coherence_matrix)) {
    if (verbose) cat("Method 1: Analyzing coherence matrix data...\n")
    
    # Handle different matrix formats
    if (is.data.frame(results$coherence_matrix) && "metaprogram" %in% colnames(results$coherence_matrix)) {
      # Enhanced format
      matrix_data <- results$coherence_matrix[, !names(results$coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- results$coherence_matrix$metaprogram
      ecotype_names <- results$coherence_matrix$metaprogram
    } else {
      # Standard matrix format
      matrix_data <- as.matrix(results$coherence_matrix)
      ecotype_names <- rownames(matrix_data)
    }
    
    # Calculate comprehensive sample statistics
    sample_rankings <- data.frame(
      sample_id = colnames(matrix_data),
      stringsAsFactors = FALSE
    )
    
    # Core coherence metrics
    sample_rankings$mean_coherence <- colMeans(matrix_data, na.rm = TRUE)
    sample_rankings$median_coherence <- apply(matrix_data, 2, median, na.rm = TRUE)
    sample_rankings$sd_coherence <- apply(matrix_data, 2, sd, na.rm = TRUE)
    sample_rankings$min_coherence <- apply(matrix_data, 2, min, na.rm = TRUE)
    sample_rankings$max_coherence <- apply(matrix_data, 2, max, na.rm = TRUE)
    sample_rankings$range_coherence <- sample_rankings$max_coherence - sample_rankings$min_coherence
    
    # Additional statistics
    sample_rankings$cv_coherence <- sample_rankings$sd_coherence / sample_rankings$mean_coherence  # Coefficient of variation
    sample_rankings$n_ecotypes_analyzed <- colSums(!is.na(matrix_data))
    sample_rankings$n_organized_ecotypes <- colSums(matrix_data >= coherence_threshold, na.rm = TRUE)
    sample_rankings$n_disorganized_ecotypes <- colSums(matrix_data < coherence_threshold, na.rm = TRUE)
    sample_rankings$prop_organized_ecotypes <- sample_rankings$n_organized_ecotypes / sample_rankings$n_ecotypes_analyzed
    
    # Quartile-based metrics
    sample_rankings$q25_coherence <- apply(matrix_data, 2, quantile, probs = 0.25, na.rm = TRUE)
    sample_rankings$q75_coherence <- apply(matrix_data, 2, quantile, probs = 0.75, na.rm = TRUE)
    sample_rankings$iqr_coherence <- sample_rankings$q75_coherence - sample_rankings$q25_coherence
    
    # Count ecotypes in different coherence ranges
    sample_rankings$n_very_high_coherence <- colSums(matrix_data >= 0.7, na.rm = TRUE)     # Very organized
    sample_rankings$n_high_coherence <- colSums(matrix_data >= 0.5 & matrix_data < 0.7, na.rm = TRUE)  # Moderately organized
    sample_rankings$n_medium_coherence <- colSums(matrix_data >= 0.3 & matrix_data < 0.5, na.rm = TRUE) # Intermediate
    sample_rankings$n_low_coherence <- colSums(matrix_data < 0.3, na.rm = TRUE)           # Disorganized
    
  } else if (!is.null(results$detailed_results)) {
    
    # ====================================================================
    # METHOD 2: FROM DETAILED RESULTS (ALTERNATIVE)
    # ====================================================================
    
    if (verbose) cat("Method 2: Analyzing detailed results data...\n")
    
    sample_rankings <- results$detailed_results %>%
      group_by(sample) %>%
      summarise(
        mean_coherence = mean(coherence_score, na.rm = TRUE),
        median_coherence = median(coherence_score, na.rm = TRUE),
        sd_coherence = sd(coherence_score, na.rm = TRUE),
        min_coherence = min(coherence_score, na.rm = TRUE),
        max_coherence = max(coherence_score, na.rm = TRUE),
        range_coherence = max_coherence - min_coherence,
        cv_coherence = sd_coherence / mean_coherence,
        n_ecotypes_analyzed = n(),
        n_organized_ecotypes = sum(coherence_score >= coherence_threshold, na.rm = TRUE),
        n_disorganized_ecotypes = sum(coherence_score < coherence_threshold, na.rm = TRUE),
        prop_organized_ecotypes = n_organized_ecotypes / n_ecotypes_analyzed,
        q25_coherence = quantile(coherence_score, 0.25, na.rm = TRUE),
        q75_coherence = quantile(coherence_score, 0.75, na.rm = TRUE),
        iqr_coherence = q75_coherence - q25_coherence,
        n_very_high_coherence = sum(coherence_score >= 0.7, na.rm = TRUE),
        n_high_coherence = sum(coherence_score >= 0.5 & coherence_score < 0.7, na.rm = TRUE),
        n_medium_coherence = sum(coherence_score >= 0.3 & coherence_score < 0.5, na.rm = TRUE),
        n_low_coherence = sum(coherence_score < 0.3, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      rename(sample_id = sample) %>%
      as.data.frame()
    
  } else {
    stop("No coherence matrix or detailed results available for analysis")
  }
  
  # ========================================================================
  # ORGANIZATION CLASSIFICATION AND RANKING
  # ========================================================================
  
  if (verbose) cat("Calculating organization classifications and rankings...\n")
  
  # Multi-level organization classification
  sample_rankings$organization_class <- case_when(
    sample_rankings$mean_coherence >= coherence_threshold + 0.10 ~ "Highly Organized",
    sample_rankings$mean_coherence >= coherence_threshold + 0.05 ~ "Moderately Organized", 
    sample_rankings$mean_coherence >= coherence_threshold ~ "Organized",
    sample_rankings$mean_coherence >= coherence_threshold - 0.05 ~ "Intermediate",
    sample_rankings$mean_coherence >= coherence_threshold - 0.10 ~ "Moderately Disorganized",
    TRUE ~ "Highly Disorganized"
  )
  
  # Binary classification
  sample_rankings$binary_organization <- ifelse(
    sample_rankings$mean_coherence >= coherence_threshold, "Organized", "Disorganized"
  )
  
  # Consistency score (low CV = more consistent organization)
  sample_rankings$organization_consistency <- case_when(
    sample_rankings$cv_coherence <= 0.2 ~ "Very Consistent",
    sample_rankings$cv_coherence <= 0.4 ~ "Consistent", 
    sample_rankings$cv_coherence <= 0.6 ~ "Moderately Variable",
    sample_rankings$cv_coherence <= 0.8 ~ "Variable",
    TRUE ~ "Highly Variable"
  )
  
  # Calculate rankings (1 = most organized)
  sample_rankings$rank_by_mean <- rank(-sample_rankings$mean_coherence, ties.method = "min")
  sample_rankings$rank_by_median <- rank(-sample_rankings$median_coherence, ties.method = "min")
  sample_rankings$rank_by_prop_organized <- rank(-sample_rankings$prop_organized_ecotypes, ties.method = "min")
  sample_rankings$rank_by_min <- rank(-sample_rankings$min_coherence, ties.method = "min")  # Worst-case ranking
  
  # Composite organization score (weighted average of multiple metrics)
  sample_rankings$composite_score <- (
    0.4 * scales::rescale(sample_rankings$mean_coherence, to = c(0, 100)) +
    0.2 * scales::rescale(sample_rankings$median_coherence, to = c(0, 100)) +
    0.2 * scales::rescale(sample_rankings$prop_organized_ecotypes, to = c(0, 100)) +
    0.1 * scales::rescale(sample_rankings$min_coherence, to = c(0, 100)) +
    0.1 * scales::rescale(-sample_rankings$cv_coherence, to = c(0, 100))  # Lower CV is better
  )
  
  sample_rankings$rank_composite <- rank(-sample_rankings$composite_score, ties.method = "min")
  
  # Calculate percentiles
  sample_rankings$percentile_rank <- percent_rank(sample_rankings$mean_coherence) * 100
  
  # Distance from threshold
  sample_rankings$distance_from_threshold <- sample_rankings$mean_coherence - coherence_threshold
  sample_rankings$abs_distance_from_threshold <- abs(sample_rankings$distance_from_threshold)
  
  # Sort by composite ranking (most organized first)
  sample_rankings <- sample_rankings[order(sample_rankings$rank_composite), ]
  
  # ========================================================================
  # ADD TREATMENT INFORMATION IF AVAILABLE
  # ========================================================================
  
  if (!is.null(results$detailed_results) && any(c("treatment", "nac_treatment") %in% colnames(results$detailed_results))) {
    if (verbose) cat("Adding treatment information...\n")
    
    treatment_col <- ifelse("treatment" %in% colnames(results$detailed_results), "treatment", "nac_treatment")
    
    treatment_info <- results$detailed_results %>%
      select(sample, !!sym(treatment_col)) %>%
      distinct() %>%
      rename(sample_id = sample, treatment_status = !!sym(treatment_col))
    
    sample_rankings <- merge(sample_rankings, treatment_info, by = "sample_id", all.x = TRUE)
    
    # Treatment-specific statistics
    if (any(!is.na(sample_rankings$treatment_status))) {
      treatment_summary <- sample_rankings %>%
        filter(!is.na(treatment_status)) %>%
        group_by(treatment_status) %>%
        summarise(
          n_samples = n(),
          mean_coherence_group = mean(mean_coherence, na.rm = TRUE),
          sd_coherence_group = sd(mean_coherence, na.rm = TRUE),
          median_coherence_group = median(mean_coherence, na.rm = TRUE),
          prop_organized = mean(binary_organization == "Organized", na.rm = TRUE),
          .groups = 'drop'
        )
      
      if (verbose) {
        cat("\nTreatment group summary:\n")
        print(treatment_summary)
      }
    }
  }
  
  # ========================================================================
  # SAVE CSV FILES
  # ========================================================================
  
  if (save_csvs) {
    if (verbose) cat("\nSaving CSV files...\n")
    
    # Main comprehensive ranking file
    readr::write_csv(sample_rankings, file.path(output_dir, "sample_organization_comprehensive_ranking.csv"))
    if (verbose) cat("  ✓ sample_organization_comprehensive_ranking.csv\n")
    
    # Simplified ranking (key metrics only)
    simplified_ranking <- sample_rankings %>%
      select(sample_id, rank_composite, mean_coherence, organization_class, binary_organization,
             prop_organized_ecotypes, n_ecotypes_analyzed, composite_score, percentile_rank,
             distance_from_threshold, organization_consistency) %>%
      arrange(rank_composite)
    
    readr::write_csv(simplified_ranking, file.path(output_dir, "sample_organization_simplified_ranking.csv"))
    if (verbose) cat("  ✓ sample_organization_simplified_ranking.csv\n")
    
    # Organization class summary
    class_summary <- sample_rankings %>%
      group_by(organization_class) %>%
      summarise(
        n_samples = n(),
        mean_coherence = mean(mean_coherence, na.rm = TRUE),
        sd_coherence = sd(mean_coherence, na.rm = TRUE),
        min_coherence = min(mean_coherence, na.rm = TRUE),
        max_coherence = max(mean_coherence, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(desc(mean_coherence))
    
    readr::write_csv(class_summary, file.path(output_dir, "organization_class_summary.csv"))
    if (verbose) cat("  ✓ organization_class_summary.csv\n")
    
    # Top and bottom samples
    top_organized <- head(sample_rankings, 10)
    bottom_organized <- tail(sample_rankings, 10)
    
    readr::write_csv(top_organized, file.path(output_dir, "top_10_most_organized_samples.csv"))
    readr::write_csv(bottom_organized, file.path(output_dir, "top_10_most_disorganized_samples.csv"))
    if (verbose) cat("  ✓ top_10_most_organized_samples.csv\n")
    if (verbose) cat("  ✓ top_10_most_disorganized_samples.csv\n")
    
    # Quartile-based groups
    quartile_groups <- sample_rankings %>%
      mutate(
        quartile = case_when(
          percentile_rank >= 75 ~ "Q4 (Most Organized)",
          percentile_rank >= 50 ~ "Q3 (Moderately Organized)",
          percentile_rank >= 25 ~ "Q2 (Moderately Disorganized)",
          TRUE ~ "Q1 (Most Disorganized)"
        )
      ) %>%
      group_by(quartile) %>%
      summarise(
        n_samples = n(),
        mean_coherence = mean(mean_coherence, na.rm = TRUE),
        samples = paste(sample_id, collapse = ", "),
        .groups = 'drop'
      )
    
    readr::write_csv(quartile_groups, file.path(output_dir, "quartile_based_organization_groups.csv"))
    if (verbose) cat("  ✓ quartile_based_organization_groups.csv\n")
    
    # Treatment comparison if available
    if ("treatment_status" %in% colnames(sample_rankings)) {
      treatment_comparison <- sample_rankings %>%
        filter(!is.na(treatment_status)) %>%
        group_by(treatment_status, organization_class) %>%
        summarise(n_samples = n(), .groups = 'drop') %>%
        tidyr::pivot_wider(names_from = organization_class, values_from = n_samples, values_fill = 0)
      
      readr::write_csv(treatment_comparison, file.path(output_dir, "treatment_organization_comparison.csv"))
      if (verbose) cat("  ✓ treatment_organization_comparison.csv\n")
    }
  }
  
  # ========================================================================
  # GENERATE VISUALIZATION
  # ========================================================================
  
  if (generate_plot) {
    if (verbose) cat("\nGenerating organization ranking plot...\n")
    
    # Prepare data for plotting
    plot_data <- sample_rankings %>%
      mutate(
        sample_id = factor(sample_id, levels = sample_id),  # Maintain ranking order
        organization_class = factor(organization_class, levels = c(
          "Highly Organized", "Moderately Organized", "Organized",
          "Intermediate", "Moderately Disorganized", "Highly Disorganized"
        ))
      )
    
    # Color palette for organization classes
    org_colors <- c(
      "Highly Organized" = "#004225",      # Dark green
      "Moderately Organized" = "#2E8B57",  # Sea green  
      "Organized" = "#66CDAA",             # Medium aquamarine
      "Intermediate" = "#FFD700",          # Gold
      "Moderately Disorganized" = "#FF6347", # Tomato
      "Highly Disorganized" = "#8B0000"   # Dark red
    )
    
    ranking_plot <- ggplot(plot_data, aes(x = rank_composite, y = mean_coherence, fill = organization_class)) +
      geom_col(alpha = 0.8, width = 0.8) +
      geom_hline(yintercept = coherence_threshold, linetype = "dashed", color = "black", size = 1) +
      geom_hline(yintercept = coherence_threshold + 0.05, linetype = "dotted", color = "darkgreen", size = 0.8) +
      geom_hline(yintercept = coherence_threshold - 0.05, linetype = "dotted", color = "darkred", size = 0.8) +
      scale_fill_manual(values = org_colors, name = "Organization Class") +
      theme_minimal() +
      theme(
        axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Sample Rank (1 = Most Organized)",
        y = "Mean Spatial Coherence Score",
        title = "Sample Organization Ranking",
        subtitle = paste0("Ranked by composite organization score | Threshold = ", coherence_threshold)
      ) +
      scale_x_continuous(breaks = seq(1, nrow(plot_data), by = max(1, floor(nrow(plot_data)/20)))) +
      scale_y_continuous(limits = c(0, max(plot_data$mean_coherence) * 1.05))
    
    if (save_csvs) {
      ggsave(file.path(output_dir, "sample_organization_ranking_plot.pdf"), 
             ranking_plot, width = max(12, nrow(plot_data) * 0.15), height = 8, dpi = 300)
      ggsave(file.path(output_dir, "sample_organization_ranking_plot.png"), 
             ranking_plot, width = max(12, nrow(plot_data) * 0.15), height = 8, dpi = 300)
      if (verbose) cat("  ✓ sample_organization_ranking_plot.pdf/.png\n")
    }
  }
  
  # ========================================================================
  # SUMMARY STATISTICS
  # ========================================================================
  
  if (verbose) {
    cat("\n============================================================================\n")
    cat("SAMPLE ORGANIZATION QUANTIFICATION SUMMARY\n")
    cat("============================================================================\n")
    
    cat("Total samples analyzed:", nrow(sample_rankings), "\n")
    cat("Coherence threshold used:", coherence_threshold, "\n\n")
    
    cat("Organization class distribution:\n")
    org_dist <- table(sample_rankings$organization_class)
    for (i in 1:length(org_dist)) {
      cat("  ", names(org_dist)[i], ":", org_dist[i], paste0("(", round(100*org_dist[i]/sum(org_dist), 1), "%)"), "\n")
    }
    
    cat("\nTop 5 most organized samples:\n")
    top_5 <- head(sample_rankings, 5)
    for (i in 1:nrow(top_5)) {
      cat("  ", i, ". ", top_5$sample_id[i], " (coherence = ", round(top_5$mean_coherence[i], 3), 
          ", class = ", top_5$organization_class[i], ")\n", sep = "")
    }
    
    cat("\nTop 5 most disorganized samples:\n")
    bottom_5 <- tail(sample_rankings, 5)
    for (i in nrow(bottom_5):1) {
      rank_pos <- nrow(sample_rankings) - i + 1
      cat("  ", rank_pos, ". ", bottom_5$sample_id[i], " (coherence = ", round(bottom_5$mean_coherence[i], 3), 
          ", class = ", bottom_5$organization_class[i], ")\n", sep = "")
    }
    
    cat("\nOverall statistics:\n")
    cat("  Mean coherence across all samples:", round(mean(sample_rankings$mean_coherence, na.rm = TRUE), 3), "\n")
    cat("  Median coherence:", round(median(sample_rankings$mean_coherence, na.rm = TRUE), 3), "\n")
    cat("  Standard deviation:", round(sd(sample_rankings$mean_coherence, na.rm = TRUE), 3), "\n")
    cat("  Range:", round(min(sample_rankings$mean_coherence, na.rm = TRUE), 3), "-", 
        round(max(sample_rankings$mean_coherence, na.rm = TRUE), 3), "\n")
    
    if (save_csvs) {
      cat("\n=== FILES CREATED ===\n")
      cat("Main analysis files:\n")
      cat("  • sample_organization_comprehensive_ranking.csv - Complete analysis\n")
      cat("  • sample_organization_simplified_ranking.csv - Key metrics only\n")
      cat("  • organization_class_summary.csv - Summary by organization class\n")
      cat("\nSpecialty files:\n") 
      cat("  • top_10_most_organized_samples.csv\n")
      cat("  • top_10_most_disorganized_samples.csv\n")
      cat("  • quartile_based_organization_groups.csv\n")
      if ("treatment_status" %in% colnames(sample_rankings)) {
        cat("  • treatment_organization_comparison.csv\n")
      }
      cat("  • sample_organization_ranking_plot.pdf/.png\n")
      
      cat(paste0("\nAll files saved in: ", output_dir, "/\n"))
    }
    
    cat("\n============================================================================\n")
  }
  
  return(list(
    rankings = sample_rankings,
    class_summary = if(exists("class_summary")) class_summary else NULL,
    treatment_summary = if(exists("treatment_summary")) treatment_summary else NULL,
    plot = if(generate_plot) ranking_plot else NULL,
    threshold_used = coherence_threshold
  ))
}

# ============================================================================
# CONVENIENCE FUNCTIONS
# ============================================================================

# Quick function to get just the top N most/least organized samples
get_top_organized_samples <- function(results, n = 10, least = FALSE) {
  analysis <- quantify_sample_organization(results, save_csvs = FALSE, generate_plot = FALSE, verbose = FALSE)
  
  if (least) {
    return(tail(analysis$rankings, n))
  } else {
    return(head(analysis$rankings, n))
  }
}

# Function to compare specific samples
compare_samples <- function(results, sample_ids) {
  analysis <- quantify_sample_organization(results, save_csvs = FALSE, generate_plot = FALSE, verbose = FALSE)
  
  comparison <- analysis$rankings %>%
    filter(sample_id %in% sample_ids) %>%
    select(sample_id, rank_composite, mean_coherence, organization_class, 
           prop_organized_ecotypes, composite_score) %>%
    arrange(rank_composite)
  
  return(comparison)
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# # Example usage (uncomment to run):
# 
# # Load your results first
# source("00_load_spatialcoherence_data.R")
# results <- load_spatialcoherence_data("your_data_directory")
# 
# # Run comprehensive quantification
# organization_analysis <- quantify_sample_organization(
#   results = results,
#   output_dir = "sample_organization_analysis",
#   coherence_threshold = 0.47,  # or NULL to use default from results
#   save_csvs = TRUE,
#   generate_plot = TRUE,
#   verbose = TRUE
# )
# 
# # Quick access to top organized samples
# top_10_organized <- get_top_organized_samples(results, n = 10, least = FALSE)
# print(top_10_organized[, c("sample_id", "mean_coherence", "organization_class")])
# 
# # Quick access to most disorganized samples  
# top_10_disorganized <- get_top_organized_samples(results, n = 10, least = TRUE)
# print(top_10_disorganized[, c("sample_id", "mean_coherence", "organization_class")])
# 
# # Compare specific samples
# sample_comparison <- compare_samples(results, c("sample1", "sample2", "sample3"))
# print(sample_comparison)

cat("Sample Organization Quantifier loaded successfully!\n")
cat("Main function: quantify_sample_organization()\n")
cat("Quick functions: get_top_organized_samples(), compare_samples()\n")
cat("Use verbose=TRUE to see detailed progress and results.\n")