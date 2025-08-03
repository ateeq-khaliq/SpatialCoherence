#!/usr/bin/env Rscript

# ============================================================================
# MASTER SpatialCoherence Plotting Script - Complete Analysis Suite
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 2.0.0
# 
# This master script incorporates all plotting functionality and adds missing
# plots identified from your analysis. It's designed to be non-hardcoded and
# flexible for different datasets.

# ============================================================================
# CONFIGURATION AND SETUP
# ============================================================================

# Set up paths (modify these as needed)
SCRIPT_DIR <- "/Users/akhaliq/Downloads/SpatialCoherence-main/plotting"
DATA_DIR <- "/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test"
OUTPUT_DIR <- "comprehensive_spatial_plots"

# Required libraries
required_packages <- c("ggplot2", "dplyr", "reshape2", "patchwork", "viridis", 
                      "RColorBrewer", "scales", "broom", "tidyr")

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
# LOAD EXISTING PLOTTING FUNCTIONS
# ============================================================================

cat("=== LOADING SPATIAL COHERENCE PLOTTING SUITE ===\n")

# Load data loader
if (file.exists(file.path(SCRIPT_DIR, "00_load_spatialcoherence_data.R"))) {
  source(file.path(SCRIPT_DIR, "00_load_spatialcoherence_data.R"))
  cat("✓ Data loader loaded\n")
} else {
  stop("Cannot find data loader script")
}

# Load plotting modules
plotting_scripts <- c(
  "01_basic_coherence_plots.R",
  "02_treatment_analysis_plots.R", 
  "03_segmentation_analysis_plots.R",
  "04_heatmap_visualization_plots.R"
)

for (script in plotting_scripts) {
  script_path <- file.path(SCRIPT_DIR, script)
  if (file.exists(script_path)) {
    source(script_path)
    cat("✓", script, "loaded\n")
  } else {
    cat("⚠️  Warning:", script, "not found\n")
  }
}

# ============================================================================
# ENHANCED PLOTTING CONFIGURATION
# ============================================================================

# Extended plot configuration
ENHANCED_PLOT_CONFIG <- list(
  # Dimensions
  width = 12,
  height = 8,
  dpi = 300,
  
  # Colors
  organized_color = "#2E8B57",      # Sea green
  disorganized_color = "#DC143C",   # Crimson
  intermediate_color = "#FFD700",   # Gold
  
  # Treatment colors
  treated_color = "#CF1C90",        # Magenta
  untreated_color = "#11A579",      # Teal
  
  # Functional colors
  functional_colors = c(
    "Hypoxic Core" = "#8B0000",
    "Proliferative" = "#4169E1", 
    "Stromal" = "#228B22",
    "Immune" = "#FF6347",
    "Vascular" = "#9370DB",
    "Normal" = "#DAA520",
    "Other" = "#808080"
  ),
  
  # Text sizes
  title_size = 16,
  subtitle_size = 12,
  axis_title_size = 14,
  axis_text_size = 12,
  legend_text_size = 10,
  
  # Output formats
  formats = c("pdf", "png")
)

# ============================================================================
# MISSING PLOT IMPLEMENTATIONS
# ============================================================================

# Plot 05: Cell Density by Treatment, Organization, and Functional Type (3-panel)
create_05_comprehensive_cell_density <- function(results, 
                                                output_dir = ".",
                                                treatment_column = "treatment",
                                                save_plots = TRUE) {
  
  cat("Creating Plot 05: Comprehensive Cell Density Analysis (3-panel)...\n")
  
  if (is.null(results$detailed_results)) {
    cat("  ! No detailed results available. Skipping comprehensive density plot.\n")
    return(NULL)
  }
  
  # Prepare data for all three panels
  detailed_data <- results$detailed_results
  
  # Panel A: Treatment comparison
  if (treatment_column %in% colnames(detailed_data)) {
    panel_a_data <- detailed_data %>%
      group_by(sample, !!sym(treatment_column)) %>%
      summarise(cell_count_per_spot = n(), .groups = 'drop') %>%
      rename(treatment_status = !!sym(treatment_column)) %>%
      filter(!is.na(treatment_status))
    
    # Statistical test
    if (length(unique(panel_a_data$treatment_status)) == 2) {
      test_a <- t.test(cell_count_per_spot ~ treatment_status, data = panel_a_data)
      p_val_a <- format.pval(test_a$p.value, digits = 3)
    } else {
      p_val_a <- "N/A"
    }
    
    panel_a <- ggplot(panel_a_data, aes(x = treatment_status, y = cell_count_per_spot)) +
      geom_boxplot(aes(fill = treatment_status), alpha = 0.7, outlier.size = 0.8) +
      geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
      scale_fill_manual(values = c("Untreated" = ENHANCED_PLOT_CONFIG$untreated_color, 
                                  "Treated" = ENHANCED_PLOT_CONFIG$treated_color)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
        axis.title.y = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
        axis.text = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size),
        legend.position = "none",
        plot.title = element_text(size = ENHANCED_PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = ENHANCED_PLOT_CONFIG$subtitle_size, hjust = 0.5)
      ) +
      labs(
        x = "Treatment Status",
        y = "cell count per spot",
        title = "Cell Density by Treatment Status",
        subtitle = paste("p =", p_val_a)
      )
  } else {
    panel_a <- ggplot() + 
      annotate("text", x = 0.5, y = 0.5, label = "No treatment data available", size = 6) +
      theme_void()
  }
  
  # Panel B: Organization zones
  # Calculate sample-level coherence and organization
  coherence_threshold <- ifelse(!is.null(results$analysis_parameters$coherence_threshold),
                               results$analysis_parameters$coherence_threshold, 0.47)
  
  sample_coherence <- detailed_data %>%
    group_by(sample) %>%
    summarise(
      mean_coherence = mean(coherence_score, na.rm = TRUE),
      cell_count_per_spot = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      organization_zone = case_when(
        mean_coherence >= coherence_threshold + 0.08 ~ "structured",
        mean_coherence >= coherence_threshold - 0.08 ~ "intermediate",
        TRUE ~ "disorganized"
      )
    )
  
  panel_b <- ggplot(sample_coherence, aes(x = organization_zone, y = cell_count_per_spot)) +
    geom_boxplot(aes(fill = organization_zone), alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = c("structured" = ENHANCED_PLOT_CONFIG$organized_color,
                                "intermediate" = ENHANCED_PLOT_CONFIG$intermediate_color,
                                "disorganized" = ENHANCED_PLOT_CONFIG$disorganized_color)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size),
      legend.position = "none",
      plot.title = element_text(size = ENHANCED_PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = ENHANCED_PLOT_CONFIG$subtitle_size, hjust = 0.5)
    ) +
    labs(
      x = "Organization Zone",
      y = "cell count per spot",
      title = "Cell Density by Organization Zone",
      subtitle = paste("Threshold =", coherence_threshold)
    )
  
  # Panel C: Functional classification
  if (!is.null(results$ecotype_annotations)) {
    detailed_data$ecotype_function <- results$ecotype_annotations[detailed_data$metaprogram]
  } else {
    # Infer functional types from ecotype names
    detailed_data$ecotype_function <- sapply(detailed_data$metaprogram, function(x) {
      x_lower <- tolower(as.character(x))
      if (grepl("10", x)) return("Hypoxic Core")
      if (grepl("07|03", x)) return("Proliferative")
      if (grepl("06|08|09", x)) return("Stromal")
      if (grepl("01|02", x)) return("Immune")
      if (grepl("05", x)) return("Vascular")
      if (grepl("04", x)) return("Normal")
      return("Other")
    })
  }
  
  functional_data <- detailed_data %>%
    group_by(sample, ecotype_function) %>%
    summarise(cell_count_per_spot = n(), .groups = 'drop') %>%
    filter(!is.na(ecotype_function) & ecotype_function != "Other")
  
  panel_c <- ggplot(functional_data, aes(x = ecotype_function, y = cell_count_per_spot)) +
    geom_boxplot(aes(fill = ecotype_function), alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = ENHANCED_PLOT_CONFIG$functional_colors) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text.x = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size),
      legend.position = "none",
      plot.title = element_text(size = ENHANCED_PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = ENHANCED_PLOT_CONFIG$subtitle_size, hjust = 0.5)
    ) +
    labs(
      x = "Functional Ecotype",
      y = "cell count per spot",
      title = "Cell Density by Functional Classification",
      subtitle = "Manuscript-based ecotypes"
    )
  
  # Combine panels
  combined_plot <- panel_a | panel_b | panel_c
  
  if (save_plots) {
    output_subdir <- file.path(output_dir, "05_comprehensive_plots")
    if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)
    
    ggsave(file.path(output_subdir, "05_comprehensive_cell_density_analysis.pdf"),
           combined_plot, width = 20, height = 6, dpi = ENHANCED_PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "05_comprehensive_cell_density_analysis.png"),
           combined_plot, width = 20, height = 6, dpi = ENHANCED_PLOT_CONFIG$dpi)
    
    cat("  ✓ Saved: 05_comprehensive_cell_density_analysis\n")
  }
  
  return(combined_plot)
}

# Plot 06: Sample-level Organization Summary with Statistics
create_06_sample_organization_summary <- function(results,
                                                 output_dir = ".",
                                                 save_plots = TRUE) {
  
  cat("Creating Plot 06: Sample Organization Summary with Statistics...\n")
  
  if (is.null(results$detailed_results)) {
    cat("  ! No detailed results available. Skipping sample organization summary.\n")
    return(NULL)
  }
  
  # Calculate sample-level statistics
  coherence_threshold <- ifelse(!is.null(results$analysis_parameters$coherence_threshold),
                               results$analysis_parameters$coherence_threshold, 0.47)
  
  sample_stats <- results$detailed_results %>%
    group_by(sample) %>%
    summarise(
      mean_coherence = mean(coherence_score, na.rm = TRUE),
      sd_coherence = sd(coherence_score, na.rm = TRUE),
      n_ecotypes = n_distinct(metaprogram),
      n_observations = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      organization_status = case_when(
        mean_coherence >= coherence_threshold + 0.08 ~ "Organized",
        mean_coherence >= coherence_threshold - 0.08 ~ "Intermediate",
        TRUE ~ "Disorganized"
      ),
      se_coherence = sd_coherence / sqrt(n_observations)
    ) %>%
    arrange(desc(mean_coherence))
  
  # Create summary statistics table
  org_summary <- sample_stats %>%
    group_by(organization_status) %>%
    summarise(
      n_samples = n(),
      mean_coherence_group = mean(mean_coherence),
      sd_coherence_group = sd(mean_coherence),
      .groups = 'drop'
    )
  
  # Plot with error bars and organization coloring
  sample_stats$sample <- factor(sample_stats$sample, levels = sample_stats$sample)
  
  p6 <- ggplot(sample_stats, aes(x = sample, y = mean_coherence, fill = organization_status)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_errorbar(aes(ymin = pmax(mean_coherence - se_coherence, 0),
                     ymax = pmin(mean_coherence + se_coherence, 1)),
                 width = 0.4, size = 0.5) +
    geom_hline(yintercept = coherence_threshold, linetype = "dashed", color = "black", size = 0.8) +
    geom_hline(yintercept = coherence_threshold + 0.08, linetype = "dotted", color = "darkgreen", size = 0.6) +
    geom_hline(yintercept = coherence_threshold - 0.08, linetype = "dotted", color = "darkred", size = 0.6) +
    scale_fill_manual(
      values = c("Organized" = ENHANCED_PLOT_CONFIG$organized_color,
                "Intermediate" = ENHANCED_PLOT_CONFIG$intermediate_color,
                "Disorganized" = ENHANCED_PLOT_CONFIG$disorganized_color),
      name = "Organization\nStatus"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size),
      axis.title = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = ENHANCED_PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = ENHANCED_PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = ENHANCED_PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = ENHANCED_PLOT_CONFIG$legend_text_size),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "Mean Coherence Score",
      title = "Sample Organization Summary with Statistics",
      subtitle = paste0("Threshold = ", coherence_threshold, " | Total Samples: ", nrow(sample_stats))
    )
  
  if (save_plots) {
    output_subdir <- file.path(output_dir, "06_sample_organization")
    if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)
    
    # Save plot
    ggsave(file.path(output_subdir, "06_sample_organization_summary_with_stats.pdf"),
           p6, width = max(16, nrow(sample_stats) * 0.3), height = 8,
           dpi = ENHANCED_PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "06_sample_organization_summary_with_stats.png"),
           p6, width = max(16, nrow(sample_stats) * 0.3), height = 8,
           dpi = ENHANCED_PLOT_CONFIG$dpi)
    
    # Save summary tables
    write.csv(sample_stats, file.path(output_subdir, "sample_statistics_detailed.csv"), row.names = FALSE)
    write.csv(org_summary, file.path(output_subdir, "organization_summary_statistics.csv"), row.names = FALSE)
    
    cat("  ✓ Saved: 06_sample_organization_summary_with_stats\n")
    cat("  ✓ Saved: sample_statistics_detailed.csv\n")
    cat("  ✓ Saved: organization_summary_statistics.csv\n")
  }
  
  return(list(plot = p6, sample_stats = sample_stats, org_summary = org_summary))
}

# Plot 07: Ecotype Ranking and Classification Overview
create_07_ecotype_ranking_overview <- function(results,
                                              output_dir = ".",
                                              save_plots = TRUE) {
  
  cat("Creating Plot 07: Ecotype Ranking and Classification Overview...\n")
  
  if (is.null(results$mean_coherence)) {
    cat("  ! No mean coherence data available. Skipping ecotype ranking.\n")
    return(NULL)
  }
  
  # Prepare ranking data
  ranking_data <- data.frame(
    ecotype = names(results$mean_coherence),
    mean_coherence = as.numeric(results$mean_coherence),
    organization = results$organization_results[names(results$mean_coherence)],
    stringsAsFactors = FALSE
  )
  
  # Add functional annotations if available
  if (!is.null(results$ecotype_annotations)) {
    ranking_data$functional_type <- results$ecotype_annotations[ranking_data$ecotype]
  } else {
    # Infer functional types
    ranking_data$functional_type <- sapply(ranking_data$ecotype, function(x) {
      x_str <- as.character(x)
      if (grepl("10", x_str)) return("Hypoxic Core")
      if (grepl("07|03", x_str)) return("Proliferative")
      if (grepl("06|08|09", x_str)) return("Stromal")
      if (grepl("01|02", x_str)) return("Immune")
      if (grepl("05", x_str)) return("Vascular")
      if (grepl("04", x_str)) return("Normal")
      return("Other")
    })
  }
  
  # Order by coherence
  ranking_data <- ranking_data[order(ranking_data$mean_coherence, decreasing = TRUE), ]
  ranking_data$rank <- 1:nrow(ranking_data)
  ranking_data$ecotype <- factor(ranking_data$ecotype, levels = ranking_data$ecotype)
  
  # Create ranking plot
  p7 <- ggplot(ranking_data, aes(x = rank, y = mean_coherence)) +
    geom_segment(aes(x = rank, xend = rank, y = 0, yend = mean_coherence, color = functional_type),
                size = 3, alpha = 0.8) +
    geom_point(aes(color = functional_type, shape = organization), size = 4) +
    geom_text(aes(label = ecotype), vjust = -0.5, hjust = 0.5, size = 3, fontface = "bold") +
    scale_color_manual(values = ENHANCED_PLOT_CONFIG$functional_colors, name = "Functional Type") +
    scale_shape_manual(values = c("Organized" = 16, "Disorganized" = 17, "Intermediate" = 15),
                      name = "Organization") +
    theme_minimal() +
    theme(
      axis.title = element_text(size = ENHANCED_PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = ENHANCED_PLOT_CONFIG$axis_text_size),
      plot.title = element_text(size = ENHANCED_PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = ENHANCED_PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = ENHANCED_PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = ENHANCED_PLOT_CONFIG$legend_text_size),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Rank (by Mean Coherence)",
      y = "Mean Spatial Coherence Score",
      title = "Ecotype Ranking and Classification Overview",
      subtitle = "Functional types and organization status by spatial coherence"
    ) +
    scale_x_continuous(breaks = ranking_data$rank, labels = ranking_data$rank) +
    scale_y_continuous(limits = c(0, max(ranking_data$mean_coherence) * 1.1))
  
  if (save_plots) {
    output_subdir <- file.path(output_dir, "07_ecotype_ranking")
    if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)
    
    ggsave(file.path(output_subdir, "07_ecotype_ranking_and_classification.pdf"),
           p7, width = 12, height = 8, dpi = ENHANCED_PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "07_ecotype_ranking_and_classification.png"),
           p7, width = 12, height = 8, dpi = ENHANCED_PLOT_CONFIG$dpi)
    
    # Save ranking table
    write.csv(ranking_data, file.path(output_subdir, "ecotype_ranking_table.csv"), row.names = FALSE)
    
    cat("  ✓ Saved: 07_ecotype_ranking_and_classification\n")
    cat("  ✓ Saved: ecotype_ranking_table.csv\n")
  }
  
  return(list(plot = p7, ranking_data = ranking_data))
}

# ============================================================================
# MASTER PLOTTING FUNCTION
# ============================================================================

create_comprehensive_spatial_analysis <- function(data_dir,
                                                  output_dir = "comprehensive_spatial_plots",
                                                  treatment_column = "treatment",
                                                  save_plots = TRUE,
                                                  generate_missing_plots = TRUE) {
  
  cat("============================================================================\n")
  cat("COMPREHENSIVE SPATIAL COHERENCE ANALYSIS - MASTER SUITE\n") 
  cat("============================================================================\n\n")
  
  # Step 1: Load data
  cat("Step 1: Loading spatial coherence data...\n")
  results <- load_spatialcoherence_data(data_dir, verbose = TRUE)
  
  if (is.null(results) || length(results) == 0) {
    stop("Failed to load data from: ", data_dir)
  }
  
  # Step 2: Create output directory structure
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  all_plots <- list()
  
  # Step 3: Generate existing plots (01-04)
  cat("\nStep 3: Generating existing plot suites...\n")
  
  # Basic plots (01)
  if (exists("create_all_basic_plots")) {
    tryCatch({
      basic_plots <- create_all_basic_plots(results, output_dir, save_plots)
      all_plots <- c(all_plots, basic_plots)
      cat("✓ Basic plots generated\n")
    }, error = function(e) {
      cat("⚠️  Error in basic plots:", e$message, "\n")
    })
  }
  
  # Treatment plots (02)
  if (exists("create_all_treatment_plots")) {
    tryCatch({
      treatment_plots <- create_all_treatment_plots(results, output_dir, save_plots, treatment_column)
      all_plots <- c(all_plots, treatment_plots)
      cat("✓ Treatment plots generated\n")
    }, error = function(e) {
      cat("⚠️  Error in treatment plots:", e$message, "\n")
    })
  }
  
  # Segmentation plots (03)
  if (exists("create_all_segmentation_plots")) {
    tryCatch({
      segmentation_plots <- create_all_segmentation_plots(results, output_dir, save_plots, treatment_column)
      all_plots <- c(all_plots, segmentation_plots)
      cat("✓ Segmentation plots generated\n")
    }, error = function(e) {
      cat("⚠️  Error in segmentation plots:", e$message, "\n")
    })
  }
  
  # Heatmap plots (04)
  if (exists("create_all_heatmap_plots")) {
    tryCatch({
      heatmap_plots <- create_all_heatmap_plots(results, output_dir, save_plots)
      all_plots <- c(all_plots, heatmap_plots)
      cat("✓ Heatmap plots generated\n")
    }, error = function(e) {
      cat("⚠️  Error in heatmap plots:", e$message, "\n")
    })
  }
  
  # Step 4: Generate missing plots (05+)
  if (generate_missing_plots) {
    cat("\nStep 4: Generating missing/additional plots...\n")
    
    # Plot 05: Comprehensive cell density (matches your image)
    tryCatch({
      plot_05 <- create_05_comprehensive_cell_density(results, output_dir, treatment_column, save_plots)
      if (!is.null(plot_05)) {
        all_plots$plot_05 <- plot_05
        cat("✓ Plot 05: Comprehensive cell density generated\n")
      }
    }, error = function(e) {
      cat("⚠️  Error in plot 05:", e$message, "\n")
    })
    
    # Plot 06: Sample organization summary with statistics
    tryCatch({
      plot_06_result <- create_06_sample_organization_summary(results, output_dir, save_plots)
      if (!is.null(plot_06_result)) {
        all_plots$plot_06 <- plot_06_result$plot
        cat("✓ Plot 06: Sample organization summary generated\n")
      }
    }, error = function(e) {
      cat("⚠️  Error in plot 06:", e$message, "\n")
    })
    
    # Plot 07: Ecotype ranking overview
    tryCatch({
      plot_07_result <- create_07_ecotype_ranking_overview(results, output_dir, save_plots)
      if (!is.null(plot_07_result)) {
        all_plots$plot_07 <- plot_07_result$plot
        cat("✓ Plot 07: Ecotype ranking overview generated\n")
      }
    }, error = function(e) {
      cat("⚠️  Error in plot 07:", e$message, "\n")
    })
  }
  
  # Step 5: Generate summary report
  cat("\nStep 5: Generating analysis summary...\n")
  
  summary_report <- list(
    data_summary = results$metadata,
    total_plots_generated = length(all_plots),
    plot_categories = c(
      "Basic coherence plots" = sum(grepl("plot_0[1-4]", names(all_plots))),
      "Treatment analysis plots" = sum(grepl("plot_0[5-8]", names(all_plots))),
      "Segmentation plots" = sum(grepl("plot_(09|10|11|12)", names(all_plots))),
      "Heatmap plots" = sum(grepl("plot_1[3-6]", names(all_plots))),
      "Additional plots" = sum(grepl("plot_0[5-7]", names(all_plots)))
    ),
    output_directory = output_dir,
    analysis_date = Sys.time()
  )
  
  # Save summary report
  if (save_plots) {
    saveRDS(summary_report, file.path(output_dir, "analysis_summary_report.rds"))
    
    # Create a text summary
    summary_text <- c(
      "=== COMPREHENSIVE SPATIAL COHERENCE ANALYSIS SUMMARY ===",
      paste("Analysis Date:", Sys.time()),
      paste("Data Directory:", data_dir),
      paste("Output Directory:", output_dir),
      "",
      "=== DATA SUMMARY ===",
      paste("Total Ecotypes:", summary_report$data_summary$n_ecotypes),
      paste("Total Samples:", ifelse(is.na(summary_report$data_summary$n_samples), "N/A", summary_report$data_summary$n_samples)),
      paste("Has Treatment Data:", summary_report$data_summary$has_treatment_data),
      paste("Has Compositional Data:", summary_report$data_summary$has_compositional_data),
      paste("Has Functional Annotations:", summary_report$data_summary$has_annotations),
      "",
      "=== PLOTS GENERATED ===",
      paste("Total Plots:", summary_report$total_plots_generated),
      paste("Basic Plots:", summary_report$plot_categories["Basic coherence plots"]),
      paste("Treatment Plots:", summary_report$plot_categories["Treatment analysis plots"]),
      paste("Segmentation Plots:", summary_report$plot_categories["Segmentation plots"]),
      paste("Heatmap Plots:", summary_report$plot_categories["Heatmap plots"]),
      paste("Additional Plots:", summary_report$plot_categories["Additional plots"]),
      "",
      "=== FILES STRUCTURE ===",
      "01_basic_plots/ - Core coherence analysis",
      "02_treatment_plots/ - Treatment effect analysis",
      "03_segmentation_plots/ - Cell density and organization",
      "04_heatmap_plots/ - Complex visualizations",
      "05_comprehensive_plots/ - Multi-panel summaries",
      "06_sample_organization/ - Sample-level statistics",
      "07_ecotype_ranking/ - Ecotype classification overview",
      "",
      "=== RECOMMENDED NEXT STEPS ===",
      "1. Review individual plots in each subdirectory",
      "2. Check CSV files for detailed statistics",
      "3. Use heatmaps for publication-quality figures",
      "4. Examine treatment effects if applicable",
      "5. Validate functional classifications"
    )
    
    writeLines(summary_text, file.path(output_dir, "ANALYSIS_SUMMARY.txt"))
  }
  
  # Step 6: Final output
  cat("\n============================================================================\n")
  cat("COMPREHENSIVE SPATIAL COHERENCE ANALYSIS COMPLETE\n")
  cat("============================================================================\n")
  cat("Total plots generated:", summary_report$total_plots_generated, "\n")
  cat("Output directory:", output_dir, "\n")
  
  if (summary_report$total_plots_generated > 0) {
    cat("\n📊 PLOT BREAKDOWN:\n")
    for (category in names(summary_report$plot_categories)) {
      if (summary_report$plot_categories[category] > 0) {
        cat("  •", category, ":", summary_report$plot_categories[category], "\n")
      }
    }
    
    cat("\n📁 KEY FILES CREATED:\n")
    if (save_plots) {
      cat("  • ANALYSIS_SUMMARY.txt - Overview of all analyses\n")
      cat("  • analysis_summary_report.rds - Detailed R object\n")
      cat("  • Individual plot directories with PDF and PNG versions\n")
      cat("  • CSV files with statistical summaries\n")
    }
    
    cat("\n✅ Analysis completed successfully!\n")
    cat("🔍 Check the ANALYSIS_SUMMARY.txt file for detailed information.\n")
  } else {
    cat("\n⚠️  No plots were generated. Please check your data and configuration.\n")
  }
  
  return(list(
    plots = all_plots,
    summary = summary_report,
    results = results
  ))
}

# ============================================================================
# ADDITIONAL UTILITY FUNCTIONS
# ============================================================================

# Function to validate data requirements
validate_data_requirements <- function(results) {
  cat("Validating data requirements...\n")
  
  requirements <- list(
    basic_plots = !is.null(results$mean_coherence),
    treatment_plots = !is.null(results$detailed_results) && 
                     any(c("treatment", "nac_treatment") %in% colnames(results$detailed_results)),
    segmentation_plots = !is.null(results$detailed_results),
    heatmap_plots = !is.null(results$coherence_matrix) || !is.null(results$detailed_results),
    comprehensive_plots = !is.null(results$detailed_results)
  )
  
  cat("Data validation results:\n")
  for (req in names(requirements)) {
    status <- ifelse(requirements[[req]], "✓", "✗")
    cat("  ", status, req, "\n")
  }
  
  return(requirements)
}

# Function to create a quick diagnostic plot
create_diagnostic_plot <- function(results, output_dir = ".") {
  cat("Creating diagnostic overview plot...\n")
  
  if (is.null(results$mean_coherence)) {
    return(NULL)
  }
  
  # Simple diagnostic showing data completeness
  completeness_data <- data.frame(
    component = c("Mean Coherence", "Organization Results", "Detailed Results", 
                 "Coherence Matrix", "Treatment Data", "Compositional Data", 
                 "Functional Annotations"),
    available = c(
      !is.null(results$mean_coherence),
      !is.null(results$organization_results),
      !is.null(results$detailed_results),
      !is.null(results$coherence_matrix),
      !is.null(results$treatment_results) || 
        (!is.null(results$detailed_results) && "treatment" %in% colnames(results$detailed_results)),
      !is.null(results$compositional_results),
      !is.null(results$ecotype_annotations)
    ),
    stringsAsFactors = FALSE
  )
  
  completeness_data$status <- ifelse(completeness_data$available, "Available", "Missing")
  
  diagnostic_plot <- ggplot(completeness_data, aes(x = component, y = 1, fill = status)) +
    geom_tile(color = "white", size = 1) +
    geom_text(aes(label = status), fontface = "bold", size = 4) +
    scale_fill_manual(values = c("Available" = "#2E8B57", "Missing" = "#DC143C"),
                     name = "Data Status") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_blank(),
      axis.title = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      panel.grid = element_blank()
    ) +
    labs(title = "Data Completeness Diagnostic") +
    coord_fixed()
  
  ggsave(file.path(output_dir, "data_diagnostic.pdf"), diagnostic_plot, 
         width = 12, height = 3, dpi = 300)
  
  return(diagnostic_plot)
}

# Function to generate a methods summary
generate_methods_summary <- function(results, output_dir = ".") {
  
  methods_text <- c(
    "=== SPATIAL COHERENCE ANALYSIS METHODS ===",
    "",
    "## Data Processing",
    paste("- Coherence threshold:", ifelse(!is.null(results$analysis_parameters$coherence_threshold), 
                                         results$analysis_parameters$coherence_threshold, "0.47 (default)")),
    paste("- Enhanced analysis:", ifelse(!is.null(results$analysis_parameters$enhanced_analysis), 
                                       results$analysis_parameters$enhanced_analysis, "FALSE (default)")),
    paste("- Minimum spots per ecotype:", ifelse(!is.null(results$analysis_parameters$min_spots_per_ecotype), 
                                               results$analysis_parameters$min_spots_per_ecotype, "5 (default)")),
    "",
    "## Organization Classification",
    "Samples are classified into organization zones based on mean spatial coherence:",
    "- Structured: Mean coherence ≥ threshold + 0.08",
    "- Intermediate: Threshold - 0.08 ≤ Mean coherence < threshold + 0.08", 
    "- Disorganized: Mean coherence < threshold - 0.08",
    "",
    "## Statistical Tests",
    "- Two-sample t-tests for treatment comparisons",
    "- Chi-square or Fisher's exact tests for categorical associations",
    "- Effect sizes calculated using Cohen's d",
    "",
    "## Functional Classifications",
    if (!is.null(results$ecotype_annotations)) {
      "- Based on provided functional annotations"
    } else {
      "- Inferred from ecotype naming patterns"
    },
    "",
    "## Plot Descriptions",
    "- Basic plots: Core coherence analysis and organization",
    "- Treatment plots: Effects of therapeutic interventions",
    "- Segmentation plots: Cell density analysis by various factors",
    "- Heatmap plots: Complex visualizations with clustering",
    "- Comprehensive plots: Multi-panel summary analyses",
    "",
    "## Software",
    paste("- R version:", R.version.string),
    "- Key packages: ggplot2, dplyr, ComplexHeatmap, patchwork",
    paste("- Analysis date:", Sys.Date())
  )
  
  writeLines(methods_text, file.path(output_dir, "METHODS_SUMMARY.txt"))
  cat("✓ Methods summary saved to METHODS_SUMMARY.txt\n")
}

# ============================================================================
# EXAMPLE USAGE AND MAIN EXECUTION
# ============================================================================

# Main execution function
main_spatial_analysis <- function(data_dir = NULL, 
                                 output_dir = "comprehensive_spatial_plots",
                                 treatment_column = "treatment",
                                 save_plots = TRUE) {
  
  # Use provided paths or defaults
  if (is.null(data_dir)) {
    data_dir <- DATA_DIR
  }
  
  cat("Starting comprehensive spatial coherence analysis...\n")
  cat("Data directory:", data_dir, "\n")
  cat("Output directory:", output_dir, "\n")
  
  # Validate inputs
  if (!dir.exists(data_dir)) {
    stop("Data directory does not exist: ", data_dir)
  }
  
  # Run comprehensive analysis
  analysis_result <- create_comprehensive_spatial_analysis(
    data_dir = data_dir,
    output_dir = output_dir,
    treatment_column = treatment_column,
    save_plots = save_plots,
    generate_missing_plots = TRUE
  )
  
  # Generate additional outputs
  if (save_plots) {
    validate_data_requirements(analysis_result$results)
    create_diagnostic_plot(analysis_result$results, output_dir)
    generate_methods_summary(analysis_result$results, output_dir)
  }
  
  return(analysis_result)
}

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

# Uncomment the following lines to run the analysis automatically
# (Modify paths as needed)

# # Set your paths
# DATA_DIR <- "/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test"
# OUTPUT_DIR <- "comprehensive_spatial_plots"
# 
# # Run the complete analysis
# results <- main_spatial_analysis(
#   data_dir = DATA_DIR,
#   output_dir = OUTPUT_DIR,
#   treatment_column = "treatment",  # or "nac_treatment" depending on your data
#   save_plots = TRUE
# )
# 
# # Display summary
# cat("Analysis complete! Check", OUTPUT_DIR, "for results.\n")

# ============================================================================
# ALTERNATIVE: STEP-BY-STEP EXECUTION (Current approach)
# ============================================================================

cat("=== SPATIAL COHERENCE MASTER PLOTTING SCRIPT LOADED ===\n")
cat("Available functions:\n")
cat("  • main_spatial_analysis() - Complete automated analysis\n")
cat("  • create_comprehensive_spatial_analysis() - Core analysis function\n") 
cat("  • create_05_comprehensive_cell_density() - Missing plot 05\n")
cat("  • create_06_sample_organization_summary() - Missing plot 06\n")
cat("  • create_07_ecotype_ranking_overview() - Missing plot 07\n")
cat("  • validate_data_requirements() - Check data completeness\n")
cat("  • create_diagnostic_plot() - Generate diagnostic overview\n")
cat("\n")
cat("To run your current workflow with added missing plots:\n")
cat("  results <- main_spatial_analysis(\n")
cat("    data_dir = '", DATA_DIR, "',\n")
cat("    output_dir = 'comprehensive_plots',\n")
cat("    treatment_column = 'treatment'\n")
cat("  )\n")
cat("\n")
cat("Or continue with your step-by-step approach:\n")
cat("  # Your existing code will work, plus you now have plots 05-07 available\n")

# ============================================================================
# END OF SCRIPT
# ============================================================================