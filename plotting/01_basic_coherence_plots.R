#!/usr/bin/env Rscript

# ============================================================================
# SpatialCoherence Plotting Suite - Part 1: Basic Coherence Plots
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script generates basic spatial coherence visualization plots
# Part 1 of 4: Core coherence analysis plots

# Required libraries
required_packages <- c("ggplot2", "dplyr", "reshape2", "RColorBrewer", "viridis")

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
# PLOTTING CONFIGURATION
# ============================================================================

# Standard plot settings
PLOT_CONFIG <- list(
  # Dimensions
  width = 12,
  height = 8,
  dpi = 300,
  
  # Colors
  organized_color = "#2E8B57",      # Sea green
  disorganized_color = "#DC143C",   # Crimson
  intermediate_color = "#FFD700",   # Gold
  
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
# UTILITY FUNCTIONS
# ============================================================================

# Generate color palette for any number of ecotypes
generate_ecotype_colors <- function(ecotype_names) {
  n_ecotypes <- length(ecotype_names)
  
  if (n_ecotypes <= 10) {
    # Use predefined palette for small numbers
    base_colors <- c("#F2B701", "#F6CF71", "#66C5CC", "#E68310", "#f97b72",
                    "#7F3C8D", "#3969AC", "#EF4868", "#4b4b8f", "#B95FBB")
    colors <- base_colors[1:n_ecotypes]
  } else if (n_ecotypes <= 20) {
    # Use RColorBrewer for medium numbers
    colors <- RColorBrewer::brewer.pal(min(n_ecotypes, 12), "Set3")
    if (n_ecotypes > 12) {
      additional_colors <- RColorBrewer::brewer.pal(n_ecotypes - 12, "Pastel1")
      colors <- c(colors, additional_colors)
    }
  } else {
    # Use viridis for large numbers
    colors <- viridis::viridis(n_ecotypes, option = "D")
  }
  
  names(colors) <- ecotype_names
  return(colors)
}

# Create output directory structure with user-defined path
create_output_structure <- function(base_dir) {
  # Ensure base directory is created first
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  dirs <- c(
    file.path(base_dir, "01_basic_plots"),
    file.path(base_dir, "01_basic_plots", "pdf"),
    file.path(base_dir, "01_basic_plots", "png")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(dirs[1])
}

# Save plot in multiple formats
save_plot_formats <- function(plot, filename, output_dir, formats = PLOT_CONFIG$formats) {
  for (format in formats) {
    filepath <- file.path(output_dir, format, paste0(filename, ".", format))
    
    if (format == "pdf") {
      ggsave(filepath, plot, width = PLOT_CONFIG$width, height = PLOT_CONFIG$height, 
             device = "pdf", dpi = PLOT_CONFIG$dpi)
    } else if (format == "png") {
      ggsave(filepath, plot, width = PLOT_CONFIG$width, height = PLOT_CONFIG$height, 
             device = "png", dpi = PLOT_CONFIG$dpi)
    }
  }
  
  cat("  ✓ Saved:", filename, "in", paste(formats, collapse = ", "), "format(s)\n")
}

# ============================================================================
# PLOT 01: MEAN COHERENCE BY ECOTYPE (LINE PLOT STYLE - MATCHING YOUR IMAGE 1)
# ============================================================================

create_01_mean_coherence_lineplot <- function(mean_coherence_data, 
                                            organization_results = NULL,
                                            output_dir = ".",
                                            coherence_threshold = 0.47,
                                            save_plots = TRUE,
                                            ecotype_annotations = NULL,
                                            coherence_matrix = NULL) {
  
  cat("Creating Plot 01: Mean Spatial Coherence by Ecotype (Line Plot)...\n")
  
  # Prepare data - calculate mean and standard deviation
  if (!is.null(coherence_matrix)) {
    if (is.data.frame(coherence_matrix) && "metaprogram" %in% colnames(coherence_matrix)) {
      # Enhanced analysis format
      matrix_data <- coherence_matrix[, !names(coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- coherence_matrix$metaprogram
      matrix_data <- as.matrix(matrix_data)
    } else {
      matrix_data <- as.matrix(coherence_matrix)
    }
    
    # Calculate from matrix
    mean_values <- rowMeans(matrix_data, na.rm = TRUE)
    sd_values <- apply(matrix_data, 1, sd, na.rm = TRUE)
  } else {
    # Use provided mean coherence data
    mean_values <- as.numeric(mean_coherence_data)
    names(mean_values) <- names(mean_coherence_data)
    sd_values <- rep(0, length(mean_values))  # No SD available
  }
  
  plot_data <- data.frame(
    metaprogram = names(mean_values),
    mean_spatial = mean_values,
    sd = sd_values,
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$mean_spatial), ]
  
  # Order by mean coherence (descending) - this matches your plot exactly
  plot_data <- plot_data[order(plot_data$mean_spatial, decreasing = TRUE), ]
  plot_data$metaprogram <- factor(plot_data$metaprogram, levels = plot_data$metaprogram)
  
  # Generate colors for ecotypes (matching your image colors)
  ecotype_colors <- generate_ecotype_colors(plot_data$metaprogram)
  
  # Create plot matching your Image exactly
  p1 <- ggplot(plot_data, aes(x = metaprogram, y = mean_spatial, group = 1)) + 
    geom_line(color = "grey", linewidth = 1.5) +
    geom_point(size = 4, aes(color = metaprogram)) +
    geom_errorbar(aes(ymin = pmax(mean_spatial - sd, 0), ymax = pmin(mean_spatial + sd, 1)), 
                  width = 0.2, linewidth = 0.8) + 
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = PLOT_CONFIG$axis_text_size, 
                                angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.position = "none",  # Remove legend as in your image
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values = ecotype_colors) +
    scale_y_continuous(
      limits = c(0, 1.0),
      breaks = seq(0, 1.0, by = 0.2),
      expand = c(0.02, 0.02)
    ) +
    labs(
      x = "metaprogram",
      y = "Mean Spatial Coherence Score", 
      title = "Mean Spatial Coherence by Ecotype (CORRECTED)"
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p1, "01_mean_coherence_by_ecotype_lineplot", output_subdir)
  }
  
  return(p1)
}

# ============================================================================
# PLOT 02: COHERENCE BY SAMPLE (TIROSH STYLE - MATCHING YOUR IMAGE 2)
# ============================================================================

create_02_coherence_by_sample_tirosh <- function(detailed_results = NULL,
                                               coherence_matrix = NULL,
                                               output_dir = ".",
                                               save_plots = TRUE,
                                               ecotype_annotations = NULL) {
  
  cat("Creating Plot 02: Spatial Coherence by Sample (Tirosh Style)...\n")
  
  # Prepare data from coherence matrix or detailed results
  if (!is.null(coherence_matrix)) {
    if (is.data.frame(coherence_matrix) && "metaprogram" %in% colnames(coherence_matrix)) {
      # Enhanced analysis format - convert to long format
      matrix_data <- coherence_matrix[, !names(coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- coherence_matrix$metaprogram
      spatial_long <- reshape2::melt(as.matrix(matrix_data))
      colnames(spatial_long) <- c("metaprogram", "sample", "spatial_score")
    } else {
      # Regular matrix format
      spatial_long <- reshape2::melt(coherence_matrix)
      colnames(spatial_long) <- c("metaprogram", "sample", "spatial_score")
    }
  } else if (!is.null(detailed_results)) {
    # Use detailed results directly
    spatial_long <- detailed_results
    colnames(spatial_long)[colnames(spatial_long) == "coherence_score"] <- "spatial_score"
  } else {
    stop("Either coherence_matrix or detailed_results must be provided")
  }
  
  # Remove NA values
  spatial_long <- spatial_long[!is.na(spatial_long$spatial_score), ]
  
  # Calculate sample statistics for ordering and error bars
  sample_stats <- spatial_long %>%
    group_by(sample) %>%
    summarise(
      mean_coherence = mean(spatial_score, na.rm = TRUE),
      sd_coherence = sd(spatial_score, na.rm = TRUE),
      n_ccs = n(),
      .groups = 'drop'
    ) %>%
    # Only calculate error bars for samples with at least 3 observations
    mutate(
      sd_coherence = ifelse(n_ccs >= 3 & !is.na(sd_coherence), sd_coherence, NA),
      # Handle cases where SD is 0 or very small
      sd_coherence = ifelse(sd_coherence < 0.01, NA, sd_coherence)
    ) %>%
    # Order samples by mean coherence (descending)
    arrange(desc(mean_coherence)) %>%
    mutate(sample = factor(sample, levels = sample))
  
  # Order the long data by sample coherence
  spatial_long$sample <- factor(spatial_long$sample, levels = levels(sample_stats$sample))
  
  # Generate colors for ecotypes (matching your image colors)
  unique_ecotypes <- sort(unique(spatial_long$metaprogram))
  ecotype_colors <- generate_ecotype_colors(unique_ecotypes)
  
  # Create plot matching your Image 2 exactly
  p2 <- ggplot() +
    # Add individual ecotype points (colored by ecotype)
    geom_point(data = spatial_long, 
               aes(x = sample, y = spatial_score, color = metaprogram),
               size = 2.5, alpha = 0.8) +
    
    # Add sample mean points (black dots)
    geom_point(data = sample_stats,
               aes(x = sample, y = mean_coherence),
               size = 1.5, color = "black", shape = 16) +
    
    # Add error bars ONLY for samples with sufficient data
    geom_errorbar(data = sample_stats %>% filter(!is.na(sd_coherence)),
                  aes(x = sample, 
                      ymin = pmax(mean_coherence - sd_coherence, 0), 
                      ymax = pmin(mean_coherence + sd_coherence, 1)),
                  width = 0.4, linewidth = 0.8, color = "black",
                  inherit.aes = FALSE) +
    
    # Styling to match your image
    theme_minimal() +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
      axis.text.y = element_text(size = 10, color = "black"),
      axis.title.x = element_text(size = 12, face = "bold", color = "black"),
      axis.title.y = element_text(size = 12, face = "bold", color = "black"),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "darkgray"),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.position = "right",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9)
    ) +
    
    # Color scale for ecotypes
    scale_color_manual(values = ecotype_colors, name = "Spatial Ecotype") +
    
    # Y-axis scale matching your image
    scale_y_continuous(
      limits = c(0.0, 0.8),
      breaks = seq(0.0, 0.8, by = 0.2),
      expand = c(0.02, 0.02)
    ) +
    
    # Labels matching your image
    labs(
      x = "sample",
      y = "spatial coherence",
      title = "Spatial coherence by sample (CORRECTED)",
      subtitle = paste("Analyzed", length(unique(spatial_long$sample)), "samples with", 
                      length(unique(spatial_long$metaprogram)), "ecotypes")
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    # Use wider format for this plot due to many samples
    ggsave(file.path(output_subdir, "pdf", "02_coherence_by_sample_tirosh_style.pdf"), 
           p2, width = 16, height = 6, device = "pdf", dpi = PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "png", "02_coherence_by_sample_tirosh_style.png"), 
           p2, width = 16, height = 6, device = "png", dpi = PLOT_CONFIG$dpi)
    cat("  ✓ Saved: 02_coherence_by_sample_tirosh_style in pdf, png format(s)\n")
  }
  
  return(p2)
}

# ============================================================================
# PLOT 03: ORGANIZATION CLASSIFICATION PIE CHART
# ============================================================================

create_03_organization_pie_chart <- function(organization_results,
                                            output_dir = ".",
                                            save_plots = TRUE) {
  
  cat("Creating Plot 03: Organization Classification Pie Chart...\n")
  
  # Prepare data
  org_counts <- table(organization_results)
  plot_data <- data.frame(
    organization = names(org_counts),
    count = as.numeric(org_counts),
    percentage = round(as.numeric(org_counts) / sum(org_counts) * 100, 1),
    stringsAsFactors = FALSE
  )
  
  # Calculate cumulative percentages for label positioning
  plot_data$cumsum <- cumsum(plot_data$count)
  plot_data$pos <- plot_data$cumsum - plot_data$count/2
  
  # Create plot
  p3 <- ggplot(plot_data, aes(x = "", y = count, fill = organization)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 2) +
    coord_polar("y", start = 0) +
    geom_text(aes(y = pos, label = paste0(organization, "\n", count, " (", percentage, "%)")),
              size = 4, fontface = "bold", color = "white") +
    scale_fill_manual(
      values = c("Organized" = PLOT_CONFIG$organized_color, 
                 "Disorganized" = PLOT_CONFIG$disorganized_color,
                 "Intermediate" = PLOT_CONFIG$intermediate_color),
      name = "Organization\nClassification"
    ) +
    theme_void() +
    theme(
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      legend.position = "bottom"
    ) +
    labs(
      title = "Ecotype Organization Classification",
      subtitle = paste("Total Ecotypes:", sum(plot_data$count))
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p3, "03_organization_classification_pie_chart", output_subdir)
  }
  
  return(p3)
}

# ============================================================================
# PLOT 04: COHERENCE DISTRIBUTION HISTOGRAM
# ============================================================================

create_04_coherence_distribution_histogram <- function(mean_coherence_data,
                                                      organization_results,
                                                      output_dir = ".",
                                                      coherence_threshold = 0.47,
                                                      save_plots = TRUE,
                                                      bins = 20) {
  
  cat("Creating Plot 04: Coherence Distribution Histogram...\n")
  
  # Prepare data
  coherence_values <- as.numeric(mean_coherence_data)
  coherence_values <- coherence_values[!is.na(coherence_values)]
  
  plot_data <- data.frame(
    coherence = coherence_values,
    organization = organization_results[names(mean_coherence_data)][!is.na(coherence_values)]
  )
  
  # Calculate statistics
  mean_coherence <- mean(plot_data$coherence)
  median_coherence <- median(plot_data$coherence)
  sd_coherence <- sd(plot_data$coherence)
  
  # Create plot
  p4 <- ggplot(plot_data, aes(x = coherence, fill = organization)) +
    geom_histogram(bins = bins, alpha = 0.7, color = "black", size = 0.3) +
    geom_vline(xintercept = coherence_threshold, 
               linetype = "dashed", color = "black", size = 1) +
    geom_vline(xintercept = mean_coherence, 
               linetype = "solid", color = "blue", size = 0.8, alpha = 0.7) +
    geom_vline(xintercept = median_coherence, 
               linetype = "dotted", color = "red", size = 0.8, alpha = 0.7) +
    scale_fill_manual(
      values = c("Organized" = PLOT_CONFIG$organized_color, 
                 "Disorganized" = PLOT_CONFIG$disorganized_color,
                 "Intermediate" = PLOT_CONFIG$intermediate_color),
      name = "Organization\nClassification"
    ) +
    theme_minimal() +
    theme(
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size)
    ) +
    labs(
      x = "Spatial Coherence Score",
      y = "Frequency",
      title = "Distribution of Spatial Coherence Scores",
      subtitle = paste0("Mean = ", round(mean_coherence, 3), 
                       " | Median = ", round(median_coherence, 3), 
                       " | SD = ", round(sd_coherence, 3))
    ) +
    # Add legend for lines
    annotate("text", x = Inf, y = Inf, 
             label = paste0("Threshold (", coherence_threshold, ")\n",
                           "Mean (", round(mean_coherence, 3), ")\n",
                           "Median (", round(median_coherence, 3), ")"),
             hjust = 1.1, vjust = 1.1, size = 3,
             color = "black", fontface = "bold")
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p4, "04_coherence_distribution_histogram", output_subdir)
  }
  
  return(p4)
}

# ============================================================================
# MAIN FUNCTION: CREATE ALL BASIC PLOTS
# ============================================================================

create_all_basic_plots <- function(results, 
                                  output_dir = "spatial_coherence_plots",
                                  save_plots = TRUE,
                                  coherence_threshold = NULL) {
  
  cat("=== SpatialCoherence Basic Plots Generation (Part 1/4) ===\n\n")
  
  # Extract threshold from results if not provided
  if (is.null(coherence_threshold)) {
    if (!is.null(results$analysis_parameters$coherence_threshold)) {
      coherence_threshold <- results$analysis_parameters$coherence_threshold
    } else {
      coherence_threshold <- 0.47  # Default
    }
  }
  
  cat("Using coherence threshold:", coherence_threshold, "\n\n")
  
  # Create only the plots that make sense
  plots <- list()
  
  # Plot 01: Mean Coherence Line Plot (the one you want)
  plots$plot_01 <- create_01_mean_coherence_lineplot(
    mean_coherence_data = results$mean_coherence,
    organization_results = results$organization_results,
    output_dir = output_dir,
    coherence_threshold = coherence_threshold,
    save_plots = save_plots,
    ecotype_annotations = results$ecotype_annotations,
    coherence_matrix = if(is.matrix(results$coherence_matrix) || 
                          (is.data.frame(results$coherence_matrix) && "metaprogram" %in% colnames(results$coherence_matrix))) 
                      results$coherence_matrix else NULL
  )
  
  # Plot 02: Coherence by Sample (Tirosh style)
  plots$plot_02 <- create_02_coherence_by_sample_tirosh(
    detailed_results = results$detailed_results,
    coherence_matrix = if(is.matrix(results$coherence_matrix) || 
                          (is.data.frame(results$coherence_matrix) && "metaprogram" %in% colnames(results$coherence_matrix))) 
                      results$coherence_matrix else NULL,
    output_dir = output_dir,
    save_plots = save_plots,
    ecotype_annotations = results$ecotype_annotations
  )
  
  # REMOVED: Plot 03 (pie chart) - not needed
  # REMOVED: Plot 04 (histogram) - not needed
  
  cat("\n=== Basic Plots Generation Complete (Part 1/4) ===\n")
  cat("Output directory:", file.path(output_dir, "01_basic_plots"), "\n")
  cat("Number of plots created:", length(plots), "\n")
  cat("Plots created:\n")
  cat("  01: Mean Coherence by Ecotype (Line Plot)\n")
  cat("  02: Coherence by Sample (Tirosh Style)\n")
  cat("Note: Pie chart and histogram plots removed as not needed for this analysis.\n")
  
  return(plots)
}

# ============================================================================
# EXAMPLE USAGE
# ============================================================================

# Example usage (commented out):
# 
# # Assuming you have results from run_spatial_analysis()
# basic_plots <- create_all_basic_plots(
#   results = your_spatial_analysis_results,
#   output_dir = "output/spatial_coherence_plots",
#   save_plots = TRUE
# )

cat("SpatialCoherence Basic Plotting Suite (Part 1/4) loaded successfully.\n")
cat("Use create_all_basic_plots() to generate all basic visualization plots.\n")