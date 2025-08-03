#!/usr/bin/env Rscript

# ============================================================================
# SpatialCoherence Plotting Suite - Part 3: Segmentation Analysis Plots
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script generates segmentation analysis visualization plots
# Part 3 of 4: Cell density and functional classification plots

# Required libraries
required_packages <- c("ggplot2", "dplyr", "reshape2", "patchwork", "viridis")

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

# Source basic plotting config if not already loaded
if (!exists("PLOT_CONFIG")) {
  PLOT_CONFIG <- list(
    width = 12, height = 8, dpi = 300,
    organized_color = "#2E8B57", disorganized_color = "#DC143C", 
    intermediate_color = "#FFD700",
    title_size = 16, subtitle_size = 12, axis_title_size = 14, 
    axis_text_size = 12, legend_text_size = 10,
    formats = c("pdf", "png")
  )
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

# Generate functional classification colors
generate_functional_colors <- function(functional_types) {
  standard_colors <- c(
    "Hypoxic" = "#8B0000", "Hypoxic Core" = "#8B0000",
    "Proliferative" = "#4169E1", 
    "Stromal" = "#228B22",
    "Immune" = "#FF6347",
    "Vascular" = "#9370DB",
    "Normal" = "#DAA520",
    "Other" = "#808080", "Unknown" = "#808080"
  )
  
  colors <- character(length(functional_types))
  names(colors) <- functional_types
  
  for (i in seq_along(functional_types)) {
    type <- functional_types[i]
    if (type %in% names(standard_colors)) {
      colors[type] <- standard_colors[type]
    } else {
      colors[type] <- rainbow(length(functional_types))[i]
    }
  }
  
  return(colors)
}

# Create output directory structure with user-defined path
create_output_structure <- function(base_dir) {
  # Ensure base directory is created first
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  dirs <- c(
    file.path(base_dir, "03_segmentation_plots"),
    file.path(base_dir, "03_segmentation_plots", "pdf"),
    file.path(base_dir, "03_segmentation_plots", "png")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(dirs[1])
}

# Save plot in multiple formats
save_plot_formats <- function(plot, filename, output_dir, width = NULL, height = NULL) {
  w <- ifelse(is.null(width), PLOT_CONFIG$width, width)
  h <- ifelse(is.null(height), PLOT_CONFIG$height, height)
  
  for (format in PLOT_CONFIG$formats) {
    filepath <- file.path(output_dir, format, paste0(filename, ".", format))
    ggsave(filepath, plot, width = w, height = h, device = format, dpi = PLOT_CONFIG$dpi)
  }
  
  cat("  ✓ Saved:", filename, "in", paste(PLOT_CONFIG$formats, collapse = ", "), "format(s)\n")
}

# ============================================================================
# PLOT 09: CELL DENSITY BY TREATMENT STATUS (MATCHING YOUR IMAGE 7 - LEFT PANEL)
# ============================================================================

create_09_cell_density_by_treatment <- function(density_data = NULL,
                                              detailed_results = NULL,
                                              treatment_column = "treatment",
                                              output_dir = ".",
                                              save_plots = TRUE) {
  
  cat("Creating Plot 09: Cell Density by Treatment Status...\n")
  
  # Prepare data
  if (!is.null(density_data)) {
    plot_data <- density_data
  } else if (!is.null(detailed_results)) {
    # Calculate sample-level cell density from detailed results
    plot_data <- detailed_results %>%
      group_by(sample) %>%
      summarise(cell_count_per_spot = n(), .groups = 'drop')
    
    # Add treatment information if available
    if (treatment_column %in% colnames(detailed_results)) {
      treatment_info <- detailed_results %>%
        select(sample, !!sym(treatment_column)) %>%
        distinct()
      
      plot_data <- merge(plot_data, treatment_info, by = "sample")
      colnames(plot_data)[colnames(plot_data) == treatment_column] <- "treatment_status"
    } else {
      stop("Treatment column not found in data")
    }
  } else {
    stop("Either density_data or detailed_results must be provided")
  }
  
  # Remove missing values
  plot_data <- plot_data[!is.na(plot_data$treatment_status), ]
  
  # Ensure proper factor ordering
  plot_data$treatment_status <- factor(plot_data$treatment_status, 
                                      levels = c("Untreated", "Treated"))
  
  # Perform statistical test
  if (length(unique(plot_data$treatment_status)) == 2) {
    test_result <- t.test(cell_count_per_spot ~ treatment_status, data = plot_data)
    p_value <- format.pval(test_result$p.value, digits = 3)
  } else {
    p_value <- "N/A"
  }
  
  # Create plot matching your Image 7 (left panel)
  p9 <- ggplot(plot_data, aes(x = treatment_status, y = cell_count_per_spot)) +
    geom_boxplot(aes(fill = treatment_status), alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = c("Untreated" = "#11A579", "Treated" = "#CF1C90")) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      legend.position = "none",
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Treatment Status",
      y = "cell count per spot",
      title = "Cell Density by Treatment Status",
      subtitle = "CORRECTED: Naive vs Treated"
    )
  
  # Add p-value annotation if available
  if (p_value != "N/A") {
    p9 <- p9 + 
      annotate("text", x = 1.5, y = max(plot_data$cell_count_per_spot) * 0.9, 
               label = paste("p =", p_value), size = 3.5, fontface = "bold")
  }
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p9, "09_cell_density_by_treatment_status", output_subdir)
  }
  
  return(p9)
}

# ============================================================================
# PLOT 10: CELL DENSITY BY ORGANIZATION ZONE (MATCHING YOUR IMAGE 7 - MIDDLE PANEL)
# ============================================================================

create_10_cell_density_by_organization <- function(density_data = NULL,
                                                  detailed_results = NULL,
                                                  coherence_threshold = 0.47,
                                                  output_dir = ".",
                                                  save_plots = TRUE) {
  
  cat("Creating Plot 10: Cell Density by Organization Zone...\n")
  
  # Prepare data
  if (!is.null(density_data)) {
    plot_data <- density_data
  } else if (!is.null(detailed_results)) {
    # Calculate sample-level organization from detailed results
    sample_coherence <- detailed_results %>%
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
    
    plot_data <- sample_coherence
  } else {
    stop("Either density_data or detailed_results must be provided")
  }
  
  # Remove missing values
  plot_data <- plot_data[!is.na(plot_data$organization_zone), ]
  
  # Ensure proper factor ordering
  plot_data$organization_zone <- factor(plot_data$organization_zone, 
                                       levels = c("structured", "intermediate", "disorganized"))
  
  # Perform statistical tests
  p_values <- list()
  if ("structured" %in% plot_data$organization_zone && "disorganized" %in% plot_data$organization_zone) {
    struct_vs_disorg <- plot_data %>%
      filter(organization_zone %in% c("structured", "disorganized"))
    
    if (nrow(struct_vs_disorg) > 0) {
      test_result <- t.test(cell_count_per_spot ~ organization_zone, data = struct_vs_disorg)
      p_values$struct_vs_disorg <- format.pval(test_result$p.value, digits = 3)
    }
  }
  
  if ("intermediate" %in% plot_data$organization_zone && "disorganized" %in% plot_data$organization_zone) {
    inter_vs_disorg <- plot_data %>%
      filter(organization_zone %in% c("intermediate", "disorganized"))
    
    if (nrow(inter_vs_disorg) > 0) {
      test_result <- t.test(cell_count_per_spot ~ organization_zone, data = inter_vs_disorg)
      p_values$inter_vs_disorg <- format.pval(test_result$p.value, digits = 3)
    }
  }
  
  # Create plot matching your Image 7 (middle panel)
  p10 <- ggplot(plot_data, aes(x = organization_zone, y = cell_count_per_spot)) +
    geom_boxplot(aes(fill = organization_zone), alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = c("structured" = PLOT_CONFIG$organized_color, 
                                "intermediate" = PLOT_CONFIG$intermediate_color, 
                                "disorganized" = PLOT_CONFIG$disorganized_color)) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      legend.position = "none",
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Organization Zone",
      y = "cell count per spot",
      title = "Cell Density by Organization Zone",
      subtitle = "CORRECTED: Coherence-based thresholds"
    )
  
  # Add p-value annotations if available
  if ("struct_vs_disorg" %in% names(p_values)) {
    p10 <- p10 + 
      annotate("text", x = 1.5, y = max(plot_data$cell_count_per_spot) * 0.95, 
               label = paste("p =", p_values$struct_vs_disorg), size = 3.5, fontface = "bold")
  }
  
  if ("inter_vs_disorg" %in% names(p_values)) {
    p10 <- p10 + 
      annotate("text", x = 2.5, y = max(plot_data$cell_count_per_spot) * 0.85, 
               label = paste("p =", p_values$inter_vs_disorg), size = 3.5, fontface = "bold")
  }
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p10, "10_cell_density_by_organization_zone", output_subdir)
  }
  
  return(p10)
}

# ============================================================================
# PLOT 11: CELL DENSITY BY FUNCTIONAL CLASSIFICATION (MATCHING YOUR IMAGE 7 - RIGHT PANEL)
# ============================================================================

create_11_cell_density_by_functional <- function(density_data = NULL,
                                                detailed_results = NULL,
                                                functional_annotations = NULL,
                                                output_dir = ".",
                                                save_plots = TRUE) {
  
  cat("Creating Plot 11: Cell Density by Functional Classification...\n")
  
  # Prepare data
  if (!is.null(density_data)) {
    plot_data <- density_data
  } else if (!is.null(detailed_results)) {
    # Add functional classification based on annotations or infer from ecotype names
    if (!is.null(functional_annotations)) {
      detailed_results$ecotype_function <- functional_annotations[detailed_results$metaprogram]
    } else {
      # Try to infer functional types from ecotype names (basic heuristic)
      detailed_results$ecotype_function <- sapply(detailed_results$metaprogram, function(x) {
        x_lower <- tolower(as.character(x))
        if (grepl("hypox|core", x_lower)) return("Hypoxic Core")
        if (grepl("prolif|cycling", x_lower)) return("Proliferative")
        if (grepl("strom|fibro|caf", x_lower)) return("Stromal")
        if (grepl("immune|tcell|bcell|macro", x_lower)) return("Immune")
        if (grepl("vasc|endo", x_lower)) return("Vascular")
        if (grepl("normal|acin", x_lower)) return("Normal")
        return("Other")
      })
    }
    
    # Calculate sample-level density by functional type
    plot_data <- detailed_results %>%
      group_by(sample, ecotype_function) %>%
      summarise(cell_count_per_spot = n(), .groups = 'drop') %>%
      filter(!is.na(ecotype_function) & ecotype_function != "Other")
    
  } else {
    stop("Either density_data or detailed_results must be provided")
  }
  
  # Remove missing values and order by biological hierarchy
  plot_data <- plot_data[!is.na(plot_data$ecotype_function), ]
  
  functional_levels <- c("Hypoxic Core", "Proliferative", "Stromal", "Immune", "Vascular", "Normal")
  available_levels <- intersect(functional_levels, unique(plot_data$ecotype_function))
  plot_data$ecotype_function <- factor(plot_data$ecotype_function, levels = available_levels)
  
  # Generate functional colors
  functional_colors <- generate_functional_colors(available_levels)
  
  # Create plot matching your Image 7 (right panel)
  p11 <- ggplot(plot_data, aes(x = ecotype_function, y = cell_count_per_spot)) +
    geom_boxplot(aes(fill = ecotype_function), alpha = 0.7, outlier.size = 0.8) +
    geom_jitter(width = 0.2, size = 0.5, alpha = 0.6) +
    scale_fill_manual(values = functional_colors) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text.x = element_text(size = PLOT_CONFIG$axis_text_size, angle = 45, hjust = 1),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      legend.position = "none",
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Functional Ecotype",
      y = "cell count per spot",
      title = "Cell Density by Functional Classification",
      subtitle = "CORRECTED: Manuscript-based ecotypes"
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p11, "11_cell_density_by_functional_classification", output_subdir, width = 10, height = 5)
  }
  
  return(p11)
}

# ============================================================================
# PLOT 12: COMBINED SEGMENTATION ANALYSIS (MATCHING YOUR IMAGE 7 - ALL PANELS)
# ============================================================================

create_12_combined_segmentation_analysis <- function(density_data = NULL,
                                                    detailed_results = NULL,
                                                    treatment_column = "treatment",
                                                    coherence_threshold = 0.47,
                                                    functional_annotations = NULL,
                                                    output_dir = ".",
                                                    save_plots = TRUE) {
  
  cat("Creating Plot 12: Combined Segmentation Analysis (Three Panels)...\n")
  
  # Create individual panels
  p9 <- create_09_cell_density_by_treatment(
    density_data = density_data,
    detailed_results = detailed_results,
    treatment_column = treatment_column,
    output_dir = output_dir,
    save_plots = FALSE
  )
  
  p10 <- create_10_cell_density_by_organization(
    density_data = density_data,
    detailed_results = detailed_results,
    coherence_threshold = coherence_threshold,
    output_dir = output_dir,
    save_plots = FALSE
  )
  
  p11 <- create_11_cell_density_by_functional(
    density_data = density_data,
    detailed_results = detailed_results,
    functional_annotations = functional_annotations,
    output_dir = output_dir,
    save_plots = FALSE
  )
  
  # Combine panels using patchwork
  combined_plot <- p9 | p10 | p11
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(combined_plot, "12_combined_segmentation_analysis", output_subdir, width = 20, height = 5)
  }
  
  return(combined_plot)
}

# ============================================================================
# MAIN FUNCTION: CREATE ALL SEGMENTATION PLOTS
# ============================================================================

create_all_segmentation_plots <- function(results, 
                                         output_dir = "spatial_coherence_plots",
                                         save_plots = TRUE,
                                         treatment_column = "treatment",
                                         coherence_threshold = NULL) {
  
  cat("=== SpatialCoherence Segmentation Analysis Plots (Part 3/4) ===\n\n")
  
  # Extract threshold from results if not provided
  if (is.null(coherence_threshold)) {
    if (!is.null(results$analysis_parameters$coherence_threshold)) {
      coherence_threshold <- results$analysis_parameters$coherence_threshold
    } else {
      coherence_threshold <- 0.47  # Default
    }
  }
  
  plots <- list()
  
  # Check if detailed results are available
  if (is.null(results$detailed_results)) {
    cat("No detailed results available. Skipping segmentation plots.\n")
    return(plots)
  }
  
  # Check if treatment column exists
  has_treatment_data <- treatment_column %in% colnames(results$detailed_results)
  
  if (!has_treatment_data) {
    cat("No treatment data available in detailed results.\n")
    cat("Skipping treatment-based segmentation plots.\n")
    cat("Available columns:", paste(colnames(results$detailed_results), collapse = ", "), "\n")
  }
  
  # Only create plots that don't require treatment data
  # Plot 10: Cell Density by Organization Zone (this doesn't need treatment data)
  tryCatch({
    plots$plot_10 <- create_10_cell_density_by_organization(
      detailed_results = results$detailed_results,
      coherence_threshold = coherence_threshold,
      output_dir = output_dir,
      save_plots = save_plots
    )
  }, error = function(e) {
    cat("  ! Skipped Plot 10 (Cell Density by Organization):", e$message, "\n")
  })
  
  # Plot 11: Cell Density by Functional Classification (this doesn't need treatment data)
  tryCatch({
    plots$plot_11 <- create_11_cell_density_by_functional(
      detailed_results = results$detailed_results,
      functional_annotations = results$ecotype_annotations,
      output_dir = output_dir,
      save_plots = save_plots
    )
  }, error = function(e) {
    cat("  ! Skipped Plot 11 (Cell Density by Functional):", e$message, "\n")
  })
  
  # Only create treatment-based plots if treatment data exists
  if (has_treatment_data) {
    # Plot 09: Cell Density by Treatment Status
    tryCatch({
      plots$plot_09 <- create_09_cell_density_by_treatment(
        detailed_results = results$detailed_results,
        treatment_column = treatment_column,
        output_dir = output_dir,
        save_plots = save_plots
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 09 (Cell Density by Treatment):", e$message, "\n")
    })
    
    # Plot 12: Combined Segmentation Analysis
    tryCatch({
      plots$plot_12 <- create_12_combined_segmentation_analysis(
        detailed_results = results$detailed_results,
        treatment_column = treatment_column,
        coherence_threshold = coherence_threshold,
        functional_annotations = results$ecotype_annotations,
        output_dir = output_dir,
        save_plots = save_plots
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 12 (Combined Segmentation):", e$message, "\n")
    })
  }
  
  if (length(plots) > 0) {
    cat("\n=== Segmentation Analysis Plots Complete (Part 3/4) ===\n")
    cat("Output directory:", file.path(output_dir, "03_segmentation_plots"), "\n")
    cat("Number of plots created:", length(plots), "\n")
    cat("Plots created:\n")
    if ("plot_09" %in% names(plots)) cat("  09: Cell Density by Treatment Status\n")
    if ("plot_10" %in% names(plots)) cat("  10: Cell Density by Organization Zone\n")
    if ("plot_11" %in% names(plots)) cat("  11: Cell Density by Functional Classification\n")
    if ("plot_12" %in% names(plots)) cat("  12: Combined Segmentation Analysis\n")
  } else {
    cat("No segmentation plots were created.\n")
  }
  
  return(plots)
}

cat("SpatialCoherence Segmentation Analysis Plotting Suite (Part 3/4) loaded successfully.\n")
cat("Use create_all_segmentation_plots() to generate all segmentation analysis plots.\n") 