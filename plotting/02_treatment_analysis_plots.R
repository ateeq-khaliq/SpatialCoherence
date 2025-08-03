#!/usr/bin/env Rscript

# ============================================================================
# SpatialCoherence Plotting Suite - Part 2: Treatment Analysis Plots
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script generates treatment analysis visualization plots
# Part 2 of 4: Treatment effects and sample organization plots

# Required libraries
required_packages <- c("ggplot2", "dplyr", "reshape2", "patchwork", "viridis", "RColorBrewer")

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
  # Standard functional colors mapping
  standard_colors <- c(
    "Hypoxic" = "#8B0000", "Hypoxic Core" = "#8B0000",
    "Proliferative" = "#4169E1", 
    "Stromal" = "#228B22",
    "Immune" = "#FF6347",
    "Vascular" = "#9370DB",
    "Normal" = "#DAA520",
    "Other" = "#808080", "Unknown" = "#808080"
  )
  
  # Use standard colors where available, generate others
  colors <- character(length(functional_types))
  names(colors) <- functional_types
  
  for (i in seq_along(functional_types)) {
    type <- functional_types[i]
    if (type %in% names(standard_colors)) {
      colors[type] <- standard_colors[type]
    } else {
      # Generate unique color for unknown types
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
    file.path(base_dir, "02_treatment_plots"),
    file.path(base_dir, "02_treatment_plots", "pdf"),
    file.path(base_dir, "02_treatment_plots", "png")
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
# PLOT 05: SAMPLE ORGANIZATION BY TREATMENT (MATCHING YOUR IMAGE 3)
# ============================================================================

create_05_sample_organization_by_treatment <- function(sample_data = NULL,
                                                     detailed_results = NULL,
                                                     treatment_column = "treatment",
                                                     coherence_threshold = 0.47,
                                                     output_dir = ".",
                                                     save_plots = TRUE) {
  
  cat("Creating Plot 05: Sample Organization Distribution by Treatment...\n")
  
  # Prepare sample-level organization data
  if (!is.null(sample_data)) {
    plot_data <- sample_data
  } else if (!is.null(detailed_results)) {
    # Calculate sample-level coherence from detailed results
    sample_coherence <- detailed_results %>%
      group_by(sample) %>%
      summarise(
        mean_coherence = mean(coherence_score, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      mutate(
        organization_zone = case_when(
          mean_coherence >= coherence_threshold + 0.08 ~ "structured",
          mean_coherence >= coherence_threshold - 0.08 ~ "intermediate", 
          TRUE ~ "disorganized"
        )
      )
    
    # Add treatment information if available
    if (treatment_column %in% colnames(detailed_results)) {
      treatment_info <- detailed_results %>%
        select(sample, !!sym(treatment_column)) %>%
        distinct()
      
      plot_data <- merge(sample_coherence, treatment_info, by = "sample")
      colnames(plot_data)[colnames(plot_data) == treatment_column] <- "treatment_status"
    } else {
      stop("Treatment column not found in data")
    }
  } else {
    stop("Either sample_data or detailed_results must be provided")
  }
  
  # Calculate proportions for the stacked bar chart
  summary_data <- plot_data %>%
    group_by(treatment_status, organization_zone) %>%
    summarise(count = n(), .groups = 'drop') %>%
    group_by(treatment_status) %>%
    mutate(
      total = sum(count),
      proportion = count / total,
      percentage = round(proportion * 100, 1)
    )
  
  # Perform statistical test
  contingency_table <- table(plot_data$treatment_status, plot_data$organization_zone)
  if (all(contingency_table >= 5)) {
    test_result <- chisq.test(contingency_table)
    p_value <- test_result$p.value
  } else {
    test_result <- fisher.test(contingency_table, simulate.p.value = TRUE)
    p_value <- test_result$p.value
  }
  
  # Create plot matching your Image 3
  p5 <- ggplot(summary_data, aes(x = treatment_status, y = proportion, fill = organization_zone)) +
    geom_bar(stat = "identity", position = "stack", alpha = 0.8, color = "white", size = 0.5) +
    geom_text(aes(label = paste0(percentage, "%")), 
              position = position_stack(vjust = 0.5), 
              size = 4, fontface = "bold", color = "white") +
    scale_fill_manual(
      values = c("structured" = PLOT_CONFIG$organized_color, 
                 "intermediate" = PLOT_CONFIG$intermediate_color, 
                 "disorganized" = PLOT_CONFIG$disorganized_color),
      name = "Organization Zone"
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = scales::percent_format()
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Treatment Status",
      y = "Proportion of Samples",
      title = "Sample Organization Distribution by Treatment (CORRECTED)",
      subtitle = paste("p =", format.pval(p_value, digits = 3))
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p5, "05_sample_organization_by_treatment", output_subdir)
  }
  
  return(p5)
}

# ============================================================================
# PLOT 06: ALL ECOTYPES COHERENCE BY TREATMENT (MATCHING YOUR IMAGE 4)
# ============================================================================

create_06_ecotypes_coherence_by_treatment <- function(detailed_results = NULL,
                                                    treatment_column = "treatment",
                                                    output_dir = ".",
                                                    save_plots = TRUE,
                                                    functional_annotations = NULL) {
  
  cat("Creating Plot 06: All Ecotypes Coherence by Treatment...\n")
  
  if (is.null(detailed_results)) {
    stop("detailed_results must be provided")
  }
  
  # Prepare data
  if (treatment_column %in% colnames(detailed_results)) {
    plot_data <- detailed_results
    colnames(plot_data)[colnames(plot_data) == treatment_column] <- "treatment_status"
  } else {
    stop("Treatment column not found in detailed_results")
  }
  
  # Remove missing values
  plot_data <- plot_data[!is.na(plot_data$treatment_status) & 
                        !is.na(plot_data$coherence_score), ]
  
  # Get all unique ecotypes and order them
  all_ecotypes <- sort(unique(plot_data$metaprogram))
  plot_data$metaprogram <- factor(plot_data$metaprogram, levels = all_ecotypes)
  
  # Perform statistical tests for each ecotype
  p_values <- sapply(all_ecotypes, function(ecotype) {
    ecotype_data <- plot_data[plot_data$metaprogram == ecotype, ]
    if (length(unique(ecotype_data$treatment_status)) == 2 && nrow(ecotype_data) >= 10) {
      test_result <- t.test(coherence_score ~ treatment_status, data = ecotype_data)
      return(test_result$p.value)
    } else {
      return(NA)
    }
  })
  
  # Create plot matching your Image 4
  p6 <- ggplot(plot_data, aes(x = metaprogram, y = coherence_score, fill = treatment_status)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.8, color = "black", size = 0.3) +
    geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.2), 
               size = 0.8, alpha = 0.6) +
    scale_fill_manual(
      values = c("Treated" = "#CF1C90", "Untreated" = "#11A579"),
      name = "Treatment Status"
    ) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Ecotype",
      y = "Spatial Coherence Score",
      title = paste("All Ecotypes Coherence by Treatment (", paste(all_ecotypes, collapse = "-"), ")", sep = ""),
      subtitle = "Statistical significance shown above each ecotype"
    )
  
  # Add p-value annotations
  y_max <- max(plot_data$coherence_score, na.rm = TRUE) * 1.15
  
  for (i in seq_along(all_ecotypes)) {
    if (!is.na(p_values[i])) {
      p_val <- format.pval(p_values[i], digits = 3)
      p_color <- ifelse(p_values[i] < 0.05, "red", "black")
      p_text <- ifelse(p_values[i] < 0.05, paste("p =", p_val, "*"), paste("p =", p_val))
      
      p6 <- p6 + 
        annotate("text", x = i, y = y_max, 
                 label = p_text, size = 2.5, fontface = "bold", color = p_color)
    }
  }
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p6, "06_all_ecotypes_coherence_by_treatment", output_subdir, width = 14, height = 8)
  }
  
  return(p6)
}

# ============================================================================
# PLOT 07: TREATMENT EFFECT MAGNITUDE (MATCHING YOUR IMAGE 5)
# ============================================================================

create_07_treatment_effect_magnitude <- function(detailed_results = NULL,
                                               treatment_column = "treatment",
                                               output_dir = ".",
                                               save_plots = TRUE,
                                               functional_annotations = NULL) {
  
  cat("Creating Plot 07: Treatment Effect Magnitude for All Ecotypes...\n")
  
  if (is.null(detailed_results)) {
    stop("detailed_results must be provided")
  }
  
  # Prepare data
  if (treatment_column %in% colnames(detailed_results)) {
    plot_data <- detailed_results
    colnames(plot_data)[colnames(plot_data) == treatment_column] <- "treatment_status"
  } else {
    stop("Treatment column not found in detailed_results")
  }
  
  # Remove missing values
  plot_data <- plot_data[!is.na(plot_data$treatment_status) & 
                        !is.na(plot_data$coherence_score), ]
  
  # Calculate effect sizes for each ecotype
  all_ecotypes <- unique(plot_data$metaprogram)
  effect_data <- data.frame()
  
  for (ecotype in all_ecotypes) {
    ecotype_subset <- plot_data[plot_data$metaprogram == ecotype, ]
    
    if (length(unique(ecotype_subset$treatment_status)) == 2 && nrow(ecotype_subset) >= 10) {
      
      test_result <- t.test(coherence_score ~ treatment_status, data = ecotype_subset)
      
      # Calculate effect size (Cohen's d)
      treated_mean <- mean(ecotype_subset$coherence_score[ecotype_subset$treatment_status == "Treated"], na.rm = TRUE)
      untreated_mean <- mean(ecotype_subset$coherence_score[ecotype_subset$treatment_status == "Untreated"], na.rm = TRUE)
      pooled_sd <- sqrt(((sum(ecotype_subset$treatment_status == "Treated") - 1) * 
                        var(ecotype_subset$coherence_score[ecotype_subset$treatment_status == "Treated"], na.rm = TRUE) + 
                        (sum(ecotype_subset$treatment_status == "Untreated") - 1) * 
                        var(ecotype_subset$coherence_score[ecotype_subset$treatment_status == "Untreated"], na.rm = TRUE)) / 
                       (nrow(ecotype_subset) - 2))
      
      cohens_d <- (treated_mean - untreated_mean) / pooled_sd
      
      # Add functional annotation if available
      if (!is.null(functional_annotations) && ecotype %in% names(functional_annotations)) {
        functional_type <- functional_annotations[ecotype]
      } else {
        functional_type <- "Other"
      }
      
      effect_data <- rbind(effect_data, data.frame(
        ecotype = ecotype,
        effect_size = cohens_d,
        p_value = test_result$p.value,
        significant = test_result$p.value < 0.05,
        functional_type = functional_type,
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # Order by effect size
  effect_data <- effect_data[order(abs(effect_data$effect_size), decreasing = TRUE), ]
  effect_data$ecotype <- factor(effect_data$ecotype, levels = effect_data$ecotype)
  
  # Generate functional colors
  functional_types <- unique(effect_data$functional_type)
  functional_colors <- generate_functional_colors(functional_types)
  
  # Create plot matching your Image 5
  p7 <- ggplot(effect_data, aes(x = ecotype, y = effect_size, fill = functional_type, alpha = significant)) +
    geom_col(width = 0.7, color = "black", size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
    geom_hline(yintercept = c(-0.2, 0.2), linetype = "dotted", color = "grey", alpha = 0.7) +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dotted", color = "grey", alpha = 0.5) +
    scale_fill_manual(values = functional_colors, name = "Functional Type") +
    scale_alpha_manual(values = c("TRUE" = 1.0, "FALSE" = 0.4), 
                      name = "Significant", labels = c("No", "Yes")) +
    theme_minimal() +
    theme(
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      plot.caption = element_text(size = 8, hjust = 0.5),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Ecotypes (Ranked by Effect Size)",
      y = "Treatment Effect Size (Cohen's d)",
      title = "Treatment Effect Magnitude for All Ecotypes",
      subtitle = "Positive = Higher coherence in treated; Negative = Lower coherence in treated",
      caption = "Dotted lines: Small (±0.2) and Medium (±0.5) effect thresholds"
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p7, "07_treatment_effect_magnitude", output_subdir, width = 14, height = 8)
  }
  
  return(p7)
}

# ============================================================================
# PLOT 08: COMPOSITIONAL DIFFERENCES (MATCHING YOUR IMAGE 6)
# ============================================================================

create_08_compositional_differences <- function(compositional_results = NULL,
                                              output_dir = ".",
                                              save_plots = TRUE) {
  
  cat("Creating Plot 08: Compositional Differences - Structured vs Disorganized Regions...\n")
  
  if (is.null(compositional_results)) {
    stop("compositional_results must be provided")
  }
  
  # Extract the differences data
  if (is.list(compositional_results) && "differences_data" %in% names(compositional_results)) {
    differences_data <- compositional_results$differences_data
  } else {
    differences_data <- compositional_results
  }
  
  # Order by relative abundance for better visualization
  differences_data <- differences_data[order(differences_data$rel_abundance), ]
  differences_data$metaprogram <- factor(differences_data$metaprogram, levels = differences_data$metaprogram)
  
  # Create color gradient
  n_colors <- nrow(differences_data)
  colors <- colorRampPalette(c("#8B1538", "#CD5C5C", "#F4A460", "#F0E68C", 
                              "#98FB98", "#66CDAA", "#4682B4", "#6A5ACD", "#9370DB"))(n_colors)
  color_mapping <- setNames(colors, differences_data$metaprogram)
  
  # Create plot matching your Image 6
  p8 <- ggplot(differences_data, aes(x = metaprogram, y = rel_abundance)) +
    geom_col(aes(fill = metaprogram), alpha = 0.8, width = 0.7, color = "black", size = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
    scale_fill_manual(values = color_mapping) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = PLOT_CONFIG$axis_text_size, angle = 45, hjust = 1, vjust = 1),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title.x = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      axis.title.y = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "metaprogram",
      y = "rel. abundance", 
      title = "Compositional Differences: Structured vs Disorganized Regions"
    )
  
  # Add significance annotations
  for (i in 1:nrow(differences_data)) {
    if (!is.na(differences_data$significance[i]) && differences_data$significance[i] != "") {
      
      # Position for annotation
      y_pos <- differences_data$rel_abundance[i]
      if (y_pos > 0) {
        y_pos <- y_pos + max(abs(differences_data$rel_abundance)) * 0.05
      } else {
        y_pos <- y_pos - max(abs(differences_data$rel_abundance)) * 0.05
      }
      
      # Add significance stars
      p8 <- p8 + annotate("text", x = i, y = y_pos, 
                         label = differences_data$significance[i], 
                         size = 4, fontface = "bold", color = "black")
      
      # Add p-value text if available
      if (!is.na(differences_data$p_value[i])) {
        p_text <- ifelse(differences_data$p_value[i] < 0.001,
                        sprintf("p=%.2e", differences_data$p_value[i]),
                        sprintf("p=%.6f", differences_data$p_value[i]))
        
        if (y_pos > 0) {
          p_y_pos <- y_pos + max(abs(differences_data$rel_abundance)) * 0.03
        } else {
          p_y_pos <- y_pos - max(abs(differences_data$rel_abundance)) * 0.03
        }
        
        p8 <- p8 + annotate("text", x = i, y = p_y_pos,
                           label = p_text, size = 3, color = "black")
      }
    }
  }
  
  # Add region labels
  max_val <- max(abs(differences_data$rel_abundance))
  p8 <- p8 + 
    annotate("text", x = nrow(differences_data) * 0.8, y = max_val * 0.8, 
             label = "struct", size = 4, fontface = "bold", color = "darkgreen") +
    annotate("text", x = nrow(differences_data) * 0.2, y = -max_val * 0.8,
             label = "disorg", size = 4, fontface = "bold", color = "darkred")
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    save_plot_formats(p8, "08_compositional_differences_structured_vs_disorganized", output_subdir)
  }
  
  return(p8)
}

# ============================================================================
# MAIN FUNCTION: CREATE ALL TREATMENT PLOTS
# ============================================================================

create_all_treatment_plots <- function(results, 
                                      output_dir = "spatial_coherence_plots",
                                      save_plots = TRUE,
                                      treatment_column = "treatment") {
  
  cat("=== SpatialCoherence Treatment Analysis Plots (Part 2/4) ===\n\n")
  
  plots <- list()
  
  # Check if treatment results are available and treatment column exists
  has_treatment_column <- FALSE
  if (!is.null(results$detailed_results)) {
    has_treatment_column <- treatment_column %in% colnames(results$detailed_results)
  }
  
  if (is.null(results$treatment_results) && !has_treatment_column) {
    cat("No treatment data available. Skipping treatment plots.\n")
    cat("Treatment analysis plots require either:\n")
    cat("  - treatment_results in the results object, OR\n")
    cat("  - '", treatment_column, "' column in detailed_results\n")
    return(plots)
  }
  
  # Only proceed if we have treatment data
  if (has_treatment_column) {
    # Plot 05: Sample Organization by Treatment
    tryCatch({
      plots$plot_05 <- create_05_sample_organization_by_treatment(
        detailed_results = results$detailed_results,
        treatment_column = treatment_column,
        output_dir = output_dir,
        save_plots = save_plots
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 05 (Sample Organization by Treatment):", e$message, "\n")
    })
    
    # Plot 06: All Ecotypes Coherence by Treatment
    tryCatch({
      plots$plot_06 <- create_06_ecotypes_coherence_by_treatment(
        detailed_results = results$detailed_results,
        treatment_column = treatment_column,
        output_dir = output_dir,
        save_plots = save_plots,
        functional_annotations = results$ecotype_annotations
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 06 (Ecotypes by Treatment):", e$message, "\n")
    })
    
    # Plot 07: Treatment Effect Magnitude
    tryCatch({
      plots$plot_07 <- create_07_treatment_effect_magnitude(
        detailed_results = results$detailed_results,
        treatment_column = treatment_column,
        output_dir = output_dir,
        save_plots = save_plots,
        functional_annotations = results$ecotype_annotations
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 07 (Treatment Effect Magnitude):", e$message, "\n")
    })
  }
  
  # Plot 08: Compositional Differences (independent of treatment column)
  if (!is.null(results$compositional_results)) {
    tryCatch({
      plots$plot_08 <- create_08_compositional_differences(
        compositional_results = results$compositional_results,
        output_dir = output_dir,
        save_plots = save_plots
      )
    }, error = function(e) {
      cat("  ! Skipped Plot 08 (Compositional Differences):", e$message, "\n")
    })
  }
  
  if (length(plots) > 0) {
    cat("\n=== Treatment Analysis Plots Complete (Part 2/4) ===\n")
    cat("Output directory:", file.path(output_dir, "02_treatment_plots"), "\n")
    cat("Number of plots created:", length(plots), "\n")
    cat("Plots created:\n")
    if ("plot_05" %in% names(plots)) cat("  05: Sample Organization by Treatment\n")
    if ("plot_06" %in% names(plots)) cat("  06: All Ecotypes Coherence by Treatment\n")
    if ("plot_07" %in% names(plots)) cat("  07: Treatment Effect Magnitude\n")
    if ("plot_08" %in% names(plots)) cat("  08: Compositional Differences\n")
  } else {
    cat("No treatment plots were created.\n")
  }
  
  return(plots)
}

cat("SpatialCoherence Treatment Analysis Plotting Suite (Part 2/4) loaded successfully.\n")
cat("Use create_all_treatment_plots() to generate all treatment analysis plots.\n")