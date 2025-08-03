#!/usr/bin/env Rscript

# ============================================================================
# SpatialCoherence Plotting Suite - Part 4: Heatmap Visualization Plots
# ============================================================================
# Author: Ateeq Khaliq (akhaliq@iu.edu)
# Institution: Indiana University
# Version: 1.0.0
# 
# This script generates heatmap visualization plots
# Part 4 of 4: Complex heatmaps and organization clustering plots

# Required libraries
required_packages <- c("ggplot2", "dplyr", "reshape2", "viridis", "RColorBrewer", "ComplexHeatmap", "circlize")

# Function to check and install packages
check_and_install <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      cat("Installing package:", pkg, "\n")
      if (pkg %in% c("ComplexHeatmap", "circlize")) {
        if (!requireNamespace("BiocManager", quietly = TRUE))
          install.packages("BiocManager")
        BiocManager::install(pkg, quiet = TRUE)
      } else {
        install.packages(pkg, quiet = TRUE)
      }
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

# Create output directory structure with user-defined path
create_output_structure <- function(base_dir) {
  # Ensure base directory is created first
  if (!dir.exists(base_dir)) {
    dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  dirs <- c(
    file.path(base_dir, "04_heatmap_plots"),
    file.path(base_dir, "04_heatmap_plots", "pdf"),
    file.path(base_dir, "04_heatmap_plots", "png")
  )
  
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    }
  }
  
  return(dirs[1])
}

# Function to order samples within groups by similarity
order_within_group <- function(matrix, samples) {
  if (length(samples) <= 1) return(samples)
  group_matrix <- matrix[, samples, drop = FALSE]
  if (ncol(group_matrix) > 1) {
    dist_matrix <- dist(t(group_matrix))
    hc <- hclust(dist_matrix, method = "complete")
    return(samples[hc$order])
  }
  return(samples)
}

# ============================================================================
# PLOT 13: ORGANIZATION HEATMAP WITH COMPLEX HEATMAP (MATCHING YOUR IMAGE 8)
# ============================================================================

create_13_organization_heatmap_complex <- function(coherence_matrix = NULL,
                                                  detailed_results = NULL,
                                                  coherence_threshold = 0.47,
                                                  output_dir = ".",
                                                  save_plots = TRUE,
                                                  cluster_samples = TRUE,
                                                  cluster_ecotypes = FALSE) {
  
  cat("Creating Plot 13: Spatial Organization Heatmap (ComplexHeatmap)...\n")
  
  # Prepare matrix data
  if (!is.null(coherence_matrix)) {
    if (is.data.frame(coherence_matrix) && "metaprogram" %in% colnames(coherence_matrix)) {
      # Enhanced analysis format
      matrix_data <- coherence_matrix[, !names(coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- coherence_matrix$metaprogram
    } else {
      matrix_data <- coherence_matrix
    }
  } else if (!is.null(detailed_results)) {
    # Convert detailed results to matrix format
    if ("metaprogram" %in% colnames(detailed_results)) {
      matrix_data <- reshape2::dcast(detailed_results, metaprogram ~ sample, 
                                   value.var = "coherence_score")
      rownames(matrix_data) <- matrix_data$metaprogram
      matrix_data <- matrix_data[, -1]
    } else {
      stop("Unable to process detailed_results format")
    }
  } else {
    stop("Either coherence_matrix or detailed_results must be provided")
  }
  
  # Convert to numeric matrix and handle NAs
  coherence_final <- as.matrix(matrix_data)
  
  # Handle NAs by imputation (use row means)
  for (i in 1:nrow(coherence_final)) {
    row_mean <- mean(coherence_final[i, ], na.rm = TRUE)
    if (!is.na(row_mean)) {
      coherence_final[i, is.na(coherence_final[i, ])] <- row_mean
    }
  }
  
  # Remove problematic rows/columns (optional filtering)
  coherence_final <- coherence_final[rowSums(!is.na(coherence_final)) > ncol(coherence_final) * 0.5, ]
  coherence_final <- coherence_final[, colSums(!is.na(coherence_final)) > nrow(coherence_final) * 0.1]
  
  # Clean up sample names (remove trailing 'b' if present)
  clean_sample_names <- gsub("b$", "", colnames(coherence_final))
  colnames(coherence_final) <- clean_sample_names
  
  cat("Heatmap dimensions:", nrow(coherence_final), "ecotypes x", ncol(coherence_final), "samples\n")
  
  # Calculate sample organization status
  sample_means <- colMeans(coherence_final, na.rm = TRUE)
  sample_organization <- ifelse(sample_means > coherence_threshold + 0.08, "Organized", 
                               ifelse(sample_means > coherence_threshold - 0.08, "Intermediate", "Disorganized"))
  
  # Calculate ecotype organization status
  ecotype_means <- rowMeans(coherence_final, na.rm = TRUE)
  ecotype_org <- ifelse(ecotype_means > coherence_threshold, "Organized", "Disorganized")
  
  # Create ordered samples by organization status
  sample_order_df <- data.frame(
    sample = names(sample_means),
    organization = sample_organization,
    mean_coherence = sample_means,
    stringsAsFactors = FALSE
  )
  
  # Order samples: Organized -> Intermediate -> Disorganized, then by coherence within each group
  sample_order_df <- sample_order_df %>%
    arrange(desc(organization), desc(mean_coherence))
  
  # Order samples within each group by similarity
  organized_samples <- sample_order_df$sample[sample_order_df$organization == "Organized"]
  intermediate_samples <- sample_order_df$sample[sample_order_df$organization == "Intermediate"] 
  disorganized_samples <- sample_order_df$sample[sample_order_df$organization == "Disorganized"]
  
  if (cluster_samples) {
    organized_ordered <- order_within_group(coherence_final, organized_samples)
    intermediate_ordered <- order_within_group(coherence_final, intermediate_samples)
    disorganized_ordered <- order_within_group(coherence_final, disorganized_samples)
  } else {
    organized_ordered <- organized_samples
    intermediate_ordered <- intermediate_samples
    disorganized_ordered <- disorganized_samples
  }
  
  # Combine ordered samples
  final_sample_order <- c(organized_ordered, intermediate_ordered, disorganized_ordered)
  coherence_final_ordered <- coherence_final[, final_sample_order]
  
  # Update sample annotation for final order
  sample_org_final <- sample_organization[final_sample_order]
  names(sample_org_final) <- final_sample_order
  
  # Color scheme for heatmap
  col_fun <- circlize::colorRamp2(
    c(0, 0.25, 0.5, 0.75, 1),
    c("#D73027", "#FC8D59", "#FEE08B", "#D9EF8B", "#66BD63")
  )
  
  # Row annotation (ecotype organization)
  row_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Status = ecotype_org,
    col = list(Status = c("Organized" = PLOT_CONFIG$organized_color, 
                         "Disorganized" = PLOT_CONFIG$disorganized_color)),
    which = "row",
    annotation_name_gp = grid::gpar(fontsize = 10),
    simple_anno_size = grid::unit(0.5, "cm")
  )
  
  # Top annotation (sample organization)
  top_annotation <- ComplexHeatmap::HeatmapAnnotation(
    Sample_Status = sample_org_final,
    col = list(Sample_Status = c("Organized" = PLOT_CONFIG$organized_color, 
                                "Intermediate" = PLOT_CONFIG$intermediate_color, 
                                "Disorganized" = PLOT_CONFIG$disorganized_color)),
    annotation_name_gp = grid::gpar(fontsize = 10),
    simple_anno_size = grid::unit(0.4, "cm")
  )
  
  # Create main heatmap
  ht_final <- ComplexHeatmap::Heatmap(
    coherence_final_ordered,
    name = "Spatial\nCoherence",
    
    # Colors
    col = col_fun,
    
    # Clustering options
    cluster_rows = cluster_ecotypes,
    cluster_columns = FALSE,  # Keep organization grouping
    
    # Show names
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = grid::gpar(fontsize = 12),
    column_names_gp = grid::gpar(fontsize = 8, angle = 90),
    
    # Titles
    row_title = "Ecotypes",
    row_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    column_title = "Samples (Organized → Intermediate → Disorganized)",
    column_title_gp = grid::gpar(fontsize = 14, fontface = "bold"),
    
    # Annotations
    left_annotation = row_annotation,
    top_annotation = top_annotation,
    
    # Legend
    heatmap_legend_param = list(
      title_position = "topcenter",
      legend_direction = "vertical",
      title_gp = grid::gpar(fontsize = 12, fontface = "bold")
    ),
    
    # Cell borders
    rect_gp = grid::gpar(col = "white", lwd = 0.1),
    
    # Add column split lines to separate groups visually
    column_split = sample_org_final,
    column_gap = grid::unit(2, "mm")
  )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    
    # Save PDF version
    pdf(file.path(output_subdir, "pdf", "13_pdac_organization_heatmap_complex.pdf"), 
        width = 20, height = 10)
    ComplexHeatmap::draw(ht_final, heatmap_legend_side = "right")
    
    grid::grid.text("Spatial Organization Heatmap - Samples by Organization Status", 
                   x = 0.5, y = 0.96, 
                   gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text("Red = Disorganized, Green = Organized", 
                   x = 0.85, y = 0.02, 
                   gp = grid::gpar(fontsize = 10))
    
    dev.off()
    
    # Save PNG version (simplified)
    png(file.path(output_subdir, "png", "13_pdac_organization_heatmap_complex.png"), 
        width = 20, height = 10, units = "in", res = PLOT_CONFIG$dpi)
    ComplexHeatmap::draw(ht_final, heatmap_legend_side = "right")
    
    grid::grid.text("Spatial Organization Heatmap - Samples by Organization Status", 
                   x = 0.5, y = 0.96, 
                   gp = grid::gpar(fontsize = 16, fontface = "bold"))
    grid::grid.text("Red = Disorganized, Green = Organized", 
                   x = 0.85, y = 0.02, 
                   gp = grid::gpar(fontsize = 10))
    
    dev.off()
    
    cat("  ✓ Saved: 13_organization_heatmap_complex in pdf, png format(s)\n")
  }
  
  # Print summary statistics
  cat("\nHeatmap Summary:\n")
  org_summary <- sample_order_df %>%
    group_by(organization) %>%
    summarise(
      count = n(),
      mean_coherence = round(mean(mean_coherence), 3),
      range = paste(round(min(mean_coherence), 3), "-", round(max(mean_coherence), 3)),
      .groups = 'drop'
    )
  print(org_summary)
  
  return(ht_final)
}

# ============================================================================
# PLOT 14: SIMPLE COHERENCE HEATMAP WITH GGPLOT (ALTERNATIVE VERSION)
# ============================================================================

create_14_coherence_heatmap_ggplot <- function(coherence_matrix = NULL,
                                             detailed_results = NULL,
                                             output_dir = ".",
                                             save_plots = TRUE,
                                             cluster_samples = TRUE,
                                             cluster_ecotypes = TRUE) {
  
  cat("Creating Plot 14: Coherence Scores Heatmap (ggplot version)...\n")
  
  # Prepare matrix data (same as complex heatmap)
  if (!is.null(coherence_matrix)) {
    if (is.data.frame(coherence_matrix) && "metaprogram" %in% colnames(coherence_matrix)) {
      matrix_data <- coherence_matrix[, !names(coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- coherence_matrix$metaprogram
    } else {
      matrix_data <- coherence_matrix
    }
  } else if (!is.null(detailed_results)) {
    if ("metaprogram" %in% colnames(detailed_results)) {
      matrix_data <- reshape2::dcast(detailed_results, metaprogram ~ sample, 
                                   value.var = "coherence_score")
      rownames(matrix_data) <- matrix_data$metaprogram
      matrix_data <- matrix_data[, -1]
    } else {
      stop("Unable to process detailed_results format")
    }
  } else {
    stop("Either coherence_matrix or detailed_results must be provided")
  }
  
  # Convert to numeric matrix
  matrix_data <- as.matrix(matrix_data)
  matrix_data[is.na(matrix_data)] <- 0
  
  # Prepare for ggplot
  plot_data <- reshape2::melt(matrix_data)
  colnames(plot_data) <- c("Ecotype", "Sample", "Coherence")
  
  # Remove zero values for cleaner visualization
  plot_data <- plot_data[plot_data$Coherence > 0, ]
  
  # Optional clustering
  if (cluster_samples && ncol(matrix_data) > 1) {
    sample_order <- colnames(matrix_data)[hclust(dist(t(matrix_data)))$order]
    plot_data$Sample <- factor(plot_data$Sample, levels = sample_order)
  }
  
  if (cluster_ecotypes && nrow(matrix_data) > 1) {
    ecotype_order <- rownames(matrix_data)[hclust(dist(matrix_data))$order]
    plot_data$Ecotype <- factor(plot_data$Ecotype, levels = ecotype_order)
  }
  
  # Create ggplot heatmap
  p14 <- ggplot(plot_data, aes(x = Sample, y = Ecotype, fill = Coherence)) +
    geom_tile(color = "white", size = 0.1) +
    scale_fill_viridis_c(name = "Spatial\nCoherence", option = "plasma") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      panel.grid = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "Ecotype",
      title = "Spatial Coherence Heatmap (ggplot version)",
      subtitle = paste0("Samples: ", length(unique(plot_data$Sample)), 
                       " | Ecotypes: ", length(unique(plot_data$Ecotype)))
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    # Use dynamic dimensions based on data size
    plot_width <- max(12, length(unique(plot_data$Sample)) * 0.3)
    plot_height <- max(8, length(unique(plot_data$Ecotype)) * 0.5)
    
    ggsave(file.path(output_subdir, "pdf", "14_coherence_heatmap_ggplot.pdf"), 
           p14, width = plot_width, height = plot_height, device = "pdf", dpi = PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "png", "14_coherence_heatmap_ggplot.png"), 
           p14, width = plot_width, height = plot_height, device = "png", dpi = PLOT_CONFIG$dpi)
    cat("  ✓ Saved: 14_coherence_heatmap_ggplot in pdf, png format(s)\n")
  }
  
  return(p14)
}

# ============================================================================
# PLOT 15: SAMPLE ORGANIZATION SUMMARY TABLE VISUALIZATION
# ============================================================================

create_15_sample_organization_summary <- function(coherence_matrix = NULL,
                                                 detailed_results = NULL,
                                                 coherence_threshold = 0.47,
                                                 output_dir = ".",
                                                 save_plots = TRUE) {
  
  cat("Creating Plot 15: Sample Organization Summary Table...\n")
  
  # Prepare data (similar to heatmap preparation)
  if (!is.null(coherence_matrix)) {
    if (is.data.frame(coherence_matrix) && "metaprogram" %in% colnames(coherence_matrix)) {
      matrix_data <- coherence_matrix[, !names(coherence_matrix) %in% "metaprogram"]
      rownames(matrix_data) <- coherence_matrix$metaprogram
    } else {
      matrix_data <- coherence_matrix
    }
  } else if (!is.null(detailed_results)) {
    if ("metaprogram" %in% colnames(detailed_results)) {
      matrix_data <- reshape2::dcast(detailed_results, metaprogram ~ sample, 
                                   value.var = "coherence_score")
      rownames(matrix_data) <- matrix_data$metaprogram
      matrix_data <- matrix_data[, -1]
    } else {
      stop("Unable to process detailed_results format")
    }
  } else {
    stop("Either coherence_matrix or detailed_results must be provided")
  }
  
  # Calculate sample summary statistics
  sample_means <- colMeans(matrix_data, na.rm = TRUE)
  sample_organization <- ifelse(sample_means > coherence_threshold + 0.08, "Organized", 
                               ifelse(sample_means > coherence_threshold - 0.08, "Intermediate", "Disorganized"))
  
  summary_table <- data.frame(
    sample = names(sample_means),
    mean_coherence = round(sample_means, 3),
    organization = sample_organization,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(organization), desc(mean_coherence))
  
  # Create visualization of the summary table
  summary_table$sample <- factor(summary_table$sample, levels = summary_table$sample)
  
  p15 <- ggplot(summary_table, aes(x = sample, y = mean_coherence, fill = organization)) +
    geom_col(alpha = 0.8, color = "black", size = 0.3) +
    geom_hline(yintercept = coherence_threshold, linetype = "dashed", color = "black", size = 0.8) +
    scale_fill_manual(
      values = c("Organized" = PLOT_CONFIG$organized_color, 
                 "Intermediate" = PLOT_CONFIG$intermediate_color, 
                 "Disorganized" = PLOT_CONFIG$disorganized_color),
      name = "Organization\nStatus"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      panel.grid.major.x = element_blank()
    ) +
    labs(
      x = "Sample",
      y = "Mean Coherence Score",
      title = "Sample Organization Summary",
      subtitle = paste0("Threshold = ", coherence_threshold, " | Total Samples: ", nrow(summary_table))
    )
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    
    # Save the plot
    ggsave(file.path(output_subdir, "pdf", "15_sample_organization_summary.pdf"), 
           p15, width = max(12, nrow(summary_table) * 0.2), height = 8, 
           device = "pdf", dpi = PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "png", "15_sample_organization_summary.png"), 
           p15, width = max(12, nrow(summary_table) * 0.2), height = 8, 
           device = "png", dpi = PLOT_CONFIG$dpi)
    
    # Save the summary table as CSV
    write.csv(summary_table, file.path(output_subdir, "sample_organization_summary.csv"), row.names = FALSE)
    
    cat("  ✓ Saved: 15_sample_organization_summary in pdf, png format(s)\n")
    cat("  ✓ Saved: sample_organization_summary.csv\n")
  }
  
  return(list(plot = p15, summary_table = summary_table))
}

# ============================================================================
# PLOT 16: ECOTYPE ORGANIZATION COMPARISON
# ============================================================================

create_16_ecotype_organization_comparison <- function(mean_coherence_data,
                                                    organization_results,
                                                    coherence_threshold = 0.47,
                                                    output_dir = ".",
                                                    save_plots = TRUE,
                                                    ecotype_annotations = NULL) {
  
  cat("Creating Plot 16: Ecotype Organization Comparison...\n")
  
  # Prepare data
  plot_data <- data.frame(
    ecotype = names(mean_coherence_data),
    mean_coherence = as.numeric(mean_coherence_data),
    organization = organization_results[names(mean_coherence_data)],
    stringsAsFactors = FALSE
  )
  
  # Remove NA values
  plot_data <- plot_data[!is.na(plot_data$mean_coherence), ]
  
  # Add annotations if provided
  if (!is.null(ecotype_annotations)) {
    plot_data$annotation <- ecotype_annotations[plot_data$ecotype]
    plot_data$display_name <- ifelse(is.na(plot_data$annotation), 
                                   plot_data$ecotype, 
                                   paste0(plot_data$ecotype, "\n(", plot_data$annotation, ")"))
  } else {
    plot_data$display_name <- plot_data$ecotype
  }
  
  # Create comparison plot showing organized vs disorganized
  p16 <- ggplot(plot_data, aes(x = reorder(display_name, mean_coherence), y = mean_coherence)) +
    geom_segment(aes(x = display_name, xend = display_name, 
                    y = coherence_threshold, yend = mean_coherence, 
                    color = organization), 
                size = 3, alpha = 0.8) +
    geom_point(aes(color = organization), size = 4) +
    geom_hline(yintercept = coherence_threshold, linetype = "dashed", color = "black", size = 0.8) +
    scale_color_manual(
      values = c("Organized" = PLOT_CONFIG$organized_color, 
                 "Disorganized" = PLOT_CONFIG$disorganized_color,
                 "Intermediate" = PLOT_CONFIG$intermediate_color),
      name = "Organization\nClassification"
    ) +
    coord_flip() +
    theme_minimal() +
    theme(
      axis.text = element_text(size = PLOT_CONFIG$axis_text_size),
      axis.title = element_text(size = PLOT_CONFIG$axis_title_size, face = "bold"),
      plot.title = element_text(size = PLOT_CONFIG$title_size, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = PLOT_CONFIG$subtitle_size, hjust = 0.5),
      legend.title = element_text(size = PLOT_CONFIG$legend_text_size, face = "bold"),
      legend.text = element_text(size = PLOT_CONFIG$legend_text_size),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Ecotype",
      y = "Mean Coherence Score",
      title = "Ecotype Organization Comparison",
      subtitle = paste0("Distance from threshold (", coherence_threshold, ") indicates organization strength")
    ) +
    annotate("text", x = 1, y = coherence_threshold + 0.05, 
             label = paste("Threshold =", coherence_threshold), 
             size = 3, hjust = 0)
  
  if (save_plots) {
    output_subdir <- create_output_structure(output_dir)
    
    ggsave(file.path(output_subdir, "pdf", "16_ecotype_organization_comparison.pdf"), 
           p16, width = 10, height = max(8, nrow(plot_data) * 0.4), 
           device = "pdf", dpi = PLOT_CONFIG$dpi)
    ggsave(file.path(output_subdir, "png", "16_ecotype_organization_comparison.png"), 
           p16, width = 10, height = max(8, nrow(plot_data) * 0.4), 
           device = "png", dpi = PLOT_CONFIG$dpi)
    
    cat("  ✓ Saved: 16_ecotype_organization_comparison in pdf, png format(s)\n")
  }
  
  return(p16)
}

# ============================================================================
# MAIN FUNCTION: CREATE ALL HEATMAP PLOTS
# ============================================================================

create_all_heatmap_plots <- function(results, 
                                    output_dir = "spatial_coherence_plots",
                                    save_plots = TRUE,
                                    coherence_threshold = NULL,
                                    cluster_samples = TRUE,
                                    cluster_ecotypes = FALSE) {
  
  cat("=== SpatialCoherence Heatmap Visualization Plots (Part 4/4) ===\n\n")
  
  # Extract threshold from results if not provided
  if (is.null(coherence_threshold)) {
    if (!is.null(results$analysis_parameters$coherence_threshold)) {
      coherence_threshold <- results$analysis_parameters$coherence_threshold
    } else {
      coherence_threshold <- 0.47  # Default
    }
  }
  
  plots <- list()
  
  # Check if coherence matrix or detailed results are available
  if (is.null(results$coherence_matrix) && is.null(results$detailed_results)) {
    cat("No coherence matrix or detailed results available. Skipping heatmap plots.\n")
    return(plots)
  }
  
  # Plot 13: Complex Organization Heatmap
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) && requireNamespace("circlize", quietly = TRUE)) {
    plots$plot_13 <- create_13_organization_heatmap_complex(
      coherence_matrix = results$coherence_matrix,
      detailed_results = results$detailed_results,
      coherence_threshold = coherence_threshold,
      output_dir = output_dir,
      save_plots = save_plots,
      cluster_samples = cluster_samples,
      cluster_ecotypes = cluster_ecotypes
    )
  } else {
    cat("ComplexHeatmap package not available. Skipping complex heatmap.\n")
  }
  
  # Plot 14: ggplot Heatmap (alternative)
  plots$plot_14 <- create_14_coherence_heatmap_ggplot(
    coherence_matrix = results$coherence_matrix,
    detailed_results = results$detailed_results,
    output_dir = output_dir,
    save_plots = save_plots,
    cluster_samples = cluster_samples,
    cluster_ecotypes = cluster_ecotypes
  )
  
  # Plot 15: Sample Organization Summary
  plots$plot_15 <- create_15_sample_organization_summary(
    coherence_matrix = results$coherence_matrix,
    detailed_results = results$detailed_results,
    coherence_threshold = coherence_threshold,
    output_dir = output_dir,
    save_plots = save_plots
  )
  
  # Plot 16: Ecotype Organization Comparison
  plots$plot_16 <- create_16_ecotype_organization_comparison(
    mean_coherence_data = results$mean_coherence,
    organization_results = results$organization_results,
    coherence_threshold = coherence_threshold,
    output_dir = output_dir,
    save_plots = save_plots,
    ecotype_annotations = results$ecotype_annotations
  )
  
  cat("\n=== Heatmap Visualization Plots Complete (Part 4/4) ===\n")
  cat("Output directory:", file.path(output_dir, "04_heatmap_plots"), "\n")
  cat("Number of plots created:", length(plots), "\n")
  cat("Plots created:\n")
  if ("plot_13" %in% names(plots)) cat("  13: Spatial Organization Heatmap (ComplexHeatmap)\n")
  if ("plot_14" %in% names(plots)) cat("  14: Coherence Heatmap (ggplot version)\n")
  if ("plot_15" %in% names(plots)) cat("  15: Sample Organization Summary\n")
  if ("plot_16" %in% names(plots)) cat("  16: Ecotype Organization Comparison\n")
  
  return(plots)
}

# ============================================================================
# MASTER FUNCTION: CREATE ALL SPATIALCOHERENCE PLOTS
# ============================================================================

create_all_spatialcoherence_plots <- function(results, 
                                             output_dir = "spatial_coherence_plots",
                                             save_plots = TRUE,
                                             treatment_column = "treatment") {
  
  cat("============================================================================\n")
  cat("SpatialCoherence Complete Plotting Suite - All 4 Parts\n")
  cat("============================================================================\n\n")
  
  all_plots <- list()
  
  # Part 1: Basic Plots
  if (exists("create_all_basic_plots")) {
    basic_plots <- create_all_basic_plots(results, output_dir, save_plots)
    all_plots <- c(all_plots, basic_plots)
  }
  
  # Part 2: Treatment Plots
  if (exists("create_all_treatment_plots")) {
    treatment_plots <- create_all_treatment_plots(results, output_dir, save_plots, treatment_column)
    all_plots <- c(all_plots, treatment_plots)
  }
  
  # Part 3: Segmentation Plots
  if (exists("create_all_segmentation_plots")) {
    segmentation_plots <- create_all_segmentation_plots(results, output_dir, save_plots, treatment_column)
    all_plots <- c(all_plots, segmentation_plots)
  }
  
  # Part 4: Heatmap Plots
  heatmap_plots <- create_all_heatmap_plots(results, output_dir, save_plots)
  all_plots <- c(all_plots, heatmap_plots)
  
  cat("\n============================================================================\n")
  cat("SpatialCoherence Complete Plotting Suite - FINISHED\n")
  cat("============================================================================\n")
  cat("Total plots created:", length(all_plots), "\n")
  cat("Output directory:", output_dir, "\n")
  cat("All plots have been saved in organized subdirectories.\n")
  
  return(all_plots)
}

cat("SpatialCoherence Heatmap Visualization Plotting Suite (Part 4/4) loaded successfully.\n")
cat("Use create_all_heatmap_plots() to generate all heatmap visualization plots.\n")
cat("Use create_all_spatialcoherence_plots() to generate ALL plots from all 4 parts.\n")