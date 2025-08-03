# Step 1: Load the data from your CSV files
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/00_load_spatialcoherence_data.R")
results <- load_spatialcoherence_data("/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test")

# Step 2: Load plotting functions
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/01_basic_coherence_plots.R")
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/02_treatment_analysis_plots.R")
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/03_segmentation_analysis_plots.R")
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/04_heatmap_visualization_plots.R")
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/05_master_plotting_script.r")

# Step 3: Create all plots
all_plots <- create_all_spatialcoherence_plots(
  results = results,
  output_dir = "/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test/plotting",
  save_plots = TRUE
)


# Generate all plots including missing ones
comprehensive_results <- create_comprehensive_spatial_analysis(
  data_dir = "/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test",
  output_dir = "/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test/plotting/missing_plots",
  treatment_column = "treatment"
)

#comprehensive script to quantify and rank samples from most organized to most disorganized based on coherence scores. 

# Load the script
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/06_sample_organization_quantifier.r")

# Load your data (using your existing loader)
source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/00_load_spatialcoherence_data.R")
results <- load_spatialcoherence_data("/Users/akhaliq/Desktop/pdac_nac/temp_coherence_results/pdac_spatial_coherence_test")

# Run comprehensive quantification
organization_analysis <- quantify_sample_organization(
  results = results,
  output_dir = "sample_organization_analysis",
  coherence_threshold = 0.47,  # or NULL to use default
  save_csvs = TRUE,
  generate_plot = TRUE,
  verbose = TRUE
)

### script that will use your organization ranking data to visualize the most organized vs most disorganized samples on H&E images. 

source("/Users/akhaliq/Downloads/SpatialCoherence-main/plotting/07_h&e_organization_plotter.r")
