# data_analysis.R
# Main script for analyzing RCT data on myo-inositol effects

# Source the setup script
source("scripts/setup_environment.R")

# Data loading and preprocessing functions
load_rct_data <- function(file_path) {
  cat("Loading data from:", file_path, "\n")
  data <- read_csv(file_path, locale = locale(encoding = "UTF-8"))
  
  # Convert character variables to factors where appropriate
  data$nombre_estudio <- as.factor(data$nombre_estudio)
  data$grupo <- as.factor(data$grupo)
  data$formulacion <- as.factor(data$formulacion)
  data$parametro_metabolico <- as.factor(data$parametro_metabolico)
  data$obesidad <- as.factor(data$obesidad)
  
  cat("Data loaded successfully. Dimensions:", dim(data), "\n")
  return(data)
}

# Function to prepare data for meta-analysis
prepare_meta_data <- function(data, parameter) {
  cat("Preparing meta-analysis data for parameter:", parameter, "\n")
  
  # Filter data for specific metabolic parameter
  param_data <- data %>%
    filter(parametro_metabolico == parameter) %>%
    arrange(nombre_estudio, grupo)
  
  # Create paired data (control vs intervention)
  studies <- unique(param_data$nombre_estudio)
  meta_data <- data.frame()
  
  for (study in studies) {
    study_data <- param_data %>% filter(nombre_estudio == study)
    
    control <- study_data %>% filter(grupo == "control")
    intervention <- study_data %>% filter(grupo == "intervencion")
    
    if (nrow(control) == 1 && nrow(intervention) == 1) {
      meta_row <- data.frame(
        study = study,
        n_control = control$n,
        mean_control = control$media,
        sd_control = control$desviacion_estandar,
        n_intervention = intervention$n,
        mean_intervention = intervention$media,
        sd_intervention = intervention$desviacion_estandar,
        age_mean = mean(c(control$edad_media, intervention$edad_media)),
        obesity = control$obesidad,
        formulation = intervention$formulacion
      )
      meta_data <- rbind(meta_data, meta_row)
    }
  }
  
  return(meta_data)
}

# Function to perform meta-analysis
perform_meta_analysis <- function(meta_data, parameter) {
  cat("Performing meta-analysis for:", parameter, "\n")
  
  # Calculate effect sizes (mean differences)
  meta_result <- metacont(
    n.e = meta_data$n_intervention,
    mean.e = meta_data$mean_intervention,
    sd.e = meta_data$sd_intervention,
    n.c = meta_data$n_control,
    mean.c = meta_data$mean_control,
    sd.c = meta_data$sd_control,
    studlab = meta_data$study,
    data = meta_data,
    sm = "MD",  # Mean Difference
    method.tau = "REML",  # Random effects model
    title = paste("Meta-analysis of", parameter)
  )
  
  return(meta_result)
}

# Main analysis function
analyze_parameter <- function(data, parameter) {
  cat("\n=== ANALYZING PARAMETER:", parameter, "===\n")
  
  # Prepare data
  meta_data <- prepare_meta_data(data, parameter)
  
  if (nrow(meta_data) == 0) {
    cat("No paired studies found for parameter:", parameter, "\n")
    return(NULL)
  }
  
  cat("Number of studies included:", nrow(meta_data), "\n")
  print(meta_data)
  
  # Perform meta-analysis
  meta_result <- perform_meta_analysis(meta_data, parameter)
  
  # Print summary
  cat("\n--- META-ANALYSIS SUMMARY ---\n")
  print(summary(meta_result))
  
  return(list(data = meta_data, result = meta_result))
}

# Export function
export_results <- function(analysis_results, parameter) {
  # Create output directory if it doesn't exist
  if (!dir.exists("output")) {
    dir.create("output")
  }
  
  # Save meta-analysis data
  write_csv(analysis_results$data, 
           file = paste0("output/", parameter, "_meta_data.csv"))
  
  # Save forest plot
  png(paste0("output/", parameter, "_forest_plot.png"), 
      width = 1200, height = 800, res = 150)
  forest(analysis_results$result, 
         main = paste("Forest Plot -", parameter),
         xlim = c(-3, 3))
  dev.off()
  
  cat("Results exported for parameter:", parameter, "\n")
}

# Main execution function
main <- function() {
  cat("=== RCT MYO-INOSITOL META-ANALYSIS ===\n")
  
  # Load data
  data <- load_rct_data("data/example_data.csv")
  
  # Get unique metabolic parameters
  parameters <- unique(data$parametro_metabolico)
  cat("Available parameters:", paste(parameters, collapse = ", "), "\n")
  
  # Analyze each parameter
  all_results <- list()
  for (param in parameters) {
    results <- analyze_parameter(data, param)
    if (!is.null(results)) {
      all_results[[param]] <- results
      export_results(results, param)
    }
  }
  
  cat("\n=== ANALYSIS COMPLETE ===\n")
  cat("Results saved in 'output' directory\n")
  
  return(all_results)
}

# Run main analysis if script is executed directly
if (!interactive()) {
  results <- main()
}