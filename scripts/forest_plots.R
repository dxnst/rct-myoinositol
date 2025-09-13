# forest_plots.R
# Advanced forest plot generation for RCT meta-analysis

# Source the main analysis script
source("scripts/data_analysis.R")

# Advanced forest plot function
create_advanced_forest_plot <- function(meta_result, parameter, obesity_subgroup = FALSE) {
  
  # Extract data for custom plotting
  study_data <- meta_result$data
  
  if (obesity_subgroup && "obesity" %in% colnames(study_data)) {
    # Create subgroup analysis by obesity status
    cat("Creating subgroup analysis by obesity status\n")
    
    # Separate by obesity groups
    obese_studies <- study_data[study_data$obesity == "si", ]
    non_obese_studies <- study_data[study_data$obesity == "no", ]
    
    # Perform separate meta-analyses
    if (nrow(obese_studies) > 0) {
      meta_obese <- update(meta_result, subset = (obesity == "si"))
    }
    
    if (nrow(non_obese_studies) > 0) {
      meta_non_obese <- update(meta_result, subset = (obesity == "no"))
    }
    
    # Create combined forest plot with subgroups
    forest(meta_result,
           subgroup = obesity,
           main = paste("Forest Plot -", parameter, "(by Obesity Status)"),
           xlab = "Mean Difference",
           xlim = c(-5, 5),
           col.diamond = "blue",
           col.diamond.lines = "darkblue")
    
  } else {
    # Standard forest plot
    forest(meta_result,
           main = paste("Forest Plot -", parameter),
           xlab = "Mean Difference (95% CI)",
           xlim = c(-5, 5),
           col.diamond = "red",
           col.diamond.lines = "darkred",
           col.square = "darkgreen",
           col.square.lines = "darkgreen")
  }
}

# Function to create comprehensive forest plots for all parameters
create_all_forest_plots <- function(data_file = "data/example_data.csv") {
  
  # Load and analyze data
  data <- load_rct_data(data_file)
  parameters <- unique(data$parametro_metabolico)
  
  # Create output directory
  if (!dir.exists("output/forest_plots")) {
    dir.create("output/forest_plots", recursive = TRUE)
  }
  
  for (param in parameters) {
    cat("\nCreating forest plots for:", param, "\n")
    
    # Prepare and analyze data
    analysis_results <- analyze_parameter(data, param)
    
    if (!is.null(analysis_results)) {
      # Standard forest plot
      png(paste0("output/forest_plots/", gsub("[^A-Za-z0-9]", "_", param), "_standard.png"),
          width = 1400, height = 1000, res = 150)
      
      create_advanced_forest_plot(analysis_results$result, param, obesity_subgroup = FALSE)
      dev.off()
      
      # Subgroup forest plot by obesity
      png(paste0("output/forest_plots/", gsub("[^A-Za-z0-9]", "_", param), "_by_obesity.png"),
          width = 1400, height = 1200, res = 150)
      
      create_advanced_forest_plot(analysis_results$result, param, obesity_subgroup = TRUE)
      dev.off()
      
      cat("Forest plots created for:", param, "\n")
    }
  }
  
  cat("\nAll forest plots created in 'output/forest_plots' directory\n")
}

# Function to create summary table
create_summary_table <- function(data_file = "data/example_data.csv") {
  
  data <- load_rct_data(data_file)
  parameters <- unique(data$parametro_metabolico)
  
  summary_results <- data.frame()
  
  for (param in parameters) {
    analysis_results <- analyze_parameter(data, param)
    
    if (!is.null(analysis_results)) {
      meta_result <- analysis_results$result
      
      summary_row <- data.frame(
        Parameter = param,
        N_Studies = meta_result$k,
        N_Participants = sum(meta_result$n.e + meta_result$n.c),
        Effect_Size = round(meta_result$TE.random, 3),
        CI_Lower = round(meta_result$lower.random, 3),
        CI_Upper = round(meta_result$upper.random, 3),
        P_Value = round(meta_result$pval.random, 4),
        I2 = round(meta_result$I2, 1),
        Tau2 = round(meta_result$tau2, 4)
      )
      
      summary_results <- rbind(summary_results, summary_row)
    }
  }
  
  # Save summary table
  write_csv(summary_results, "output/meta_analysis_summary.csv")
  
  cat("\nSummary table created in 'output/meta_analysis_summary.csv'\n")
  print(summary_results)
  
  return(summary_results)
}

# Comprehensive analysis function
run_comprehensive_analysis <- function(data_file = "data/example_data.csv") {
  cat("=== COMPREHENSIVE RCT META-ANALYSIS ===\n")
  
  # Create all forest plots
  create_all_forest_plots(data_file)
  
  # Create summary table
  summary_table <- create_summary_table(data_file)
  
  cat("\n=== COMPREHENSIVE ANALYSIS COMPLETE ===\n")
  cat("Check 'output' directory for all results\n")
  
  return(summary_table)
}

# Run comprehensive analysis if script is executed directly
if (!interactive()) {
  results <- run_comprehensive_analysis()
}