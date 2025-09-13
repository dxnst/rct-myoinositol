#!/usr/bin/env python3
"""
Python alternative for basic RCT data analysis
Provides data validation and basic statistical summaries
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys

def load_and_validate_data(file_path):
    """Load CSV data and validate structure"""
    try:
        data = pd.read_csv(file_path)
        print(f"Data loaded successfully from {file_path}")
        print(f"Data shape: {data.shape}")
        
        # Expected columns
        expected_cols = [
            'nombre_estudio', 'grupo', 'formulacion', 'parametro_metabolico',
            'n', 'media', 'desviacion_estandar', 'edad_media', 'obesidad'
        ]
        
        # Check if all expected columns are present
        missing_cols = [col for col in expected_cols if col not in data.columns]
        if missing_cols:
            print(f"Warning: Missing columns: {missing_cols}")
        
        print("\nData summary:")
        print(data.describe())
        
        return data
        
    except Exception as e:
        print(f"Error loading data: {e}")
        return None

def summarize_studies(data):
    """Create summary statistics by study and parameter"""
    
    print("\n=== STUDY SUMMARY ===")
    
    # Studies overview
    studies = data['nombre_estudio'].unique()
    print(f"Number of studies: {len(studies)}")
    print(f"Studies: {', '.join(studies)}")
    
    # Parameters overview
    parameters = data['parametro_metabolico'].unique()
    print(f"\nNumber of parameters: {len(parameters)}")
    print(f"Parameters: {', '.join(parameters)}")
    
    # Group distribution
    print(f"\nGroup distribution:")
    print(data['grupo'].value_counts())
    
    # Obesity distribution
    print(f"\nObesity distribution:")
    print(data['obesidad'].value_counts())
    
    # Sample sizes
    print(f"\nSample size statistics:")
    print(f"Total participants: {data['n'].sum()}")
    print(f"Mean sample size per group: {data['n'].mean():.1f}")
    print(f"Sample size range: {data['n'].min()} - {data['n'].max()}")

def create_basic_plots(data, output_dir="output"):
    """Create basic visualization plots"""
    
    # Create output directory
    Path(output_dir).mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('seaborn-v0_8')
    sns.set_palette("husl")
    
    # Plot 1: Sample sizes by study and group
    plt.figure(figsize=(12, 6))
    pivot_n = data.pivot_table(values='n', index='nombre_estudio', 
                               columns='grupo', aggfunc='sum', fill_value=0)
    pivot_n.plot(kind='bar', ax=plt.gca())
    plt.title('Sample Sizes by Study and Group')
    plt.ylabel('Number of Participants')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/sample_sizes_by_study.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 2: Mean values by parameter and group
    fig, axes = plt.subplots(2, 2, figsize=(15, 10))
    axes = axes.ravel()
    
    parameters = data['parametro_metabolico'].unique()[:4]  # First 4 parameters
    
    for i, param in enumerate(parameters):
        param_data = data[data['parametro_metabolico'] == param]
        
        pivot_mean = param_data.pivot_table(values='media', index='nombre_estudio', 
                                            columns='grupo', aggfunc='mean')
        
        if not pivot_mean.empty:
            pivot_mean.plot(kind='bar', ax=axes[i])
            axes[i].set_title(f'Mean Values - {param}')
            axes[i].set_ylabel('Mean Value')
            axes[i].tick_params(axis='x', rotation=45)
    
    plt.tight_layout()
    plt.savefig(f"{output_dir}/mean_values_by_parameter.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # Plot 3: Effect sizes (simple difference)
    effect_sizes = []
    
    for study in data['nombre_estudio'].unique():
        for param in data['parametro_metabolico'].unique():
            study_param_data = data[
                (data['nombre_estudio'] == study) & 
                (data['parametro_metabolico'] == param)
            ]
            
            if len(study_param_data) == 2:  # Has both control and intervention
                control = study_param_data[study_param_data['grupo'] == 'control']
                intervention = study_param_data[study_param_data['grupo'] == 'intervencion']
                
                if len(control) == 1 and len(intervention) == 1:
                    effect_size = intervention['media'].iloc[0] - control['media'].iloc[0]
                    effect_sizes.append({
                        'study': study,
                        'parameter': param,
                        'effect_size': effect_size,
                        'obesity': control['obesidad'].iloc[0]
                    })
    
    if effect_sizes:
        effect_df = pd.DataFrame(effect_sizes)
        
        plt.figure(figsize=(12, 8))
        for param in effect_df['parameter'].unique():
            param_effects = effect_df[effect_df['parameter'] == param]
            plt.scatter(param_effects['parameter'], param_effects['effect_size'], 
                       label=param, s=100, alpha=0.7)
        
        plt.axhline(y=0, color='red', linestyle='--', alpha=0.5)
        plt.title('Effect Sizes by Parameter\n(Intervention - Control)')
        plt.ylabel('Effect Size (Mean Difference)')
        plt.xlabel('Metabolic Parameter')
        plt.xticks(rotation=45)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(f"{output_dir}/effect_sizes.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        # Save effect sizes to CSV
        effect_df.to_csv(f"{output_dir}/effect_sizes.csv", index=False)
    
    print(f"\nPlots saved to '{output_dir}' directory")

def main():
    """Main analysis function"""
    
    # Default data file
    data_file = "data/example_data.csv"
    
    # Check if custom data file provided
    if len(sys.argv) > 1:
        data_file = sys.argv[1]
    
    print("=== RCT DATA ANALYSIS (Python Version) ===")
    print(f"Analyzing data from: {data_file}")
    
    # Load and validate data
    data = load_and_validate_data(data_file)
    
    if data is None:
        print("Failed to load data. Exiting.")
        return
    
    # Summarize studies
    summarize_studies(data)
    
    # Create plots
    create_basic_plots(data)
    
    print("\n=== ANALYSIS COMPLETE ===")
    print("Note: For comprehensive meta-analysis and forest plots, use the R scripts.")
    print("This Python script provides basic summaries and visualizations.")

if __name__ == "__main__":
    main()