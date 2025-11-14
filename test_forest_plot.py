#!/usr/bin/env python3
"""
Test script for the enhanced forest plot functionality
"""

import pandas as pd
import numpy as np
import sys
import os

# Import from the original analysis file to get the data processing functions
from analisis_inositol import *

# Import the new forest plot functionality
from analisis_inositol_forest_plot import crear_forest_plot_combinado, crear_tabla_resumen_con_abs

def test_forest_plot_combinado():
    """
    Test the new combined forest plot functionality
    """
    print("Testing forest plot combinado functionality...")
    
    # Use the same data processing from the original script
    # Load the data
    df = pd.read_csv('data/example_data.csv')
    
    # Process data using existing functions from analisis_inositol.py
    resultados = []
    
    # Get unique parameters
    parametros_unicos = df['metabolic_parameter'].unique()
    
    # Process a few parameters for testing
    test_parametros = ['bmi_kg_m2', 'glucose_mg_dl', 'homa_ir']
    
    for parametro in test_parametros:
        if parametro not in parametros_unicos:
            continue
            
        print(f"Processing {parametro}...")
        
        # Get data for this parameter
        param_df = df[df['metabolic_parameter'] == parametro]
        
        # Process intervention groups
        interventions = param_df[param_df['group_type'] == 'intervention']
        controls = param_df[param_df['group_type'] == 'control']
        
        # Group by study for comparisons
        for study_id in param_df['study_id'].unique():
            study_int = interventions[interventions['study_id'] == study_id]
            study_ctrl = controls[controls['study_id'] == study_id]
            
            if len(study_int) == 0 or len(study_ctrl) == 0:
                continue
                
            # Take first entry for each group (simplified)
            int_row = study_int.iloc[0]
            ctrl_row = study_ctrl.iloc[0]
            
            # Calculate Hedges' g
            g, se, ci_lower, ci_upper, p_value, mean_diff = calcular_hedges_g_entre_grupos(
                int_row['post_mean'], ctrl_row['post_mean'],
                int_row['post_sd'], ctrl_row['post_sd'],
                int_row['sample_size'], ctrl_row['sample_size'],
                parametro
            )
            
            # Calculate changes
            int_change = int_row['post_mean'] - int_row['baseline_mean']
            ctrl_change = ctrl_row['post_mean'] - ctrl_row['baseline_mean']
            int_pct_change = (int_change / int_row['baseline_mean']) * 100
            ctrl_pct_change = (ctrl_change / ctrl_row['baseline_mean']) * 100
            
            # Determine treatment type
            treatment_formulation = int_row['treatment_formulation']
            tipo_tratamiento = 'MI'  # Default
            if 'myo_inositol' in treatment_formulation and 'd_chiro_inositol' in treatment_formulation:
                tipo_tratamiento = 'MI+DCI'
            elif 'd_chiro_inositol' in treatment_formulation:
                tipo_tratamiento = 'DCI'
            elif 'metformin' in treatment_formulation:
                tipo_tratamiento = 'MI+MET'
            
            # Determine control type
            control_formulation = ctrl_row['treatment_formulation']
            tipo_control = 'Control'
            if 'metformin' in control_formulation:
                tipo_control = 'Metformina'
            elif 'diet' in control_formulation:
                tipo_control = 'Dieta'
            
            result_data = {
                'estudio': int_row['study_id'],
                'pais': int_row['country'],
                'parametro': parametro,
                'parametro_es': PARAMETROS.get(parametro, parametro),
                'tratamiento': treatment_formulation,
                'tipo_tratamiento': tipo_tratamiento,
                'tipo_control': tipo_control,
                'duracion_semanas': int_row['duration_weeks'],
                'obesidad': 'Sí' if int_row['has_obesity'] == 'yes' else 'No',
                'n': int_row['sample_size'],
                'n_control': ctrl_row['sample_size'],
                'valor_inicial': int_row['baseline_mean'],
                'valor_final': int_row['post_mean'],
                'cambio_absoluto': int_change,
                'cambio_porcentual': int_pct_change,
                'valor_inicial_control': ctrl_row['baseline_mean'],
                'valor_final_control': ctrl_row['post_mean'],
                'cambio_absoluto_control': ctrl_change,
                'cambio_porcentual_control': ctrl_pct_change,
                'efecto': g,
                'error_estandar': se,
                'ci_inferior': ci_lower,
                'ci_superior': ci_upper,
                'p_valor': p_value,
                'vs_control': True
            }
            
            resultados.append(result_data)
    
    # Convert to DataFrame
    resultados_df = pd.DataFrame(resultados)
    
    if len(resultados_df) == 0:
        print("No data processed for testing")
        return False
    
    # Test the new forest plot function
    print("Creating test forest plots...")
    
    for parametro in test_parametros:
        param_data = resultados_df[resultados_df['parametro'] == parametro]
        
        if len(param_data) < 2:
            print(f"Insufficient data for {parametro}")
            continue
            
        # Check if we have both obesity groups
        has_normal = len(param_data[param_data['obesidad'] == 'No']) > 0
        has_obese = len(param_data[param_data['obesidad'] == 'Sí']) > 0
        
        if not (has_normal or has_obese):
            print(f"No obesity group data for {parametro}")
            continue
        
        print(f"Creating combined forest plot for {parametro}...")
        
        # Create the combined forest plot
        fig = crear_forest_plot_combinado(
            resultados_df, 
            parametro, 
            "Test Forest Plot Combinado",
            comparacion=True,
            filename=f"test_output/test_forest_plot_combinado_{parametro}.pdf"
        )
        
        if fig is not None:
            print(f"✓ Successfully created forest plot for {parametro}")
        else:
            print(f"✗ Failed to create forest plot for {parametro}")
    
    # Test the summary table function
    print("Testing summary table with absolute values...")
    
    for parametro in test_parametros:
        param_data = resultados_df[resultados_df['parametro'] == parametro]
        
        if len(param_data) < 3:
            continue
            
        summary_table = crear_tabla_resumen_con_abs(resultados_df, parametro, vs_control=True)
        
        if summary_table is not None:
            print(f"✓ Successfully created summary table for {parametro}")
            print(f"  Table shape: {summary_table.shape}")
            
            # Save to CSV
            summary_table.to_csv(f"test_output/test_tabla_{parametro}_abs.csv", index=False)
            print(f"  Saved to test_output/test_tabla_{parametro}_abs.csv")
        else:
            print(f"✗ Failed to create summary table for {parametro}")
    
    print("Test completed!")
    return True

if __name__ == "__main__":
    # Create test directory
    os.makedirs('test_output', exist_ok=True)
    
    # Run tests
    success = test_forest_plot_combinado()
    
    if success:
        print("\n✓ All tests passed!")
        sys.exit(0)
    else:
        print("\n✗ Some tests failed!")
        sys.exit(1)