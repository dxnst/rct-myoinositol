#!/usr/bin/env python3
"""
Ejemplo de uso del forest plot combinado mejorado
Este script demuestra el uso de la nueva funcionalidad del forest plot
que muestra grupos de obesidad separados con mejor espaciado y valores absolutos.
"""

import pandas as pd
import numpy as np
import os
from matplotlib.backends.backend_pdf import PdfPages

# Import the original analysis functions
from analisis_inositol import *

# Import the new enhanced forest plot functionality
from analisis_inositol_forest_plot import crear_forest_plot_combinado, crear_tabla_resumen_con_abs

def main():
    """
    Ejemplo principal de uso del forest plot combinado
    """
    print("=== Ejemplo de Forest Plot Combinado Mejorado ===")
    print("Este ejemplo demuestra las mejoras implementadas:")
    print("1. Valores absolutos de Hedges' g para mayor claridad")
    print("2. Mejor espaciado entre grupos de obesidad")
    print("3. Preservación de signos en cambios")
    print("4. Tablas CSV actualizadas con valores absolutos")
    print()
    
    # Crear directorio de salida
    output_dir = 'ejemplo_forest_plots'
    os.makedirs(output_dir, exist_ok=True)
    
    # Cargar los datos procesados del análisis original
    print("Cargando datos...")
    
    # Usar los datos ya procesados por el script principal
    # (Para este ejemplo, reutilizamos la lógica de procesamiento)
    df = pd.read_csv('data/example_data.csv')
    
    # Procesar datos usando las funciones existentes
    resultados = []
    
    parametros_disponibles = df['metabolic_parameter'].unique()
    print(f"Parámetros disponibles: {list(parametros_disponibles)}")
    
    # Procesar cada parámetro
    for parametro in parametros_disponibles:
        param_df = df[df['metabolic_parameter'] == parametro]
        
        # Obtener grupos de intervención y control
        interventions = param_df[param_df['group_type'] == 'intervention']
        controls = param_df[param_df['group_type'] == 'control']
        
        # Procesar por estudio para comparaciones
        for study_id in param_df['study_id'].unique():
            study_int = interventions[interventions['study_id'] == study_id]
            study_ctrl = controls[controls['study_id'] == study_id]
            
            if len(study_int) == 0 or len(study_ctrl) == 0:
                continue
                
            # Tomar primera entrada de cada grupo
            int_row = study_int.iloc[0]
            ctrl_row = study_ctrl.iloc[0]
            
            # Calcular g de Hedges
            g, se, ci_lower, ci_upper, p_value, mean_diff = calcular_hedges_g_entre_grupos(
                int_row['post_mean'], ctrl_row['post_mean'],
                int_row['post_sd'], ctrl_row['post_sd'],
                int_row['sample_size'], ctrl_row['sample_size'],
                parametro
            )
            
            # Calcular cambios
            int_change = int_row['post_mean'] - int_row['baseline_mean']
            ctrl_change = ctrl_row['post_mean'] - ctrl_row['baseline_mean']
            int_pct_change = (int_change / int_row['baseline_mean']) * 100
            ctrl_pct_change = (ctrl_change / ctrl_row['baseline_mean']) * 100
            
            # Determinar tipo de tratamiento
            treatment_formulation = int_row['treatment_formulation']
            tipo_tratamiento = 'MI'  # Default
            if 'myo_inositol' in treatment_formulation and 'd_chiro_inositol' in treatment_formulation:
                tipo_tratamiento = 'MI+DCI'
            elif 'd_chiro_inositol' in treatment_formulation:
                tipo_tratamiento = 'DCI'
            elif 'metformin' in treatment_formulation:
                tipo_tratamiento = 'MI+MET'
            
            # Determinar tipo de control
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
                'tipo_tratamiento': tipo_tratamiento,
                'tipo_control': tipo_control,
                'obesidad': 'Sí' if int_row['has_obesity'] == 'yes' else 'No',
                'n': int_row['sample_size'],
                'efecto': g,
                'error_estandar': se,
                'ci_inferior': ci_lower,
                'ci_superior': ci_upper,
                'p_valor': p_value,
                'cambio_absoluto': int_change,
                'cambio_porcentual': int_pct_change,
                'cambio_absoluto_control': ctrl_change,
                'cambio_porcentual_control': ctrl_pct_change,
                'vs_control': True
            }
            
            resultados.append(result_data)
    
    # Convertir a DataFrame
    resultados_df = pd.DataFrame(resultados)
    
    if len(resultados_df) == 0:
        print("No se encontraron datos para procesar")
        return
    
    print(f"Datos procesados: {len(resultados_df)} comparaciones")
    print()
    
    # Crear forest plots combinados para cada parámetro
    parametros_principales = ['bmi_kg_m2', 'glucose_mg_dl', 'homa_ir', 'insulin_uU_ml', 'testosterone_ng_ml']
    
    # Crear un PDF con todos los forest plots
    with PdfPages(f'{output_dir}/forest_plots_combinados.pdf') as pdf:
        for parametro in parametros_principales:
            param_data = resultados_df[resultados_df['parametro'] == parametro]
            
            if len(param_data) < 2:
                print(f"Datos insuficientes para {parametro}")
                continue
            
            # Verificar que tenemos datos de grupos de obesidad
            has_normal = len(param_data[param_data['obesidad'] == 'No']) > 0
            has_obese = len(param_data[param_data['obesidad'] == 'Sí']) > 0
            
            if not (has_normal or has_obese):
                print(f"Sin datos de grupos de obesidad para {parametro}")
                continue
            
            print(f"Creando forest plot combinado para {PARAMETROS.get(parametro, parametro)}...")
            
            # Crear el forest plot combinado mejorado
            fig = crear_forest_plot_combinado(
                resultados_df,
                parametro,
                "Forest Plot Combinado Mejorado",
                comparacion=True,
                color_by='tipo_tratamiento'
            )
            
            if fig is not None:
                pdf.savefig(fig, bbox_inches='tight')
                print(f"✓ Forest plot creado para {parametro}")
                
                # También guardar como imagen individual
                fig.savefig(f'{output_dir}/forest_plot_combinado_{parametro}.png', 
                           dpi=300, bbox_inches='tight')
            else:
                print(f"✗ Error al crear forest plot para {parametro}")
    
    print()
    print("=== Creando Tablas Resumen con Valores Absolutos ===")
    
    # Crear tablas resumen con valores absolutos
    for parametro in parametros_principales:
        param_data = resultados_df[resultados_df['parametro'] == parametro]
        
        if len(param_data) < 3:
            continue
        
        print(f"Creando tabla resumen para {PARAMETROS.get(parametro, parametro)}...")
        
        # Crear tabla con valores absolutos de Hedges' g
        tabla_abs = crear_tabla_resumen_con_abs(resultados_df, parametro, vs_control=True)
        
        if tabla_abs is not None:
            # Guardar tabla
            filename = f'{output_dir}/tabla_abs_{parametro}.csv'
            tabla_abs.to_csv(filename, index=False)
            print(f"✓ Tabla guardada: {filename}")
            
            # Mostrar muestra de la tabla
            print("  Muestra de la tabla:")
            print(tabla_abs.head().to_string(index=False))
            print()
        else:
            print(f"✗ Error al crear tabla para {parametro}")
    
    print("=== Resumen de Mejoras Implementadas ===")
    print("✓ Valores absolutos |g| mostrados para mayor claridad")
    print("✓ Signos preservados en cambios absolutos y porcentuales")
    print("✓ Mejor espaciado vertical entre grupos de obesidad")
    print("✓ Separación clara con encabezados de grupo")
    print("✓ Posicionamiento mejorado de anotaciones")
    print("✓ Tablas CSV actualizadas con valores absolutos")
    print()
    print(f"Resultados guardados en: {output_dir}/")
    print("- forest_plots_combinados.pdf: Todos los forest plots en un PDF")
    print("- forest_plot_combinado_*.png: Imágenes individuales")
    print("- tabla_abs_*.csv: Tablas resumen con valores absolutos")

if __name__ == "__main__":
    main()