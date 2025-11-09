#!/usr/bin/env python3
"""
Script para análisis de efectos del inositol en variables metabólicas
Analiza cambios pre-post en grupos de intervención con inositol
Calcula tamaños del efecto usando g de Hedges y realiza meta-análisis
"""

import pandas as pd
import numpy as np
import scipy.stats as stats
from scipy.stats import t
import warnings
warnings.filterwarnings('ignore')

def calculate_hedges_g(baseline_mean, baseline_sd, post_mean, post_sd, n):
    """
    Calcula g de Hedges para cambios pre-post dentro del mismo grupo
    
    Parameters:
    -----------
    baseline_mean : float
        Media basal
    baseline_sd : float  
        Desviación estándar basal
    post_mean : float
        Media post-tratamiento
    post_sd : float
        Desviación estándar post-tratamiento
    n : int
        Tamaño de muestra
    
    Returns:
    --------
    tuple: (hedges_g, variance, se)
    """
    # Diferencia de medias
    mean_diff = post_mean - baseline_mean
    
    # Desviación estándar pooled para medidas repetidas
    # Asumiendo correlación moderada r=0.5 entre pre y post
    r = 0.5
    pooled_sd = np.sqrt(baseline_sd**2 + post_sd**2 - 2*r*baseline_sd*post_sd)
    
    # Cohen's d
    cohens_d = mean_diff / pooled_sd
    
    # Corrección de Hedges para muestras pequeñas
    correction_factor = 1 - (3 / (4 * (n - 1) - 1))
    hedges_g = cohens_d * correction_factor
    
    # Varianza de g de Hedges para medidas repetidas
    variance = (1/n) + (hedges_g**2 / (2*n)) + (2*(1-r)/n)
    se = np.sqrt(variance)
    
    return hedges_g, variance, se

def meta_analysis_random_effects(effect_sizes, variances):
    """
    Realiza meta-análisis de efectos aleatorios usando método DerSimonian-Laird
    
    Parameters:
    -----------
    effect_sizes : array-like
        Tamaños del efecto individuales
    variances : array-like
        Varianzas de los tamaños del efecto
    
    Returns:
    --------
    dict: Resultados del meta-análisis
    """
    effect_sizes = np.array(effect_sizes)
    variances = np.array(variances)
    weights = 1 / variances
    
    # Efecto combinado con efectos fijos
    pooled_effect_fixed = np.sum(weights * effect_sizes) / np.sum(weights)
    
    # Estadístico Q para heterogeneidad
    Q = np.sum(weights * (effect_sizes - pooled_effect_fixed)**2)
    df = len(effect_sizes) - 1
    
    # Tau-squared (varianza entre estudios)
    if Q > df:
        tau_squared = (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights))
        tau_squared = max(0, tau_squared)
    else:
        tau_squared = 0
    
    # Pesos para efectos aleatorios
    weights_random = 1 / (variances + tau_squared)
    
    # Efecto combinado con efectos aleatorios
    pooled_effect = np.sum(weights_random * effect_sizes) / np.sum(weights_random)
    pooled_variance = 1 / np.sum(weights_random)
    pooled_se = np.sqrt(pooled_variance)
    
    # Intervalo de confianza del 95%
    df_pooled = len(effect_sizes) - 1
    t_crit = t.ppf(0.975, df_pooled)
    ci_lower = pooled_effect - t_crit * pooled_se
    ci_upper = pooled_effect + t_crit * pooled_se
    
    # Valor p (prueba bilateral)
    t_stat = pooled_effect / pooled_se
    p_value = 2 * (1 - t.cdf(abs(t_stat), df_pooled))
    
    return {
        'pooled_effect': pooled_effect,
        'pooled_se': pooled_se,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'p_value': p_value,
        'Q': Q,
        'df': df,
        'tau_squared': tau_squared,
        'I_squared': max(0, (Q - df) / Q * 100) if Q > 0 else 0
    }

def interpret_effect_size(effect_size):
    """
    Interpreta el tamaño del efecto según criterios de Cohen
    """
    abs_effect = abs(effect_size)
    if abs_effect < 0.2:
        return "Sin efecto"
    elif abs_effect < 0.5:
        return "Pequeño"
    elif abs_effect < 0.8:
        return "Moderado"
    else:
        return "Grande"

def analyze_inositol_effects(csv_file_path):
    """
    Función principal para analizar efectos del inositol
    """
    print("🔍 Cargando y procesando datos...")
    
    # Cargar datos
    df = pd.read_csv(csv_file_path)
    
    # Filtrar solo grupos de intervención con inositol
    inositol_data = df[
        (df['group_type'] == 'intervention') & 
        (df['treatment_formulation'].str.contains('inositol', case=False, na=False))
    ].copy()
    
    print(f"📊 Encontrados {len(inositol_data)} registros de intervención con inositol")
    print(f"📋 Variables metabólicas encontradas: {inositol_data['metabolic_parameter'].nunique()}")
    
    # Agrupar por variable metabólica
    results = []
    
    for variable in inositol_data['metabolic_parameter'].unique():
        var_data = inositol_data[inositol_data['metabolic_parameter'] == variable].copy()
        
        # Calcular g de Hedges para cada estudio
        effect_sizes = []
        variances = []
        valid_studies = []
        total_sample_size = 0
        mean_differences = []
        
        for idx, row in var_data.iterrows():
            try:
                hedges_g, variance, se = calculate_hedges_g(
                    row['baseline_mean'], row['baseline_sd'],
                    row['post_mean'], row['post_sd'],
                    row['sample_size']
                )
                
                effect_sizes.append(hedges_g)
                variances.append(variance)
                valid_studies.append(row['study_id'])
                total_sample_size += row['sample_size']
                mean_differences.append(row['post_mean'] - row['baseline_mean'])
                
            except Exception as e:
                print(f"⚠️  Error procesando {row['study_id']} - {variable}: {e}")
                continue
        
        if not effect_sizes:
            continue
            
        # Realizar meta-análisis si hay múltiples estudios
        if len(effect_sizes) > 1:
            meta_results = meta_analysis_random_effects(effect_sizes, variances)
            pooled_effect = meta_results['pooled_effect']
            ci_lower = meta_results['ci_lower']
            ci_upper = meta_results['ci_upper']
            p_value = meta_results['p_value']
        else:
            # Un solo estudio
            pooled_effect = effect_sizes[0]
            se = np.sqrt(variances[0])
            t_crit = t.ppf(0.975, total_sample_size - 1)
            ci_lower = pooled_effect - t_crit * se
            ci_upper = pooled_effect + t_crit * se
            t_stat = pooled_effect / se
            p_value = 2 * (1 - t.cdf(abs(t_stat), total_sample_size - 1))
        
        # Calcular diferencia de medias promedio ponderada por tamaño de muestra
        weights = [var_data.iloc[i]['sample_size'] for i in range(len(mean_differences))]
        weighted_mean_diff = np.average(mean_differences, weights=weights)
        
        results.append({
            'Variable': variable,
            'Num_Estudios': len(effect_sizes),
            'Tamaño_Muestra_Total': total_sample_size,
            'Diferencia_Medias_Promedio': weighted_mean_diff,
            'Hedges_g': pooled_effect,
            'IC_95_Inferior': ci_lower,
            'IC_95_Superior': ci_upper,
            'Valor_P': p_value,
            'Interpretacion': interpret_effect_size(pooled_effect),
            'Estudios': ', '.join(valid_studies)
        })
    
    # Crear DataFrame con resultados
    results_df = pd.DataFrame(results)
    
    # Ordenar por tamaño del efecto (valor absoluto) de mayor a menor
    results_df['Abs_Effect'] = results_df['Hedges_g'].abs()
    results_df = results_df.sort_values('Abs_Effect', ascending=False).drop('Abs_Effect', axis=1)
    
    return results_df

def format_results_table(results_df):
    """
    Formatea la tabla de resultados para presentación
    """
    formatted_df = results_df.copy()
    
    # Redondear valores numéricos
    formatted_df['Diferencia_Medias_Promedio'] = formatted_df['Diferencia_Medias_Promedio'].round(3)
    formatted_df['Hedges_g'] = formatted_df['Hedges_g'].round(3)
    formatted_df['IC_95_Inferior'] = formatted_df['IC_95_Inferior'].round(3)
    formatted_df['IC_95_Superior'] = formatted_df['IC_95_Superior'].round(3)
    
    # Formatear valor p
    formatted_df['Valor_P'] = formatted_df['Valor_P'].apply(
        lambda x: f"{x:.4f}" if x >= 0.001 else "<0.001"
    )
    
    # Crear columna de IC 95% combinada
    formatted_df['IC_95%'] = formatted_df.apply(
        lambda row: f"[{row['IC_95_Inferior']}, {row['IC_95_Superior']}]", axis=1
    )
    
    # Seleccionar y reordenar columnas finales
    final_df = formatted_df[[
        'Variable', 'Num_Estudios', 'Tamaño_Muestra_Total', 
        'Diferencia_Medias_Promedio', 'Hedges_g', 'Interpretacion',
        'IC_95%', 'Valor_P'
    ]].copy()
    
    # Renombrar columnas para presentación
    final_df.columns = [
        'Variable Analizada',
        'N° Estudios', 
        'Tamaño Muestra Total',
        'Diferencia de Medias',
        'g de Hedges',
        'Interpretación del Efecto',
        'Intervalo Confianza 95%',
        'Valor P'
    ]
    
    return final_df

def main():
    """
    Función principal del script
    """
    print("=" * 80)
    print("📊 ANÁLISIS DE EFECTOS DEL INOSITOL EN VARIABLES METABÓLICAS")
    print("=" * 80)
    print()
    
    csv_file = "data/example_data.csv"
    
    try:
        # Realizar análisis
        results_df = analyze_inositol_effects(csv_file)
        
        if results_df.empty:
            print("❌ No se encontraron datos válidos para analizar")
            return
        
        # Formatear resultados
        formatted_results = format_results_table(results_df)
        
        print("\n" + "=" * 80)
        print("📋 TABLA RESUMEN DE EFECTOS DEL INOSITOL")
        print("=" * 80)
        print("\nNOTA: Análisis de cambios pre-post en grupos de intervención con inositol")
        print("Ordenado por tamaño del efecto (valor absoluto) de mayor a menor\n")
        
        # Mostrar tabla
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', 30)
        print(formatted_results.to_string(index=False))
        
        # Guardar resultados
        output_file = "resultados_efectos_inositol.csv"
        formatted_results.to_csv(output_file, index=False, encoding='utf-8')
        
        print(f"\n✅ Resultados guardados en: {output_file}")
        
        # Resumen estadístico
        print("\n" + "=" * 80)
        print("📊 RESUMEN ESTADÍSTICO")
        print("=" * 80)
        
        total_variables = len(formatted_results)
        efectos_significativos = len(formatted_results[formatted_results['Valor P'].str.replace('<', '').astype(float) < 0.05])
        efectos_grandes = len(formatted_results[formatted_results['Interpretación del Efecto'] == 'Grande'])
        efectos_moderados = len(formatted_results[formatted_results['Interpretación del Efecto'] == 'Moderado'])
        efectos_pequeños = len(formatted_results[formatted_results['Interpretación del Efecto'] == 'Pequeño'])
        sin_efecto = len(formatted_results[formatted_results['Interpretación del Efecto'] == 'Sin efecto'])
        
        print(f"📈 Total de variables analizadas: {total_variables}")
        print(f"🎯 Efectos estadísticamente significativos (p<0.05): {efectos_significativos}")
        print(f"🔥 Efectos grandes: {efectos_grandes}")
        print(f"🟡 Efectos moderados: {efectos_moderados}")
        print(f"🟢 Efectos pequeños: {efectos_pequeños}")
        print(f"⚪ Sin efecto: {sin_efecto}")
        
        print("\n" + "=" * 80)
        print("✨ ANÁLISIS COMPLETADO EXITOSAMENTE")
        print("=" * 80)
        
    except FileNotFoundError:
        print(f"❌ Error: No se encontró el archivo {csv_file}")
        print("   Asegúrate de que el archivo existe en la ruta especificada")
    except Exception as e:
        print(f"❌ Error inesperado: {e}")
        print("   Revisa el formato de los datos y vuelve a intentar")

if __name__ == "__main__":
    main()
