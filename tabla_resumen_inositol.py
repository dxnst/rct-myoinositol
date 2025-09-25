#!/usr/bin/env python3
"""
Script para analizar los efectos del inositol en diferentes variables metabólicas.
Se enfoca únicamente en el análisis pre-post de los grupos de intervención con inositol.

Autor: Script generado para análisis de efectos del inositol
Fecha: 2025
"""

import pandas as pd
import numpy as np
from scipy import stats
import os

# Diccionario para traducir parámetros metabólicos
PARAMETROS = {
    'bmi_kg_m2': 'IMC (kg/m²)',
    'glucose_mg_dl': 'Glucosa (mg/dL)',
    'insulin_uU_ml': 'Insulina (μU/mL)',
    'homa_ir': 'HOMA-IR',
    'testosterone_ng_ml': 'Testosterona (ng/mL)',
    'ferriman_gallwey': 'Puntuación Ferriman-Gallwey',
    'ovarian_volume_ml': 'Volumen ovárico (mL)',
    'ovarian_volume_right_ml': 'Volumen ovárico derecho (mL)',
    'ovarian_volume_left_ml': 'Volumen ovárico izquierdo (mL)',
    'menstrual_days': 'Días del ciclo menstrual'
}

# Parámetros donde valores más bajos indican mejora
MEJORA_CON_REDUCCION = [
    'bmi_kg_m2', 'glucose_mg_dl', 'insulin_uU_ml', 
    'homa_ir', 'testosterone_ng_ml', 'ovarian_volume_ml',
    'ovarian_volume_right_ml', 'ovarian_volume_left_ml', 'ferriman_gallwey'
]

def clasificar_tratamiento(tratamiento):
    """Clasifica el tipo de tratamiento basado en la formulación"""
    if pd.isna(tratamiento):
        return 'Otro'
    tratamiento = tratamiento.lower()
    if 'myo_inositol' in tratamiento and 'd_chiro_inositol' in tratamiento:
        return 'MI+DCI'
    elif 'myo_inositol' in tratamiento and 'metformin' in tratamiento:
        return 'MI+MET'
    elif 'myo_inositol' in tratamiento:
        return 'MI'
    elif 'd_chiro_inositol' in tratamiento:
        return 'DCI'
    elif 'metformin' in tratamiento:
        return 'Control'
    elif 'placebo' in tratamiento or 'diet' in tratamiento:
        return 'Control'
    else:
        return 'Otro'

def es_tratamiento_inositol(tipo_tratamiento):
    """Determina si un tratamiento contiene inositol"""
    return tipo_tratamiento in ['MI', 'DCI', 'MI+DCI', 'MI+MET']

def ajustar_direccion_efecto(efecto, parametro):
    """Ajusta la dirección del efecto según el parámetro"""
    if parametro in MEJORA_CON_REDUCCION:
        return efecto
    else:
        return -efecto  # Invertir para parámetros donde el aumento es mejora

def calcular_hedges_g(pre_mean, post_mean, pre_sd, post_sd, n, parametro):
    """Calcula el tamaño del efecto (g de Hedges) con dirección ajustada"""
    # Calcular diferencia de medias
    mean_diff = pre_mean - post_mean
    
    # Calcular desviación estándar combinada
    pooled_sd = np.sqrt((pre_sd**2 + post_sd**2) / 2)
    
    # Si la SD combinada es 0 o muy pequeña, usar un valor mínimo
    if pooled_sd < 0.0001:
        pooled_sd = 0.0001
    
    # Calcular d de Cohen
    d = mean_diff / pooled_sd
    
    # Aplicar corrección para muestras pequeñas (g de Hedges)
    correction = 1 - (3 / (4 * (n - 1) - 1))
    g = d * correction
    
    # Ajustar la dirección del efecto según el parámetro
    g_adjusted = ajustar_direccion_efecto(g, parametro)
    
    # Error estándar de g de Hedges
    se = np.sqrt((2/n) + ((g**2) / (2*n)))
    
    # Calcular intervalo de confianza del 95%
    ci_lower = g_adjusted - 1.96 * se
    ci_upper = g_adjusted + 1.96 * se
    
    # Calcular valor p
    z = g_adjusted / se
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return g_adjusted, se, ci_lower, ci_upper, p_value, mean_diff

def interpretar_efecto(g):
    """Interpreta el tamaño del efecto según la convención de Cohen"""
    abs_g = abs(g)
    if abs_g < 0.2:
        return "Sin efecto"
    elif abs_g < 0.5:
        return "Pequeño"
    elif abs_g < 0.8:
        return "Moderado"
    else:
        return "Grande"

def analizar_datos_inositol(archivo_csv):
    """Función principal para analizar los datos del inositol"""
    
    # Leer los datos
    print(f"Leyendo datos de {archivo_csv}...")
    df = pd.read_csv(archivo_csv)
    
    print(f"Total de registros: {len(df)}")
    print(f"Columnas disponibles: {list(df.columns)}")
    
    # Clasificar tipos de tratamiento
    df['tipo_tratamiento'] = df['treatment_formulation'].apply(clasificar_tratamiento)
    
    # Filtrar solo grupos de intervención con inositol
    df_inositol = df[
        (df['group_type'] == 'intervention') & 
        (df['tipo_tratamiento'].apply(es_tratamiento_inositol))
    ].copy()
    
    print(f"Registros con intervención de inositol: {len(df_inositol)}")
    print(f"Parámetros únicos: {df_inositol['metabolic_parameter'].unique()}")
    print(f"Tipos de tratamiento con inositol: {df_inositol['tipo_tratamiento'].unique()}")
    
    # Preparar lista para resultados
    resultados = []
    
    # Análisis por parámetro metabólico
    for parametro in sorted(df_inositol['metabolic_parameter'].unique()):
        param_data = df_inositol[df_inositol['metabolic_parameter'] == parametro]
        
        print(f"\nAnalizando {parametro}: {len(param_data)} estudios")
        
        # Análisis por estudio
        efectos_estudios = []
        tamaños_muestra = []
        diferencias_medias = []
        
        for _, row in param_data.iterrows():
            # Calcular efecto pre-post
            g, se, ci_lower, ci_upper, p_value, mean_diff = calcular_hedges_g(
                row['baseline_mean'], row['post_mean'],
                row['baseline_sd'], row['post_sd'],
                row['sample_size'], row['metabolic_parameter']
            )
            
            efectos_estudios.append(g)
            tamaños_muestra.append(row['sample_size'])
            diferencias_medias.append(mean_diff)
            
            print(f"  Estudio {row['study_id']}: g={g:.3f}, n={row['sample_size']}, diff={mean_diff:.2f}")
        
        # Calcular efecto combinado ponderado por tamaño de muestra
        if efectos_estudios:
            weights = np.array(tamaños_muestra)
            efecto_ponderado = np.average(efectos_estudios, weights=weights)
            diferencia_ponderada = np.average(diferencias_medias, weights=weights)
            
            # Error estándar combinado (aproximación simple)
            se_combinado = np.sqrt(np.average([(2/n) + (g**2/(2*n)) for g, n in zip(efectos_estudios, tamaños_muestra)], weights=weights))
            
            # Intervalo de confianza combinado
            ci_lower_comb = efecto_ponderado - 1.96 * se_combinado
            ci_upper_comb = efecto_ponderado + 1.96 * se_combinado
            
            # Agregar a resultados
            resultados.append({
                'Variable': PARAMETROS.get(parametro, parametro),
                'Numero_Estudios': len(param_data),
                'Tamaño_Muestra_Total': sum(tamaños_muestra),
                'Diferencia_Medias': diferencia_ponderada,
                'Tamaño_Efecto_Hedges_g': efecto_ponderado,
                'Interpretacion': interpretar_efecto(efecto_ponderado),
                'IC_95_Inferior': ci_lower_comb,
                'IC_95_Superior': ci_upper_comb,
                'Parametro_Original': parametro
            })
    
    # Crear DataFrame de resultados
    df_resultados = pd.DataFrame(resultados)
    
    # Ordenar por tamaño del efecto (valor absoluto, descendente)
    df_resultados['Abs_Efecto'] = abs(df_resultados['Tamaño_Efecto_Hedges_g'])
    df_resultados = df_resultados.sort_values('Abs_Efecto', ascending=False)
    df_resultados = df_resultados.drop('Abs_Efecto', axis=1)
    
    return df_resultados

def generar_tabla_resumen(df_resultados, guardar_csv=True):
    """Genera y muestra la tabla resumen"""
    
    print("\n" + "="*120)
    print("TABLA RESUMEN: EFECTOS DEL INOSITOL EN VARIABLES METABÓLICAS")
    print("Análisis pre-post de grupos de intervención con inositol únicamente")
    print("="*120)
    
    # Crear tabla formateada para mostrar
    tabla_display = df_resultados.copy()
    tabla_display['Diferencia_Medias'] = tabla_display['Diferencia_Medias'].round(2)
    tabla_display['Tamaño_Efecto_Hedges_g'] = tabla_display['Tamaño_Efecto_Hedges_g'].round(3)
    tabla_display['IC_95_Inferior'] = tabla_display['IC_95_Inferior'].round(3)
    tabla_display['IC_95_Superior'] = tabla_display['IC_95_Superior'].round(3)
    tabla_display['Intervalo_Confianza_95%'] = tabla_display.apply(
        lambda row: f"[{row['IC_95_Inferior']:.3f}, {row['IC_95_Superior']:.3f}]", axis=1
    )
    
    # Seleccionar columnas para display
    columnas_display = [
        'Variable', 'Numero_Estudios', 'Tamaño_Muestra_Total', 
        'Diferencia_Medias', 'Tamaño_Efecto_Hedges_g', 'Interpretacion', 
        'Intervalo_Confianza_95%'
    ]
    
    tabla_final = tabla_display[columnas_display].copy()
    tabla_final.columns = [
        'Variable Analizada', 'Nº Estudios', 'Tamaño Muestra', 
        'Diferencia Medias Pre-Post', 'Tamaño Efecto (g Hedges)', 
        'Interpretación Efecto', 'Intervalo Confianza 95%'
    ]
    
    # Mostrar tabla
    print(tabla_final.to_string(index=False, max_colwidth=30))
    
    # Guardar en CSV si se solicita
    if guardar_csv:
        # Crear directorio output si no existe
        os.makedirs('output', exist_ok=True)
        
        archivo_csv = 'output/tabla_resumen_efectos_inositol.csv'
        df_resultados.to_csv(archivo_csv, index=False)
        print(f"\n✓ Resultados guardados en: {archivo_csv}")
        
        # También guardar la tabla formateada
        archivo_csv_formatted = 'output/tabla_resumen_efectos_inositol_formatted.csv'
        tabla_final.to_csv(archivo_csv_formatted, index=False)
        print(f"✓ Tabla formateada guardada en: {archivo_csv_formatted}")
    
    print("\n" + "="*120)
    print("RESUMEN DEL ANÁLISIS:")
    print(f"• Total de variables analizadas: {len(df_resultados)}")
    print(f"• Variables con efecto grande: {sum(df_resultados['Interpretacion'] == 'Grande')}")
    print(f"• Variables con efecto moderado: {sum(df_resultados['Interpretacion'] == 'Moderado')}")
    print(f"• Variables con efecto pequeño: {sum(df_resultados['Interpretacion'] == 'Pequeño')}")
    print(f"• Variables sin efecto: {sum(df_resultados['Interpretacion'] == 'Sin efecto')}")
    print("="*120)
    
    return tabla_final, df_resultados

def main():
    """Función principal"""
    try:
        # Archivo de datos
        archivo_csv = 'data/example_data.csv'
        
        if not os.path.exists(archivo_csv):
            print(f"Error: No se encontró el archivo {archivo_csv}")
            return
        
        # Analizar datos
        df_resultados = analizar_datos_inositol(archivo_csv)
        
        if df_resultados.empty:
            print("No se encontraron datos para analizar.")
            return
        
        # Generar tabla resumen
        tabla_final, df_completo = generar_tabla_resumen(df_resultados)
        
        return tabla_final, df_completo
        
    except Exception as e:
        print(f"Error durante el análisis: {str(e)}")
        import traceback
        traceback.print_exc()
        return None, None

if __name__ == "__main__":
    main()