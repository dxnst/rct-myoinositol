import pandas as pd
import numpy as np
from scipy import stats
import math

# Función para calcular g de Hedges
def calcular_hedges_g(mean_diff, sd_pooled, n1, n2):
    # Calcular d de Cohen
    d = mean_diff / sd_pooled
    
    # Factor de corrección para g de Hedges
    j = 1 - (3 / (4 * (n1 + n2 - 2) - 1))
    g = j * d
    
    # Error estándar
    se = np.sqrt((n1 + n2) / (n1 * n2) + (g**2) / (2 * (n1 + n2 - 2)))
    
    # Intervalo de confianza del 95%
    ci_lower = g - 1.96 * se
    ci_upper = g + 1.96 * se
    
    return g, ci_lower, ci_upper

# Función para interpretar el tamaño del efecto
def interpretar_efecto(g):
    if abs(g) < 0.2:
        return "Sin efecto"
    elif abs(g) < 0.5:
        return "Pequeño"
    elif abs(g) < 0.8:
        return "Moderado"
    else:
        return "Grande"

# Cargar los datos
df = pd.read_csv('data/example_data.csv')

# Traducir los nombres de los parámetros
parametros_dict = {
    'bmi_kg_m2': 'Índice de masa corporal (kg/m2)',
    'homa_ir': 'HOMA-IR (índice de resistencia a la insulina)',
    'insulin_uU_ml': 'Insulina en ayunas (µU/mL)',
    'glucose_mg_dl': 'Glucosa en ayunas (mg/dL)',
    'testosterone_ng_ml': 'Testosterona total (ng/dL)',
    'ovarian_volume_ml': 'Volumen ovárico (mL)',
    'ovarian_volume_left_ml': 'Volumen ovárico izquierdo (mL)',
    'ovarian_volume_right_ml': 'Volumen ovárico derecho (mL)',
    'menstrual_days': 'Ciclo menstrual (ciclos por trimestre)',
    'ferriman_gallwey': 'Puntuación Ferriman-Gallwey'
}

# Inicializar lista para almacenar los resultados
resultados = []

# Procesar cada parámetro metabólico
for parametro in df['metabolic_parameter'].unique():
    # Filtrar datos para este parámetro
    param_data = df[df['metabolic_parameter'] == parametro]
    
    # Separar grupos de intervención y control
    intervenciones = param_data[param_data['group_type'] == 'intervention']
    controles = param_data[param_data['group_type'] == 'control']
    
    # Contar estudios únicos
    estudios_unicos = param_data['study_id'].unique()
    num_estudios = len(estudios_unicos)
    
    # Si no hay suficientes datos, continuar con el siguiente parámetro
    if len(intervenciones) == 0 or len(controles) == 0:
        continue
    
    # Calcular tamaño de muestra total
    n_intervencion = intervenciones['sample_size'].sum()
    n_control = controles['sample_size'].sum()
    
    # Calcular estadísticas para cada estudio y luego promediar
    efectos = []
    diferencias_medias = []
    
    for estudio in estudios_unicos:
        estudio_int = intervenciones[intervenciones['study_id'] == estudio]
        estudio_ctrl = controles[controles['study_id'] == estudio]
        
        # Si no hay datos de intervención o control para este estudio, omitirlo
        if len(estudio_int) == 0 or len(estudio_ctrl) == 0:
            continue
        
        # Para simplificar, tomamos solo la primera entrada si hay múltiples
        int_row = estudio_int.iloc[0]
        ctrl_row = estudio_ctrl.iloc[0]
        
        # Calcular cambios para cada grupo
        cambio_int = int_row['post_mean'] - int_row['baseline_mean']
        cambio_ctrl = ctrl_row['post_mean'] - ctrl_row['baseline_mean']
        
        # Calcular diferencia de medias
        diff_medias = cambio_int - cambio_ctrl
        diferencias_medias.append(diff_medias)
        
        # Calcular desviaciones estándar de los cambios
        sd_cambio_int = np.sqrt(int_row['baseline_sd']**2 + int_row['post_sd']**2)
        sd_cambio_ctrl = np.sqrt(ctrl_row['baseline_sd']**2 + ctrl_row['post_sd']**2)
        
        # Calcular desviación estándar combinada
        n1 = int_row['sample_size']
        n2 = ctrl_row['sample_size']
        sd_pooled = np.sqrt(((n1-1)*sd_cambio_int**2 + (n2-1)*sd_cambio_ctrl**2) / (n1+n2-2))
        
        # Calcular g de Hedges
        g, ci_lower, ci_upper = calcular_hedges_g(diff_medias, sd_pooled, n1, n2)
        
        efectos.append({
            'g': g,
            'ci_lower': ci_lower,
            'ci_upper': ci_upper,
            'n_int': n1,
            'n_ctrl': n2
        })
    
    # Si no hay suficientes datos para calcular efectos, continuar con el siguiente parámetro
    if len(efectos) == 0:
        continue
    
    # Calcular promedios ponderados por tamaño de muestra
    total_n = sum(e['n_int'] + e['n_ctrl'] for e in efectos)
    g_promedio = sum(e['g'] * (e['n_int'] + e['n_ctrl']) for e in efectos) / total_n
    
    # Calcular intervalo de confianza combinado usando meta-análisis simple
    # (este es un enfoque simplificado)
    se_combinado = np.sqrt(1 / total_n)
    ci_lower_combinado = g_promedio - 1.96 * se_combinado
    ci_upper_combinado = g_promedio + 1.96 * se_combinado
    
    # Calcular diferencia de medias promedio
    diff_medias_promedio = np.mean(diferencias_medias)
    
    # Determinar interpretación del efecto
    interpretacion = interpretar_efecto(g_promedio)
    
    # Agregar resultado a la lista
    resultados.append({
        'Variable': parametros_dict.get(parametro, parametro),
        'Número de estudios': num_estudios,
        'Grupo control (n)': n_control,
        'Grupo intervención (n)': n_intervencion,
        'Diferencia de medias entre grupos': round(diff_medias_promedio, 2),
        'Tamaño del efecto (g de Hedges)': round(g_promedio, 2),
        'Interpretación del efecto': interpretacion,
        'Intervalo de Confianza 95%': f"({round(ci_lower_combinado, 2)}, {round(ci_upper_combinado, 2)})"
    })

# Convertir a DataFrame
resultados_df = pd.DataFrame(resultados)

# Guardar como CSV
resultados_df.to_csv('resultados_inositol.csv', index=False)

# Imprimir la tabla
print(resultados_df.to_string(index=False))
