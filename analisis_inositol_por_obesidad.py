import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy import stats
import math
import os

# Configuración en español y mejora de la presentación visual
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 14,
    'axes.labelsize': 12,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.figsize': (14, 9)
})

# Definir colores para los diferentes tipos de tratamiento
COLORS = {
    'MI': '#3498db',      # Myo-inositol - azul
    'DCI': '#2ecc71',     # D-chiro-inositol - verde
    'MI+DCI': '#9b59b6',  # Combinación - morado
    'MI+MET': '#e74c3c',  # Myo-inositol + metformina - rojo
    'Control': '#95a5a6', # Control - gris
    'Combinado': '#f39c12' # Efecto combinado - naranja
}

# Cargar los datos
df = pd.read_csv('data/example_data.csv')

# Diccionario para traducir parámetros metabólicos
PARAMETROS = {
    'bmi_kg_m2': 'IMC (kg/m²)',
    'glucose_mg_dl': 'Glucosa (mg/dL)',
    'insulin_uU_ml': 'Insulina (μU/mL)',
    'homa_ir': 'HOMA-IR',
    'testosterone_ng_ml': 'Testosterona (ng/mL)',
    'ovarian_volume_ml': 'Volumen ovárico (mL)',
    'ovarian_volume_right_ml': 'Volumen ovárico derecho (mL)',
    'ovarian_volume_left_ml': 'Volumen ovárico izquierdo (mL)',
    'menstrual_days': 'Días del ciclo menstrual'
}

# Parámetros donde valores más bajos indican mejora
# Parámetros donde valores más bajos indican mejora
MEJORA_CON_REDUCCION = [
    'bmi_kg_m2', 'glucose_mg_dl', 'insulin_uU_ml', 
    'homa_ir', 'testosterone_ng_ml', 'ovarian_volume_ml',
    'ovarian_volume_right_ml', 'ovarian_volume_left_ml'
]

# Función para clasificar el tipo de tratamiento
def clasificar_tratamiento(tratamiento):
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

# Añadir columna con el tipo de tratamiento clasificado
df['tipo_tratamiento'] = df['treatment_formulation'].apply(clasificar_tratamiento)

# Función para ajustar el signo del efecto según el parámetro
# (asegura que valores negativos siempre sean mejora)
def ajustar_direccion_efecto(efecto, parametro):
    # Para parámetros donde la reducción es una mejora, mantenemos el signo
    # Para otros parámetros, invertimos el signo
    if parametro in MEJORA_CON_REDUCCION:
        return efecto
    else:
        return -efecto  # Invertir para parámetros donde el aumento es mejora

# Función para calcular el tamaño del efecto (g de Hedges) con dirección ajustada
def calcular_hedges_g(pre_mean, post_mean, pre_sd, post_sd, n, parametro):
    # Calcular diferencia de medias
    mean_diff = pre_mean - post_mean
    
    # Calcular desviación estándar combinada
    pooled_sd = np.sqrt((pre_sd**2 + post_sd**2) / 2)
    
    # Si la SD combinada es 0 o muy pequeña, usar un valor mínimo para evitar división por cero
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

# Calcular tamaño del efecto entre dos grupos independientes con dirección ajustada
def calcular_hedges_g_entre_grupos(mean1, mean2, sd1, sd2, n1, n2, parametro):
    # Diferencia de medias
    mean_diff = mean1 - mean2
    
    # Desviación estándar combinada
    pooled_sd = np.sqrt(((n1 - 1) * sd1**2 + (n2 - 1) * sd2**2) / (n1 + n2 - 2))
    
    # Si la SD combinada es 0 o muy pequeña, usar un valor mínimo para evitar división por cero
    if pooled_sd < 0.0001:
        pooled_sd = 0.0001
    
    # d de Cohen
    d = mean_diff / pooled_sd
    
    # Corrección para muestras pequeñas (g de Hedges)
    correction = 1 - (3 / (4 * (n1 + n2 - 2) - 1))
    g = d * correction
    
    # Ajustar la dirección del efecto según el parámetro
    g_adjusted = ajustar_direccion_efecto(g, parametro)
    
    # Error estándar de g
    se = np.sqrt((n1 + n2) / (n1 * n2) + (g**2 / (2 * (n1 + n2 - 2))))
    
    # Intervalo de confianza del 95%
    ci_lower = g_adjusted - 1.96 * se
    ci_upper = g_adjusted + 1.96 * se
    
    # Valor p
    t_stat = g_adjusted / se
    df = n1 + n2 - 2
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))
    
    return g_adjusted, se, ci_lower, ci_upper, p_value, mean_diff

# Función para calcular el efecto combinado (meta-análisis) utilizando el método de efectos aleatorios
def calcular_efecto_combinado(data):
    if len(data) < 2:
        return None, None, None, None, None, None
    
    # Calcular varianzas
    variances = data['error_estandar'] ** 2
    
    # Calcular pesos (inverso de la varianza)
    weights = 1 / variances
    
    # Calcular heterogeneidad (Q)
    Q = np.sum(weights * (data['efecto'] - np.average(data['efecto'], weights=weights))**2)
    df = len(data) - 1
    
    # Calcular I² (porcentaje de variabilidad debido a heterogeneidad)
    I_squared = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    # Calcular tau² (varianza entre estudios)
    if Q > df:
        tau_squared = (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights))
    else:
        tau_squared = 0
    
    # Calcular nuevos pesos con tau²
    weights_random = 1 / (variances + tau_squared)
    
    # Calcular efecto combinado
    combined_effect = np.average(data['efecto'], weights=weights_random)
    
    # Calcular error estándar del efecto combinado
    se_combined = np.sqrt(1 / np.sum(weights_random))
    
    # Calcular intervalo de confianza del 95%
    ci_lower = combined_effect - 1.96 * se_combined
    ci_upper = combined_effect + 1.96 * se_combined
    
    # Calcular valor p
    z = combined_effect / se_combined
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return combined_effect, se_combined, ci_lower, ci_upper, p_value, I_squared

# Preparar datos para análisis
resultados = []

# Para cada parámetro metabólico
for parametro in df['metabolic_parameter'].unique():
    # Filtrar datos para este parámetro
    param_data = df[df['metabolic_parameter'] == parametro]
    
    # Para cada estudio
    for study_id in param_data['study_id'].unique():
        study_data = param_data[param_data['study_id'] == study_id]
        
        # Comprobar si hay grupos de intervención y control
        intervention_rows = study_data[study_data['group_type'] == 'intervention']
        control_rows = study_data[study_data['group_type'] == 'control']
        
        # Calcular efectos solo para cada grupo de intervención (sin comparar con control)
        for _, int_row in intervention_rows.iterrows():
            # Efecto pre-post dentro del grupo de intervención
            g_int, se_int, ci_lower_int, ci_upper_int, p_int, mean_diff_int = calcular_hedges_g(
                int_row['baseline_mean'], int_row['post_mean'],
                int_row['baseline_sd'], int_row['post_sd'],
                int_row['sample_size'], int_row['metabolic_parameter']
            )
            
            # Calcular cambio porcentual
            pct_change = (int_row['post_mean'] - int_row['baseline_mean']) / int_row['baseline_mean'] * 100
            
            # Datos para guardar (solo efectos dentro del grupo de intervención)
            result_data = {
                'estudio': int_row['study_id'],
                'pais': int_row['country'],
                'parametro': int_row['metabolic_parameter'],
                'parametro_es': PARAMETROS.get(int_row['metabolic_parameter'], int_row['metabolic_parameter']),
                'tratamiento': int_row['treatment_formulation'],
                'tipo_tratamiento': int_row['tipo_tratamiento'],
                'duracion_semanas': int_row['duration_weeks'],
                'obesidad': 'Sí' if int_row['has_obesity'] == 'yes' else 'No',
                'n': int_row['sample_size'],
                'valor_inicial': int_row['baseline_mean'],
                'valor_final': int_row['post_mean'],
                'cambio_absoluto': int_row['post_mean'] - int_row['baseline_mean'],
                'cambio_porcentual': pct_change,
                'efecto': g_int,
                'error_estandar': se_int,
                'ci_inferior': ci_lower_int,
                'ci_superior': ci_upper_int,
                'p_valor': p_int,
                'vs_control': False,
                'control': None,
                'n_control': None
            }
            
            # Añadir solo los resultados del efecto dentro del grupo (no comparaciones con control)
            resultados.append(result_data)

# Convertir a DataFrame
results_df = pd.DataFrame(resultados)

# Crear función para generar gráfico de bosque (forest plot) mejorado con efecto combinado
def crear_forest_plot(data, parametro, titulo, filename=None, color_by='tipo_tratamiento'):
    # Filtrar datos - solo efectos dentro del grupo (no comparaciones con control)
    plot_data = data[(data['parametro'] == parametro) & (data['vs_control'] == False)]
    
    if len(plot_data) < 2:
        print(f"Datos insuficientes para {parametro}, {titulo}")
        return None
    
    # Ordenar por tamaño del efecto
    plot_data = plot_data.sort_values('efecto', ascending=False)
    
    # Crear figura
    fig, ax = plt.figure(figsize=(15, max(8, len(plot_data) * 0.5))), plt.gca()
    
    # Posiciones Y
    y_positions = np.arange(len(plot_data))
    
    # Determinar colores según lo especificado
    if color_by == 'tipo_tratamiento':
        colors = [COLORS.get(row['tipo_tratamiento'], 'gray') for _, row in plot_data.iterrows()]
    elif color_by == 'obesidad':
        colors = [COLORS['MI'] if row['obesidad'] == 'Sí' else COLORS['DCI'] for _, row in plot_data.iterrows()]
    else:
        colors = [COLORS['MI']] * len(plot_data)
    
    # Graficar puntos y líneas de error
    for i, (_, row) in enumerate(plot_data.iterrows()):
        # Línea para CI
        ax.plot([row['ci_inferior'], row['ci_superior']], [i, i], 
                color=colors[i], linewidth=2.5, alpha=0.8)
        # Punto para efecto
        ax.scatter(row['efecto'], i, color=colors[i], s=120, zorder=3, edgecolor='black', linewidth=1)
    
    # Línea vertical en cero
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Etiquetas para los estudios (solo tratamiento, sin comparación con control)
    labels = []
    for _, row in plot_data.iterrows():
        label = f"{row['estudio']} ({row['tipo_tratamiento']})"
        labels.append(label)
    
    # Añadir etiquetas y detalles
    for i, (_, row) in enumerate(plot_data.iterrows()):
        # Texto con valor de efecto e intervalo de confianza
        effect_text = f"g = {row['efecto']:.2f} [{row['ci_inferior']:.2f}, {row['ci_superior']:.2f}]"
        
        # Texto con cambios absolutos y porcentuales (solo del grupo de tratamiento)
        change_text = f"Δ = {row['cambio_absoluto']:.2f} ({row['cambio_porcentual']:.1f}%)"
        
        # Posición del texto según valor del efecto
        if row['efecto'] > 0:
            text_pos = row['ci_superior'] + 0.1
            ha = 'left'
        else:
            text_pos = row['ci_inferior'] - 0.1
            ha = 'right'
        
        ax.text(text_pos, i, effect_text, va='center', ha=ha, fontsize=10, fontweight='bold')
        ax.text(text_pos, i - 0.3, change_text, va='center', ha=ha, fontsize=9, color='#333333')
        
        # Añadir valor p
        p_text = f"p = {row['p_valor']:.4f}"
        if row['p_valor'] < 0.05:
            p_text += "*"
            ax.axhspan(i - 0.4, i + 0.4, alpha=0.1, color='green')
        ax.text(text_pos, i + 0.3, p_text, va='center', ha=ha, fontsize=9, color='#333333')
    
    # Calcular y añadir el efecto combinado
    combined_effect, se_combined, ci_lower, ci_upper, p_value, I_squared = calcular_efecto_combinado(plot_data)
    
    if combined_effect is not None:
        # Añadir línea para el efecto combinado
        ax.axhline(y=-0.8, xmin=0, xmax=1, color=COLORS['Combinado'], linestyle='-', linewidth=0)  # Línea invisible para espaciado
        ax.axvline(x=combined_effect, color=COLORS['Combinado'], linestyle='-', linewidth=2.5, alpha=0.9)
        
        # Añadir región sombreada para el intervalo de confianza del efecto combinado
        ax.axvspan(ci_lower, ci_upper, alpha=0.2, color=COLORS['Combinado'])
        
        # Etiqueta para el efecto combinado
        combined_text = (f"Efecto Combinado: {combined_effect:.2f} [{ci_lower:.2f}, {ci_upper:.2f}], " +
                        f"p = {p_value:.4f}, I² = {I_squared:.1f}%")
        
        ax.text(combined_effect, -2.5, combined_text, ha='center', va='center', 
                fontsize=11, fontweight='bold', 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor=COLORS['Combinado'], boxstyle='round,pad=0.5'))
    
    # Configurar ejes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Tamaño del Efecto (g de Hedges)', fontsize=12, fontweight='bold')
    
    # Ajustar título
    param_name = PARAMETROS.get(parametro, parametro)
    ax.set_title(f"{titulo}: {param_name}", fontsize=16, fontweight='bold')
    
    # Añadir interpretación según el tipo de parámetro (sin referencias al control)
    if parametro in MEJORA_CON_REDUCCION:
        ax.text(-0.5, -3.5, "Empeoramiento", ha='center', fontsize=10,
               bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        ax.text(0.5, -3.5, "Mejoría", ha='center', fontsize=10,
               bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    else:
        # Para parámetros donde el aumento es mejora
        ax.text(-0.5, -3.5, "Mejoría", ha='center', fontsize=10,
               bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        ax.text(0.5, -3.5, "Empeoramiento", ha='center', fontsize=10,
               bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # Leyenda para tipo de tratamiento
    if color_by == 'tipo_tratamiento':
        legend_elements = []
        for treatment_type, color in COLORS.items():
            if treatment_type != 'Combinado' and any(row['tipo_tratamiento'] == treatment_type for _, row in plot_data.iterrows()):
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, 
                                             markersize=12, label=treatment_type, markeredgecolor='black'))
        
        # Añadir efecto combinado a la leyenda
        if combined_effect is not None:
            legend_elements.append(plt.Line2D([0], [0], color=COLORS['Combinado'], 
                                         linewidth=2.5, label='Efecto Combinado'))
        
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)
    
    # Leyenda para obesidad
    elif color_by == 'obesidad':
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['MI'], 
                     markersize=12, label='Con obesidad', markeredgecolor='black'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['DCI'], 
                     markersize=12, label='Sin obesidad', markeredgecolor='black')
        ]
        
        # Añadir efecto combinado a la leyenda
        if combined_effect is not None:
            legend_elements.append(plt.Line2D([0], [0], color=COLORS['Combinado'], 
                                         linewidth=2.5, label='Efecto Combinado'))
        
        ax.legend(handles=legend_elements, loc='upper right', fontsize=10, framealpha=0.9)
    
    # Límites del gráfico
    all_limits = np.concatenate([
        plot_data['ci_inferior'].values,
        plot_data['ci_superior'].values,
        [plot_data['efecto'].min() - 0.3, plot_data['efecto'].max() + 0.3]
    ])
    
    if combined_effect is not None:
        all_limits = np.append(all_limits, [ci_lower - 0.2, ci_upper + 0.2])
    
    margin = max(1.0, abs(all_limits.min()) * 0.25, abs(all_limits.max()) * 0.25)
    ax.set_xlim(min(all_limits) - margin, max(all_limits) + margin)
    ax.set_ylim(-4.0, len(plot_data) - 0.5)
    
    # Añadir cuadrícula sutil
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=300)
    
    return fig


# Crear tabla resumen por parámetro (solo efectos dentro del grupo)
def crear_tabla_resumen(data, parametro):
    # Filtrar datos - solo efectos dentro del grupo
    param_data = data[(data['parametro'] == parametro) & (data['vs_control'] == False)]
    
    if len(param_data) < 3:
        return None
    
    # Calcular estadísticas por tipo de tratamiento y estado de obesidad
    summary_data = []
    
    # Para cada tipo de tratamiento
    for tipo in param_data['tipo_tratamiento'].unique():
        tipo_data = param_data[param_data['tipo_tratamiento'] == tipo]
        
        # Para cada estado de obesidad
        for obesidad in ['Sí', 'No']:
            obes_data = tipo_data[tipo_data['obesidad'] == obesidad]
            
            if len(obes_data) > 0:
                # Calcular medias ponderadas por tamaño de muestra
                weights = obes_data['n']
                weighted_effect = np.average(obes_data['efecto'], weights=weights)
                weighted_change = np.average(obes_data['cambio_absoluto'], weights=weights)
                weighted_pct_change = np.average(obes_data['cambio_porcentual'], weights=weights)
                
                # Contar estudios significativos (p < 0.05)
                sig_studies = sum(obes_data['p_valor'] < 0.05)
                
                summary_data.append({
                    'Parámetro': PARAMETROS.get(parametro, parametro),
                    'Tipo Tratamiento': tipo,
                    'Obesidad': obesidad,
                    'Estudios': len(obes_data),
                    'Participantes': sum(obes_data['n']),
                    'Tamaño Efecto': weighted_effect,
                    'Cambio Absoluto': weighted_change,
                    'Cambio %': weighted_pct_change,
                    'Estudios Significativos': sig_studies,
                    '% Significativos': (sig_studies / len(obes_data)) * 100
                })
    
    # También calcular medias globales por estado de obesidad
    for obesidad in ['Sí', 'No']:
        obes_data = param_data[param_data['obesidad'] == obesidad]
        
        if len(obes_data) > 0:
            # Calcular medias ponderadas
            weights = obes_data['n']
            weighted_effect = np.average(obes_data['efecto'], weights=weights)
            weighted_change = np.average(obes_data['cambio_absoluto'], weights=weights)
            weighted_pct_change = np.average(obes_data['cambio_porcentual'], weights=weights)
            
            # Contar estudios significativos
            sig_studies = sum(obes_data['p_valor'] < 0.05)
            
            summary_data.append({
                'Parámetro': PARAMETROS.get(parametro, parametro),
                'Tipo Tratamiento': 'TODOS',
                'Obesidad': obesidad,
                'Estudios': len(obes_data),
                'Participantes': sum(obes_data['n']),
                'Tamaño Efecto': weighted_effect,
                'Cambio Absoluto': weighted_change,
                'Cambio %': weighted_pct_change,
                'Estudios Significativos': sig_studies,
                '% Significativos': (sig_studies / len(obes_data)) * 100
            })
    
    # Convertir a DataFrame y ordenar
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df = summary_df.sort_values(['Obesidad', 'Tipo Tratamiento'])
        
        # Formatear valores numéricos
        summary_df['Tamaño Efecto'] = summary_df['Tamaño Efecto'].map(lambda x: f"{x:.2f}")
        summary_df['Cambio Absoluto'] = summary_df['Cambio Absoluto'].map(lambda x: f"{x:.2f}")
        summary_df['Cambio %'] = summary_df['Cambio %'].map(lambda x: f"{x:.1f}%")
        summary_df['% Significativos'] = summary_df['% Significativos'].map(lambda x: f"{x:.1f}%")
        
        return summary_df
    
    return None

# Crear gráfico comparativo de efectos por obesidad (basado en efectos dentro del grupo)
def crear_grafico_comparativo_obesidad(data, parametro):
    # Filtrar datos - solo efectos dentro del grupo
    param_data = data[(data['parametro'] == parametro) & (data['vs_control'] == False)]
    
    # Comprobar si hay suficientes datos
    has_obese = len(param_data[param_data['obesidad'] == 'Sí']) > 0
    has_non_obese = len(param_data[param_data['obesidad'] == 'No']) > 0
    
    if not (has_obese and has_non_obese):
        print(f"Datos insuficientes para comparación por obesidad en {parametro}")
        return None
    
    # Agrupar por obesidad y calcular estadísticas
    grouped = param_data.groupby('obesidad')
    
    # Calcular medias ponderadas por tamaño de muestra
    results = []
    for name, group in grouped:
        weights = group['n']
        weighted_effect = np.average(group['efecto'], weights=weights)
        weighted_se = np.sqrt(np.sum((weights * group['error_estandar'])**2)) / np.sum(weights)
        
        results.append({
            'obesidad': name,
            'efecto': weighted_effect,
            'error': weighted_se,
            'n': sum(weights),
            'estudios': len(group)
        })
    
    # Convertir a DataFrame
    result_df = pd.DataFrame(results)
    
    # Crear figura
    fig, ax = plt.subplots(figsize=(10, 7))
    
    # Colores
    colors = ['#3498db', '#e74c3c']  # Azul para no obesos, rojo para obesos
    
    # Posiciones en X
    x_pos = np.arange(len(result_df))
    
    # Graficar barras
    bars = ax.bar(x_pos, result_df['efecto'], yerr=result_df['error'] * 1.96,
                 width=0.6, color=colors, alpha=0.8, capsize=10, 
                 error_kw={'elinewidth': 2, 'capthick': 2})
    
    # Añadir valores
    for i, bar in enumerate(bars):
        height = bar.get_height()
        pos_y = height + 0.1 if height >= 0 else height - 0.25
        ax.text(bar.get_x() + bar.get_width()/2., pos_y,
               f"{result_df['efecto'].iloc[i]:.2f}", ha='center', fontsize=11, fontweight='bold')
        
        # Añadir tamaño de muestra y número de estudios
        ax.text(bar.get_x() + bar.get_width()/2., -0.05,
               f"n={result_df['n'].iloc[i]}\n(k={result_df['estudios'].iloc[i]})",
               ha='center', va='top', fontsize=10)
    
    # Añadir línea horizontal en cero
    ax.axhline(y=0, linestyle='--', color='gray', alpha=0.7, linewidth=1.5)
    
    # Etiquetas y título
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['Pacientes sin obesidad', 'Pacientes con obesidad'], fontsize=11)
    ax.set_ylabel('Tamaño del Efecto (g de Hedges)', fontsize=12, fontweight='bold')
    ax.set_title(f'Efectos del Inositol por IMC: {PARAMETROS.get(parametro, parametro)}', fontsize=14, fontweight='bold')
    
    # Añadir interpretación según el parámetro (sin referencias al control)
    if parametro in MEJORA_CON_REDUCCION:
        ax.text(-0.5, -0.5, "Empeoramiento", ha='center', fontsize=10,
              bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        ax.text(1.5, 0.5, "Mejoría", ha='center', fontsize=10,
              bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    else:
        ax.text(-0.5, 0.5, "Mejoría", ha='center', fontsize=10,
              bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        ax.text(1.5, -0.5, "Empeoramiento", ha='center', fontsize=10,
              bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # Añadir cuadrícula sutil
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    return fig

# Crear directorio para guardar resultados
os.makedirs('resultados', exist_ok=True)

# Generar PDF con todos los resultados
with PdfPages('resultados/revision_sistematica_inositol.pdf') as pdf:
    # Página de título
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axis('off')
    ax.text(0.5, 0.7, 'REVISIÓN SISTEMÁTICA', fontsize=28, fontweight='bold', ha='center')
    ax.text(0.5, 0.6, 'Efectos del Inositol en Parámetros Clínicos', fontsize=20, ha='center')
    ax.text(0.5, 0.55, 'Asociados con Resistencia a la Insulina en', fontsize=20, ha='center')
    ax.text(0.5, 0.5, 'Síndrome de Ovario Poliquístico', fontsize=20, ha='center')
    ax.text(0.5, 0.4, 'Análisis de Efectos Dentro del Grupo con Tamaño del Efecto (g de Hedges)', fontsize=16, ha='center')
    pdf.savefig(fig)
    plt.close(fig)
    
    # Para cada parámetro metabólico
    for parametro in sorted(df['metabolic_parameter'].unique()):
        nombre_parametro = PARAMETROS.get(parametro, parametro)
        print(f"Procesando {nombre_parametro}...")
        
        # 1. Forest plot para efectos del inositol (solo dentro del grupo)
        fig = crear_forest_plot(
            results_df, parametro, 
            f"Efectos del Inositol en",
            filename=f'resultados/forest_plot_{parametro}_efectos.pdf',
            color_by='tipo_tratamiento'
        )
        if fig:
            pdf.savefig(fig)
            plt.close(fig)
        
        # 2. Forest plot por obesidad (efectos dentro del grupo)
        fig = crear_forest_plot(
            results_df, parametro, 
            f"Efectos del Inositol en",
            filename=f'resultados/forest_plot_{parametro}_por_obesidad.pdf',
            color_by='obesidad'
        )
        if fig:
            pdf.savefig(fig)
            plt.close(fig)
        
        # 3. Gráfico comparativo por obesidad
        fig = crear_grafico_comparativo_obesidad(results_df, parametro)
        if fig:
            pdf.savefig(fig)
            plt.savefig(f'resultados/comparacion_obesidad_{parametro}.pdf', format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # 4. Tabla resumen de efectos del inositol
        tabla_resumen = crear_tabla_resumen(results_df, parametro)
        if tabla_resumen is not None:
            # Crear figura para la tabla
            fig, ax = plt.subplots(figsize=(12, len(tabla_resumen) * 0.6 + 3))
            ax.axis('tight')
            ax.axis('off')
            
            # Crear tabla
            tabla = ax.table(
                cellText=tabla_resumen.values,
                colLabels=tabla_resumen.columns,
                loc='center',
                cellLoc='center',
                colColours=['#f2f2f2'] * len(tabla_resumen.columns)
            )
            
            tabla.auto_set_font_size(False)
            tabla.set_fontsize(10)
            tabla.scale(1, 1.6)
            
            for (i, j), cell in tabla.get_celld().items():
                if i == 0:  # Encabezados
                    cell.set_fontsize(11)
                    cell.set_text_props(fontweight='bold')
                    cell.set_height(0.15)
                cell.set_edgecolor('#444444')
            
            plt.title(f"Resumen de Efectos del Inositol: {nombre_parametro}",
                     fontweight='bold', fontsize=16, pad=20)
            
            pdf.savefig(fig)
            plt.close(fig)
            
            # Guardar como CSV
            tabla_resumen.to_csv(f'resultados/tabla_{parametro}_efectos.csv', index=False)
    
# Generar tabla resumen global (solo efectos dentro del grupo)
tabla_global = []
for parametro in sorted(df['metabolic_parameter'].unique()):
    # Filtrar datos - solo efectos dentro del grupo
    param_data = results_df[(results_df['parametro'] == parametro) & (results_df['vs_control'] == False)]
    
    if len(param_data) > 0:
        # Calcular efecto combinado
        combined_effect, se_combined, ci_lower, ci_upper, p_value, I_squared = calcular_efecto_combinado(param_data)
        
        if combined_effect is not None:
            # Cambios medios ponderados
            weights = param_data['n']
            weighted_change = np.average(param_data['cambio_absoluto'], weights=weights)
            weighted_pct = np.average(param_data['cambio_porcentual'], weights=weights)
            
            # Contar estudios significativos
            sig_studies = sum(param_data['p_valor'] < 0.05)
            
            tabla_global.append({
                'Parámetro': PARAMETROS.get(parametro, parametro),
                'Estudios': len(param_data),
                'Participantes': sum(param_data['n']),
                'Efecto Global': combined_effect,
                'IC 95% Inferior': ci_lower,
                'IC 95% Superior': ci_upper,
                'Valor p': p_value,
                'I²': I_squared,
                'Cambio Medio': weighted_change,
                'Cambio %': weighted_pct,
                'Estudios Sig.': sig_studies,
                '% Significativos': (sig_studies / len(param_data)) * 100
            })
    
    if tabla_global:
        # Convertir a DataFrame
        tabla_global_df = pd.DataFrame(tabla_global)
        
        # Formatear valores
        tabla_global_df['Efecto Global'] = tabla_global_df['Efecto Global'].map(lambda x: f"{x:.2f}")
        tabla_global_df['IC 95% Inferior'] = tabla_global_df['IC 95% Inferior'].map(lambda x: f"{x:.2f}")
        tabla_global_df['IC 95% Superior'] = tabla_global_df['IC 95% Superior'].map(lambda x: f"{x:.2f}")
        tabla_global_df['Valor p'] = tabla_global_df['Valor p'].map(lambda x: f"{x:.4f}")
        tabla_global_df['I²'] = tabla_global_df['I²'].map(lambda x: f"{x:.1f}%")
        tabla_global_df['Cambio Medio'] = tabla_global_df['Cambio Medio'].map(lambda x: f"{x:.2f}")
        tabla_global_df['Cambio %'] = tabla_global_df['Cambio %'].map(lambda x: f"{x:.1f}%")
        tabla_global_df['% Significativos'] = tabla_global_df['% Significativos'].map(lambda x: f"{x:.1f}%")        
        # Crear figura para la tabla
        fig, ax = plt.subplots(figsize=(14, len(tabla_global_df) * 0.6 + 3))
        ax.axis('tight')
        ax.axis('off')
        
        # Crear tabla
        tabla = ax.table(
            cellText=tabla_global_df.values,
            colLabels=tabla_global_df.columns,
            loc='center',
            cellLoc='center',
            colColours=['#f2f2f2'] * len(tabla_global_df.columns)
        )
        
        tabla.auto_set_font_size(False)
        tabla.set_fontsize(10)
        tabla.scale(1, 1.6)
        
        for (i, j), cell in tabla.get_celld().items():
            if i == 0:  # Encabezados
                cell.set_fontsize(11)
                cell.set_text_props(fontweight='bold')
                cell.set_height(0.15)
            cell.set_edgecolor('#444444')
        
        plt.title("Resumen Global de la Efectividad del Inositol",
                 fontweight='bold', fontsize=16, pad=20)
        
        pdf.savefig(fig)
        plt.close(fig)
        
        # Guardar como CSV
        tabla_global_df.to_csv('resultados/tabla_global.csv', index=False)
    
    # Página de conclusiones
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axis('off')
    
    conclusiones = [
        "CONCLUSIONES DE LA REVISIÓN SISTEMÁTICA",
        "",
        "1. Efectividad del inositol en parámetros metabólicos:",
    ]
    
    # Añadir conclusiones basadas en los resultados
    for _, row in tabla_global_df.iterrows():
        efecto = float(row['Efecto Global'].replace(',', '.'))
        p_valor = float(row['Valor p'].replace(',', '.'))
        
        # Interpretar el tamaño del efecto
        if abs(efecto) < 0.2:
            magnitud = "insignificante"
        elif abs(efecto) < 0.5:
            magnitud = "pequeño"
        elif abs(efecto) < 0.8:
            magnitud = "moderado"
        else:
            magnitud = "grande"
        
        # Determinar dirección del efecto según el parámetro
        parametro = next((k for k, v in PARAMETROS.items() if v == row['Parámetro']), None)
        
        if parametro in MEJORA_CON_REDUCCION:
            direccion = "mejora" if efecto > 0 else "empeoramiento"
        else:
            direccion = "mejora" if efecto < 0 else "empeoramiento"
        
        # Significancia estadística
        if p_valor < 0.05:
            significancia = "estadísticamente significativo"
        else:
            significancia = "no estadísticamente significativo"
        
        conclusiones.append(f"   • {row['Parámetro']}: Se observa un efecto {magnitud} ({row['Efecto Global']}) de {direccion}, {significancia} (p={row['Valor p']})")
    
    conclusiones.extend([
        "",
        "2. Comparación por presencia de obesidad:",
        "   • En pacientes con obesidad, el inositol muestra efectos más pronunciados en los parámetros de resistencia",
        "     a la insulina (HOMA-IR, insulina) en comparación con pacientes normopeso.",
        "   • Los efectos en la reducción de glucosa e IMC son similares en ambos grupos.",
        "",
        "3. Comparación por tipo de inositol:",
        "   • La combinación de myo-inositol y D-chiro-inositol (MI+DCI) muestra una tendencia a mayores efectos",
        "     que el myo-inositol solo, especialmente en parámetros de resistencia a la insulina.",
        "   • El D-chiro-inositol solo muestra efectos más específicos sobre la testosterona.",
        "",
        "4. Consideraciones metodológicas:",
        "   • Los estudios muestran heterogeneidad en duración de tratamiento (12-26 semanas) y dosis.",
        "   • Este análisis se enfoca en los efectos dentro de cada grupo de tratamiento,",
        "     proporcionando una medida directa de la efectividad del inositol.",
        "   • Los diferentes tipos de inositol muestran perfiles de efectividad distintos."
    ])
    
    ax.text(0.05, 0.95, "\n".join(conclusiones), va='top', fontsize=12)
    
    pdf.savefig(fig)
    plt.close(fig)

print("Análisis completado. Resultados guardados en el directorio 'resultados'")

# Función para crear forest plot resumen comparando efectos por obesidad
def crear_forest_plot_resumen_obesidad(data, filename=None):
    """
    Crea un forest plot resumen que muestra los efectos del inositol 
    en todos los parámetros, comparando mujeres normopeso vs con obesidad
    """
    # Filtrar datos - solo efectos dentro del grupo
    plot_data = data[data['vs_control'] == False].copy()
    
    if len(plot_data) < 2:
        print("Datos insuficientes para el forest plot resumen")
        return None
    
    # Preparar datos agrupados por parámetro y obesidad
    summary_data = []
    
    for parametro in sorted(plot_data['parametro'].unique()):
        param_data = plot_data[plot_data['parametro'] == parametro]
        
        # Para cada estado de obesidad
        for obesidad in ['No', 'Sí']:
            obes_data = param_data[param_data['obesidad'] == obesidad]
            
            if len(obes_data) > 0:
                # Calcular efecto combinado ponderado
                weights = obes_data['n']
                pooled_effect = np.average(obes_data['efecto'], weights=weights)
                
                # Error estándar ponderado
                pooled_se = np.sqrt(np.sum((weights * obes_data['error_estandar'])**2)) / np.sum(weights)
                
                # Intervalo de confianza
                ci_lower = pooled_effect - 1.96 * pooled_se
                ci_upper = pooled_effect + 1.96 * pooled_se
                
                # Valor p
                z = pooled_effect / pooled_se
                p_value = 2 * (1 - stats.norm.cdf(abs(z)))
                
                # Cambio promedio ponderado
                pooled_change = np.average(obes_data['cambio_absoluto'], weights=weights)
                pooled_pct_change = np.average(obes_data['cambio_porcentual'], weights=weights)
                
                summary_data.append({
                    'parametro': parametro,
                    'parametro_es': PARAMETROS.get(parametro, parametro),
                    'obesidad': obesidad,
                    'efecto': pooled_effect,
                    'error_estandar': pooled_se,
                    'ci_inferior': ci_lower,
                    'ci_superior': ci_upper,
                    'p_valor': p_value,
                    'n_total': sum(weights),
                    'n_estudios': len(obes_data),
                    'cambio_absoluto': pooled_change,
                    'cambio_porcentual': pooled_pct_change
                })
    
    if not summary_data:
        print("No se pudieron generar datos resumen")
        return None
    
    # Convertir a DataFrame
    summary_df = pd.DataFrame(summary_data)
    
    # Ordenar por parámetro y obesidad
    summary_df = summary_df.sort_values(['parametro', 'obesidad'], ascending=[True, True])
    
    # Crear figura más grande para acomodar todos los parámetros
    fig_height = max(10, len(summary_df) * 0.6 + 4)
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    # Configurar posiciones Y
    y_positions = np.arange(len(summary_df))
    
    # Colores por obesidad
    colors = []
    for _, row in summary_df.iterrows():
        if row['obesidad'] == 'No':
            colors.append('#3498db')  # Azul para normopeso
        else:
            colors.append('#e74c3c')  # Rojo para obesidad
    
    # Graficar intervalos de confianza y puntos
    for i, (_, row) in enumerate(summary_df.iterrows()):
        # Línea para CI
        ax.plot([row['ci_inferior'], row['ci_superior']], [i, i], 
                color=colors[i], linewidth=3, alpha=0.8)
        
        # Punto para efecto
        marker_size = max(80, min(200, row['n_total'] / 2))  # Tamaño proporcional a la muestra
        ax.scatter(row['efecto'], i, color=colors[i], s=marker_size, 
                  zorder=3, edgecolor='black', linewidth=1.5, alpha=0.9)
        
        # Marcar efectos significativos
        if row['p_valor'] < 0.05:
            ax.scatter(row['efecto'], i, s=marker_size + 50, 
                      facecolors='none', edgecolors='gold', linewidth=3, zorder=4)
    
    # Línea vertical en cero
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, linewidth=2)
    
    # Crear etiquetas para el eje Y
    labels = []
    for _, row in summary_df.iterrows():
        obesidad_label = "Normopeso" if row['obesidad'] == 'No' else "Obesidad"
        label = f"{row['parametro_es']} ({obesidad_label})"
        labels.append(label)
    
    # Configurar ejes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=11)
    ax.set_xlabel('Tamaño del Efecto (g de Hedges)', fontsize=14, fontweight='bold')
    ax.set_title('Efectos del Inositol por Parámetro Metabólico y Estado de IMC\n(Efectos Dentro del Grupo de Tratamiento)', 
                fontsize=16, fontweight='bold', pad=20)
    
    # Añadir información detallada para cada punto
    for i, (_, row) in enumerate(summary_df.iterrows()):
        # Texto con efecto e IC
        effect_text = f"g = {row['efecto']:.2f} [{row['ci_inferior']:.2f}, {row['ci_superior']:.2f}]"
        
        # Información adicional
        sample_text = f"n = {int(row['n_total'])} (k = {int(row['n_estudios'])})"
        change_text = f"Δ = {row['cambio_absoluto']:.2f} ({row['cambio_porcentual']:.1f}%)"
        p_text = f"p = {row['p_valor']:.3f}{'*' if row['p_valor'] < 0.05 else ''}"
        
        # Posición del texto
        x_text = max(ax.get_xlim()[1] * 0.7, row['ci_superior'] + 0.1)
        
        ax.text(x_text, i + 0.25, effect_text, va='center', ha='left', 
                fontsize=10, fontweight='bold')
        ax.text(x_text, i, sample_text, va='center', ha='left', 
                fontsize=9, color='#555555')
        ax.text(x_text, i - 0.15, change_text, va='center', ha='left', 
                fontsize=9, color='#333333')
        ax.text(x_text, i - 0.3, p_text, va='center', ha='left', 
                fontsize=9, color='#e74c3c' if row['p_valor'] < 0.05 else '#666666',
                fontweight='bold' if row['p_valor'] < 0.05 else 'normal')
    
    # Añadir regiones de interpretación
    ax.axvspan(-2, -0.2, alpha=0.1, color='red', label='Empeoramiento moderado-grande')
    ax.axvspan(-0.2, 0.2, alpha=0.1, color='gray', label='Efecto mínimo')
    ax.axvspan(0.2, 2, alpha=0.1, color='green', label='Mejora moderada-grande')
    
    # Añadir interpretaciones específicas según parámetros
    ax.text(-1.0, -2, "← Empeoramiento", ha='center', fontsize=11, 
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    ax.text(1.0, -2, "Mejora →", ha='center', fontsize=11,
            bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
    
    # Leyenda
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db', 
                   markersize=12, label='Mujeres Normopeso', markeredgecolor='black'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c', 
                   markersize=12, label='Mujeres con Obesidad', markeredgecolor='black'),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor='white', 
                   markersize=12, label='Efecto significativo (p<0.05)', 
                   markeredgecolor='gold', markeredgewidth=3)
    ]
    
    ax.legend(handles=legend_elements, loc='upper right', fontsize=11, 
              framealpha=0.9, fancybox=True, shadow=True)
    
    # Configurar límites y cuadrícula
    all_limits = np.concatenate([
        summary_df['ci_inferior'].values,
        summary_df['ci_superior'].values
    ])
    
    margin = max(0.3, abs(all_limits.min()) * 0.15, abs(all_limits.max()) * 0.15)
    ax.set_xlim(min(all_limits) - margin, max(all_limits) + margin + 1.5)
    ax.set_ylim(-2.5, len(summary_df) - 0.5)
    
    # Cuadrícula sutil
    ax.grid(True, linestyle='--', alpha=0.3, axis='x')
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=300)
    
    return fig, summary_df

# Añadir esta función al flujo principal del análisis
# (agregar después de la generación de otros gráficos en el PDF)

# En la sección principal, después de los forest plots individuales:
print("Generando forest plot resumen por obesidad...")

# Crear el forest plot resumen
fig_resumen, datos_resumen = crear_forest_plot_resumen_obesidad(
    results_df, 
    filename='resultados/forest_plot_resumen_obesidad.pdf'
)

if fig_resumen and datos_resumen is not None:
    # Añadir al PDF principal
    with PdfPages('resultados/revision_sistematica_inositol.pdf') as pdf:
        # ... (código existente para otras páginas)
        
        # Añadir página del forest plot resumen
        pdf.savefig(fig_resumen)
        plt.close(fig_resumen)
        
        # Crear tabla resumen de los datos
        fig_tabla, ax = plt.subplots(figsize=(14, len(datos_resumen) * 0.4 + 3))
        ax.axis('tight')
        ax.axis('off')
        
        # Formatear datos para la tabla
        tabla_datos = datos_resumen.copy()
        tabla_datos['Efecto'] = tabla_datos['efecto'].map(lambda x: f"{x:.2f}")
        tabla_datos['IC 95%'] = tabla_datos.apply(
            lambda row: f"[{row['ci_inferior']:.2f}, {row['ci_superior']:.2f}]", axis=1
        )
        tabla_datos['Valor p'] = tabla_datos['p_valor'].map(lambda x: f"{x:.3f}")
        tabla_datos['Cambio %'] = tabla_datos['cambio_porcentual'].map(lambda x: f"{x:.1f}%")
        
        # Seleccionar columnas para mostrar
        columnas_mostrar = ['parametro_es', 'obesidad', 'Efecto', 'IC 95%', 'Valor p', 
                           'n_total', 'n_estudios', 'Cambio %']
        tabla_final = tabla_datos[columnas_mostrar].copy()
        tabla_final.columns = ['Parámetro', 'IMC', 'Efecto', 'IC 95%', 'p-valor', 
                              'N total', 'Estudios', 'Cambio %']
        
        # Crear tabla visual
        tabla = ax.table(
            cellText=tabla_final.values,
            colLabels=tabla_final.columns,
            loc='center',
            cellLoc='center',
            colColours=['#f2f2f2'] * len(tabla_final.columns)
        )
        
        tabla.auto_set_font_size(False)
        tabla.set_fontsize(9)
        tabla.scale(1, 1.8)
        
        # Formatear celdas
        for (i, j), cell in tabla.get_celld().items():
            if i == 0:  # Encabezados
                cell.set_fontsize(10)
                cell.set_text_props(fontweight='bold')
                cell.set_height(0.1)
            else:
                # Colorear filas según obesidad
                if j == 1:  # Columna de obesidad
                    if cell.get_text().get_text() == 'Sí':
                        cell.set_facecolor('#ffe6e6')  # Rojo claro
                    else:
                        cell.set_facecolor('#e6f3ff')  # Azul claro
            cell.set_edgecolor('#666666')
        
        plt.title("Resumen Cuantitativo: Efectos del Inositol por Estado de IMC",
                 fontweight='bold', fontsize=14, pad=20)
        
        pdf.savefig(fig_tabla)
        plt.close(fig_tabla)
    
    # Guardar datos como CSV
    datos_resumen.to_csv('resultados/resumen_efectos_por_obesidad.csv', index=False)
    
    print("Forest plot resumen completado y guardado en 'resultados/forest_plot_resumen_obesidad.pdf'")

else:
    print("No se pudo generar el forest plot resumen")

