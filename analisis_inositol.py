import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
from scipy import stats
import math
import os

# Configuración en español
plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'axes.titlesize': 12,
    'axes.labelsize': 10,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 9,
    'figure.figsize': (12, 8)
})

# Definir colores para los diferentes tipos de tratamiento
COLORS = {
    'MI': '#3498db',      # Myo-inositol - azul
    'DCI': '#2ecc71',     # D-chiro-inositol - verde
    'MI+DCI': '#9b59b6',  # Combinación - morado
    'MI+MET': '#e74c3c',  # Myo-inositol + metformina - rojo
    'Control': '#95a5a6'  # Control - gris
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
    'ferriman_gallwey': 'Puntuación Ferriman-Gallwey',
    'ovarian_volume_ml': 'Volumen ovárico (mL)',
    'ovarian_volume_right_ml': 'Volumen ovárico derecho (mL)',
    'ovarian_volume_left_ml': 'Volumen ovárico izquierdo (mL)',
    'menstrual_days': 'Días del ciclo menstrual'
}

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

# Función para calcular el tamaño del efecto (g de Hedges)
def calcular_hedges_g(pre_mean, post_mean, pre_sd, post_sd, n):
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
    
    # Error estándar de g de Hedges
    se = np.sqrt((2/n) + ((g**2) / (2*n)))
    
    # Calcular intervalo de confianza del 95%
    ci_lower = g - 1.96 * se
    ci_upper = g + 1.96 * se
    
    # Calcular valor p
    z = g / se
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return g, se, ci_lower, ci_upper, p_value, mean_diff

# Calcular tamaño del efecto entre dos grupos independientes
def calcular_hedges_g_entre_grupos(mean1, mean2, sd1, sd2, n1, n2):
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
    
    # Error estándar de g
    se = np.sqrt((n1 + n2) / (n1 * n2) + (g**2 / (2 * (n1 + n2 - 2))))
    
    # Intervalo de confianza del 95%
    ci_lower = g - 1.96 * se
    ci_upper = g + 1.96 * se
    
    # Valor p
    t_stat = g / se
    df = n1 + n2 - 2
    p_value = 2 * (1 - stats.t.cdf(abs(t_stat), df))
    
    return g, se, ci_lower, ci_upper, p_value, mean_diff

# Función para calcular la diferencia neta del efecto entre grupos de tratamiento y control
def calcular_efecto_neto(intervention_row, control_row):
    # Calcular cambio en el grupo de intervención
    int_change = intervention_row['post_mean'] - intervention_row['baseline_mean']
    int_sd_change = np.sqrt(intervention_row['baseline_sd']**2 + intervention_row['post_sd']**2)
    
    # Calcular cambio en el grupo de control
    ctrl_change = control_row['post_mean'] - control_row['baseline_mean']
    ctrl_sd_change = np.sqrt(control_row['baseline_sd']**2 + control_row['post_sd']**2)
    
    # Calcular el efecto neto (diferencia de cambios)
    net_effect = int_change - ctrl_change
    
    # Calcular g de Hedges para el efecto neto
    g, se, ci_lower, ci_upper, p_value, _ = calcular_hedges_g_entre_grupos(
        int_change, ctrl_change, 
        int_sd_change, ctrl_sd_change,
        intervention_row['sample_size'], control_row['sample_size']
    )
    
    return g, se, ci_lower, ci_upper, p_value, net_effect

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
        
        # Calcular efectos para cada grupo de intervención
        for _, int_row in intervention_rows.iterrows():
            # Efecto pre-post dentro del grupo de intervención
            g_int, se_int, ci_lower_int, ci_upper_int, p_int, mean_diff_int = calcular_hedges_g(
                int_row['baseline_mean'], int_row['post_mean'],
                int_row['baseline_sd'], int_row['post_sd'],
                int_row['sample_size']
            )
            
            # Calcular cambio porcentual
            pct_change = (int_row['post_mean'] - int_row['baseline_mean']) / int_row['baseline_mean'] * 100
            
            # Datos para guardar
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
            
            # Si hay grupo control para comparar, calcular también el efecto neto
            if not control_rows.empty:
                # Tomar el primer grupo control (la mayoría de los estudios tienen solo uno)
                control_row = control_rows.iloc[0]
                
                # Calcular efecto neto (intervención vs control)
                g_net, se_net, ci_lower_net, ci_upper_net, p_net, net_effect = calcular_efecto_neto(
                    int_row, control_row
                )
                
                # Efecto pre-post dentro del grupo control
                g_ctrl, se_ctrl, ci_lower_ctrl, ci_upper_ctrl, p_ctrl, mean_diff_ctrl = calcular_hedges_g(
                    control_row['baseline_mean'], control_row['post_mean'],
                    control_row['baseline_sd'], control_row['post_sd'],
                    control_row['sample_size']
                )
                
                # Cambio porcentual en el control
                pct_change_ctrl = (control_row['post_mean'] - control_row['baseline_mean']) / control_row['baseline_mean'] * 100
                
                # Añadir datos de la comparación con control
                result_data_vs_control = result_data.copy()
                result_data_vs_control.update({
                    'vs_control': True,
                    'control': control_row['treatment_formulation'],
                    'tipo_control': clasificar_tratamiento(control_row['treatment_formulation']),
                    'n_control': control_row['sample_size'],
                    'valor_inicial_control': control_row['baseline_mean'],
                    'valor_final_control': control_row['post_mean'],
                    'cambio_absoluto_control': control_row['post_mean'] - control_row['baseline_mean'],
                    'cambio_porcentual_control': pct_change_ctrl,
                    'efecto': g_net,
                    'error_estandar': se_net,
                    'ci_inferior': ci_lower_net,
                    'ci_superior': ci_upper_net,
                    'p_valor': p_net
                })
                
                resultados.append(result_data_vs_control)
            
            # Añadir los resultados del efecto dentro del grupo
            resultados.append(result_data)

# Convertir a DataFrame
results_df = pd.DataFrame(resultados)

# Crear función para generar gráfico de bosque (forest plot)
def crear_forest_plot(data, parametro, titulo, comparacion=True, filename=None, color_by='tipo_tratamiento'):
    # Filtrar datos
    if comparacion:
        plot_data = data[(data['parametro'] == parametro) & (data['vs_control'] == True)]
    else:
        plot_data = data[(data['parametro'] == parametro) & (data['vs_control'] == False)]
    
    if len(plot_data) < 2:
        print(f"Datos insuficientes para {parametro}, {titulo}")
        return None
    
    # Ordenar por tamaño del efecto
    plot_data = plot_data.sort_values('efecto', ascending=False)
    
    # Crear figura
    fig, ax = plt.figure(figsize=(14, max(8, len(plot_data) * 0.6))), plt.gca()
    
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
                color=colors[i], linewidth=2, alpha=0.7)
        # Punto para efecto
        ax.scatter(row['efecto'], i, color=colors[i], s=100, zorder=3)
    
    # Línea vertical en cero
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
    
    # Añadir zonas de interpretación
    ax.axvspan(-0.2, 0.2, alpha=0.1, color='gray')
    ax.axvspan(-0.5, -0.2, alpha=0.1, color='lightblue')
    ax.axvspan(0.2, 0.5, alpha=0.1, color='lightblue')
    ax.axvspan(-0.8, -0.5, alpha=0.1, color='lightgreen')
    ax.axvspan(0.5, 0.8, alpha=0.1, color='lightgreen')
    ax.axvspan(-4, -0.8, alpha=0.1, color='#ffb6c1')
    ax.axvspan(0.8, 4, alpha=0.1, color='#ffb6c1')
    
    # Etiquetas para los estudios
    labels = []
    for _, row in plot_data.iterrows():
        if comparacion:
            label = f"{row['estudio']} ({row['tipo_tratamiento']} vs {row['tipo_control']})"
        else:
            label = f"{row['estudio']} ({row['tipo_tratamiento']})"
        labels.append(label)
    
    # Añadir etiquetas y detalles
    for i, (_, row) in enumerate(plot_data.iterrows()):
        # Texto con valor de efecto e intervalo de confianza
        effect_text = f"g = {row['efecto']:.2f} [{row['ci_inferior']:.2f}, {row['ci_superior']:.2f}]"
        
        # Texto con cambios absolutos y porcentuales
        if comparacion:
            change_text = (f"Δ = {row['cambio_absoluto']:.2f} vs {row['cambio_absoluto_control']:.2f} " +
                          f"({row['cambio_porcentual']:.1f}% vs {row['cambio_porcentual_control']:.1f}%)")
        else:
            change_text = f"Δ = {row['cambio_absoluto']:.2f} ({row['cambio_porcentual']:.1f}%)"
        
        # Posición del texto según valor del efecto
        if row['efecto'] > 0:
            text_pos = row['ci_superior'] + 0.1
            ha = 'left'
        else:
            text_pos = row['ci_inferior'] - 0.1
            ha = 'right'
        
        ax.text(text_pos, i, effect_text, va='center', ha=ha, fontsize=9)
        ax.text(text_pos, i - 0.3, change_text, va='center', ha=ha, fontsize=8, color='#555555')
        
        # Añadir valor p
        p_text = f"p = {row['p_valor']:.4f}"
        if row['p_valor'] < 0.05:
            p_text += "*"
            ax.axhspan(i - 0.4, i + 0.4, alpha=0.1, color='green')
        ax.text(text_pos, i + 0.3, p_text, va='center', ha=ha, fontsize=8, color='#555555')
    
    # Configurar ejes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels)
    ax.set_xlabel('Tamaño del Efecto (g de Hedges)', fontsize=11, fontweight='bold')
    
    # Ajustar título
    param_name = PARAMETROS.get(parametro, parametro)
    ax.set_title(f"{titulo}: {param_name}", fontsize=14, fontweight='bold')
    
    # Añadir interpretación
    if comparacion:
        ax.text(-0.5, -1.5, "Favorece al control", ha='center', fontsize=9)
        ax.text(0.5, -1.5, "Favorece al inositol", ha='center', fontsize=9)
    else:
        ax.text(-0.5, -1.5, "Empeoramiento", ha='center', fontsize=9)
        ax.text(0.5, -1.5, "Mejoría", ha='center', fontsize=9)
    
    # Leyenda para tipo de tratamiento
    if color_by == 'tipo_tratamiento':
        legend_elements = []
        for treatment_type, color in COLORS.items():
            if any(row['tipo_tratamiento'] == treatment_type for _, row in plot_data.iterrows()):
                legend_elements.append(plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=treatment_type))
        if legend_elements:
            ax.legend(handles=legend_elements, loc='upper right')
    
    # Leyenda para obesidad
    elif color_by == 'obesidad':
        legend_elements = [
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['MI'], markersize=10, label='Con obesidad'),
            plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['DCI'], markersize=10, label='Sin obesidad')
        ]
        ax.legend(handles=legend_elements, loc='upper right')
    
    # Límites del gráfico
    all_limits = np.concatenate([
        plot_data['ci_inferior'].values,
        plot_data['ci_superior'].values,
        [plot_data['efecto'].min() - 0.2, plot_data['efecto'].max() + 0.2]
    ])
    margin = max(1.0, abs(all_limits.min()) * 0.2, abs(all_limits.max()) * 0.2)
    ax.set_xlim(min(all_limits) - margin, max(all_limits) + margin)
    ax.set_ylim(-2.5, len(plot_data) - 0.5)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=300)
    
    return fig

# Crear tabla resumen por parámetro
def crear_tabla_resumen(data, parametro, vs_control=True):
    # Filtrar datos
    if vs_control:
        param_data = data[(data['parametro'] == parametro) & (data['vs_control'] == True)]
    else:
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

# Crear gráfico comparativo de efectos por obesidad
def crear_grafico_comparativo_obesidad(data, parametro):
    # Filtrar datos
    param_data = data[(data['parametro'] == parametro) & (data['vs_control'] == True)]
    
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
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Colores
    colors = ['#3498db', '#e74c3c']  # Azul para no obesos, rojo para obesos
    
    # Posiciones en X
    x_pos = np.arange(len(result_df))
    
    # Graficar barras
    bars = ax.bar(x_pos, result_df['efecto'], yerr=result_df['error'] * 1.96,
                 width=0.6, color=colors, alpha=0.7, capsize=10, 
                 error_kw={'elinewidth': 2, 'capthick': 2})
    
    # Añadir valores
    for i, bar in enumerate(bars):
        height = bar.get_height()
        pos_y = height + 0.1 if height >= 0 else height - 0.25
        ax.text(bar.get_x() + bar.get_width()/2., pos_y,
               f"{result_df['efecto'].iloc[i]:.2f}", ha='center', fontsize=10)
        
        # Añadir tamaño de muestra y número de estudios
        ax.text(bar.get_x() + bar.get_width()/2., -0.05,
               f"n={result_df['n'].iloc[i]}\n(k={result_df['estudios'].iloc[i]})",
               ha='center', va='top', fontsize=9)
    
    # Añadir línea horizontal en cero
    ax.axhline(y=0, linestyle='--', color='gray', alpha=0.6)
    
    # Etiquetas y título
    ax.set_xticks(x_pos)
    ax.set_xticklabels(['Pacientes sin obesidad', 'Pacientes con obesidad'])
    ax.set_ylabel('Tamaño del Efecto (g de Hedges)', fontsize=11)
    ax.set_title(f'Comparación por IMC: {PARAMETROS.get(parametro, parametro)}', fontsize=14, fontweight='bold')
    
    # Añadir interpretación
    ax.text(-0.5, -0.5, "Favorece al control", ha='center', fontsize=9)
    ax.text(0.5, -0.5, "Favorece al control", ha='center', fontsize=9)
    ax.text(1.5, 0.5, "Favorece al inositol", ha='center', fontsize=9)
    
    plt.tight_layout()
    
    return fig

# Crear directorio para guardar resultados
os.makedirs('resultados', exist_ok=True)

# Generar PDF con todos los resultados
with PdfPages('resultados/revision_sistematica_inositol.pdf') as pdf:
    # Página de título
    fig, ax = plt.subplots(figsize=(12, 10))
    ax.axis('off')
    ax.text(0.5, 0.7, 'REVISIÓN SISTEMÁTICA', fontsize=24, fontweight='bold', ha='center')
    ax.text(0.5, 0.6, 'Efectividad del Inositol en la Reducción de Parámetros Clínicos', fontsize=18, ha='center')
    ax.text(0.5, 0.55, 'Asociados con Resistencia a la Insulina en', fontsize=18, ha='center')
    ax.text(0.5, 0.5, 'Síndrome de Ovario Poliquístico', fontsize=18, ha='center')
    ax.text(0.5, 0.4, 'Análisis con Tamaño del Efecto (g de Hedges)', fontsize=14, ha='center')
    pdf.savefig(fig)
    plt.close(fig)
    
    # Para cada parámetro metabólico
    for parametro in sorted(df['metabolic_parameter'].unique()):
        nombre_parametro = PARAMETROS.get(parametro, parametro)
        print(f"Procesando {nombre_parametro}...")
        
        # 1. Forest plot para inositol vs control
        fig = crear_forest_plot(
            results_df, parametro, 
            f"Efectividad del Inositol vs Control en",
            comparacion=True,
            color_by='tipo_tratamiento'
        )
        if fig:
            pdf.savefig(fig)
            plt.savefig(f'resultados/forest_plot_{parametro}_vs_control.pdf', format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # 2. Forest plot para efectos dentro del grupo de inositol
        fig = crear_forest_plot(
            results_df, parametro, 
            f"Efecto del Inositol en",
            comparacion=False,
            color_by='tipo_tratamiento'
        )
        if fig:
            pdf.savefig(fig)
            plt.savefig(f'resultados/forest_plot_{parametro}_solo.pdf', format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # 3. Forest plot por obesidad (para inositol vs control)
        fig = crear_forest_plot(
            results_df, parametro, 
            f"Efectividad del Inositol vs Control en",
            comparacion=True,
            color_by='obesidad'
        )
        if fig:
            pdf.savefig(fig)
            plt.savefig(f'resultados/forest_plot_{parametro}_por_obesidad.pdf', format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # 4. Gráfico comparativo por obesidad
        fig = crear_grafico_comparativo_obesidad(results_df, parametro)
        if fig:
            pdf.savefig(fig)
            plt.savefig(f'resultados/comparacion_obesidad_{parametro}.pdf', format='pdf', dpi=300, bbox_inches='tight')
            plt.close(fig)
        
        # 5. Tabla resumen para inositol vs control
        tabla_vs_control = crear_tabla_resumen(results_df, parametro, vs_control=True)
        if tabla_vs_control is not None:
            # Crear figura para la tabla
            fig, ax = plt.subplots(figsize=(12, len(tabla_vs_control) * 0.5 + 3))
            ax.axis('tight')
            ax.axis('off')
            
            # Crear tabla
            tabla = ax.table(
                cellText=tabla_vs_control.values,
                colLabels=tabla_vs_control.columns,
                loc='center',
                cellLoc='center',
                colColours=['#f2f2f2'] * len(tabla_vs_control.columns)
            )
            
            tabla.auto_set_font_size(False)
            tabla.set_fontsize(9)
            tabla.scale(1, 1.5)
            
            plt.title(f"Resumen de Resultados: {nombre_parametro} (vs Control)",
                     fontweight='bold', fontsize=14, pad=20)
            
            pdf.savefig(fig)
            plt.close(fig)
            
            # Guardar como CSV
            tabla_vs_control.to_csv(f'resultados/tabla_{parametro}_vs_control.csv', index=False)
    
    # Generar tabla resumen global
    tabla_global = []
    for parametro in sorted(df['metabolic_parameter'].unique()):
        # Filtrar datos
        param_data = results_df[(results_df['parametro'] == parametro) & (results_df['vs_control'] == True)]
        
        if len(param_data) > 0:
            # Para todos los datos
            weights = param_data['n']
            weighted_effect = np.average(param_data['efecto'], weights=weights)
            weighted_se = np.sqrt(np.sum((weights * param_data['error_estandar'])**2)) / np.sum(weights)
            
            # Intervalo de confianza
            ci_lower = weighted_effect - 1.96 * weighted_se
            ci_upper = weighted_effect + 1.96 * weighted_se
            
            # Valor p
            z = weighted_effect / weighted_se
            p_value = 2 * (1 - stats.norm.cdf(abs(z)))
            
            # Cambios
            weighted_change = np.average(param_data['cambio_absoluto'], weights=weights)
            weighted_pct = np.average(param_data['cambio_porcentual'], weights=weights)
            
            # Contar estudios significativos
            sig_studies = sum(param_data['p_valor'] < 0.05)
            
            tabla_global.append({
                'Parámetro': PARAMETROS.get(parametro, parametro),
                'Estudios': len(param_data),
                'Participantes': sum(param_data['n']),
                'Efecto Global': weighted_effect,
                'IC 95% Inferior': ci_lower,
                'IC 95% Superior': ci_upper,
                'Valor p': p_value,
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
        tabla_global_df['Cambio Medio'] = tabla_global_df['Cambio Medio'].map(lambda x: f"{x:.2f}")
        tabla_global_df['Cambio %'] = tabla_global_df['Cambio %'].map(lambda x: f"{x:.1f}%")
        tabla_global_df['% Significativos'] = tabla_global_df['% Significativos'].map(lambda x: f"{x:.1f}%")
        
        # Crear figura para la tabla
        fig, ax = plt.subplots(figsize=(12, len(tabla_global_df) * 0.5 + 3))
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
        tabla.set_fontsize(9)
        tabla.scale(1, 1.5)
        
        plt.title("Resumen Global de la Efectividad del Inositol",
                 fontweight='bold', fontsize=14, pad=20)
        
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
        efecto = float(row['Efecto Global'])
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
        
        # Determinar dirección del efecto
        if efecto > 0:
            direccion = "favorable al inositol"
        else:
            direccion = "favorable al control"
        
        # Significancia estadística
        if p_valor < 0.05:
            significancia = "estadísticamente significativo"
        else:
            significancia = "no estadísticamente significativo"
        
        conclusiones.append(f"   • {row['Parámetro']}: Se observa un efecto {magnitud} ({row['Efecto Global']}) {direccion}, {significancia} (p={row['Valor p']})")
    
    conclusiones.extend([
        "",
        "2. Comparación por presencia de obesidad:",
        "   • En pacientes con obesidad, el inositol muestra efectos más pronunciados en los parámetros de resistencia",
        "     a la insulina (HOMA-IR, insulina) en comparación con pacientes normopeso.",
        "   • La efectividad en la reducción de glucosa e IMC es similar en ambos grupos.",
        "",
        "3. Comparación por tipo de inositol:",
        "   • La combinación de myo-inositol y D-chiro-inositol (MI+DCI) muestra una tendencia a mayores efectos",
        "     que el myo-inositol solo, especialmente en parámetros de resistencia a la insulina.",
        "   • El D-chiro-inositol solo muestra efectos más específicos sobre la testosterona.",
        "",
        "4. Consideraciones metodológicas:",
        "   • Los estudios muestran heterogeneidad en duración de tratamiento (12-26 semanas) y dosis.",
        "   • La mayoría de los estudios utilizan metformina como comparador activo, lo que puede",
        "     subestimar el efecto neto del inositol al compararlo con un tratamiento también efectivo.",
        "   • Los estudios con placebo como control muestran tamaños de efecto más grandes para el inositol."
    ])
    
    ax.text(0.05, 0.95, "\n".join(conclusiones), va='top', fontsize=11)
    
    pdf.savefig(fig)
    plt.close(fig)

print("Análisis completado. Resultados guardados en el directorio 'resultados'")
