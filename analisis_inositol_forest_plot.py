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
    'figure.figsize': (16, 10)
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
    'bmi_kg_m2', 'glucose_mg_dl', 'insulin_uU_ml', 'homa_ir', 
    'ferriman_gallwey', 'ovarian_volume_ml', 'ovarian_volume_right_ml', 
    'ovarian_volume_left_ml', 'menstrual_days'
]

def calcular_efecto_combinado(data):
    """
    Calcula el efecto combinado usando meta-análisis de efectos aleatorios
    """
    if len(data) < 2:
        return None, None, None, None, None, None
    
    # Extraer valores
    effects = data['efecto'].values
    se_values = data['error_estandar'].values
    
    # Calcular pesos para efectos fijos
    weights_fixed = 1 / (se_values ** 2)
    
    # Efecto combinado de efectos fijos
    fixed_effect = np.sum(weights_fixed * effects) / np.sum(weights_fixed)
    
    # Calcular Q estadístico para heterogeneidad
    Q = np.sum(weights_fixed * (effects - fixed_effect) ** 2)
    df = len(effects) - 1
    
    # Calcular tau² (varianza entre estudios)
    if Q > df:
        C = np.sum(weights_fixed) - np.sum(weights_fixed ** 2) / np.sum(weights_fixed)
        tau_sq = (Q - df) / C
    else:
        tau_sq = 0
    
    # Pesos para efectos aleatorios
    weights_random = 1 / (se_values ** 2 + tau_sq)
    
    # Efecto combinado de efectos aleatorios
    combined_effect = np.sum(weights_random * effects) / np.sum(weights_random)
    se_combined = np.sqrt(1 / np.sum(weights_random))
    
    # Intervalo de confianza
    ci_lower = combined_effect - 1.96 * se_combined
    ci_upper = combined_effect + 1.96 * se_combined
    
    # Valor p
    z = combined_effect / se_combined
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    # I² estadístico
    I_squared = max(0, (Q - df) / Q * 100) if Q > 0 else 0
    
    return combined_effect, se_combined, ci_lower, ci_upper, p_value, I_squared

def crear_forest_plot_combinado(data, parametro, titulo, comparacion=True, filename=None, color_by='tipo_tratamiento'):
    """
    Crea un forest plot combinado que muestra tanto grupos de peso normal como sobrepeso/obesidad
    con mejor espaciado vertical y valores absolutos de Hedges' g para mayor claridad.
    """
    # Filtrar datos
    if comparacion:
        plot_data = data[(data['parametro'] == parametro) & (data['vs_control'] == True)]
    else:
        plot_data = data[(data['parametro'] == parametro) & (data['vs_control'] == False)]
    
    if len(plot_data) < 2:
        print(f"Datos insuficientes para {parametro}, {titulo}")
        return None
    
    # Separar datos por grupo de obesidad
    data_normal = plot_data[plot_data['obesidad'] == 'No'].copy()
    data_obeso = plot_data[plot_data['obesidad'] == 'Sí'].copy()
    
    # Verificar que tenemos datos de ambos grupos
    if len(data_normal) == 0 and len(data_obeso) == 0:
        print(f"Sin datos de grupos de obesidad para {parametro}")
        return None
    
    # Ordenar por tamaño del efecto absoluto para mejor visualización
    if len(data_normal) > 0:
        data_normal = data_normal.sort_values('efecto', key=abs, ascending=False)
    if len(data_obeso) > 0:
        data_obeso = data_obeso.sort_values('efecto', key=abs, ascending=False)
    
    # Calcular dimensiones del gráfico
    total_studies = len(data_normal) + len(data_obeso)
    group_spacing = 1.5  # Espaciado adicional entre grupos
    fig_height = max(10, total_studies * 0.8 + group_spacing * 2)
    
    # Crear figura
    fig, ax = plt.subplots(figsize=(16, fig_height))
    
    # Posiciones Y con mejor espaciado
    y_positions = []
    labels = []
    plot_data_combined = []
    
    current_y = 0
    
    # Procesar grupo de peso normal
    if len(data_normal) > 0:
        # Añadir encabezado del grupo
        ax.text(-0.5, current_y + 0.3, "Peso normal", fontsize=12, fontweight='bold', 
                bbox=dict(facecolor='lightblue', alpha=0.3, boxstyle='round,pad=0.3'))
        
        for i, (_, row) in enumerate(data_normal.iterrows()):
            y_positions.append(current_y)
            if comparacion:
                label = f"{row['estudio']} ({row['tipo_tratamiento']} vs {row['tipo_control']})"
            else:
                label = f"{row['estudio']} ({row['tipo_tratamiento']})"
            labels.append(label)
            plot_data_combined.append(row)
            current_y -= 1
        
        current_y -= group_spacing  # Espaciado entre grupos
    
    # Procesar grupo de sobrepeso/obesidad
    if len(data_obeso) > 0:
        # Añadir encabezado del grupo
        ax.text(-0.5, current_y + 0.3, "Sobrepeso/Obesidad", fontsize=12, fontweight='bold',
                bbox=dict(facecolor='lightcoral', alpha=0.3, boxstyle='round,pad=0.3'))
        
        for i, (_, row) in enumerate(data_obeso.iterrows()):
            y_positions.append(current_y)
            if comparacion:
                label = f"{row['estudio']} ({row['tipo_tratamiento']} vs {row['tipo_control']})"
            else:
                label = f"{row['estudio']} ({row['tipo_tratamiento']})"
            labels.append(label)
            plot_data_combined.append(row)
            current_y -= 1
    
    # Convertir a DataFrame para compatibilidad
    plot_data_combined = pd.DataFrame(plot_data_combined)
    
    # Determinar colores según lo especificado
    if color_by == 'tipo_tratamiento':
        colors = [COLORS.get(row['tipo_tratamiento'], 'gray') for _, row in plot_data_combined.iterrows()]
    elif color_by == 'obesidad':
        colors = [COLORS['MI'] if row['obesidad'] == 'Sí' else COLORS['DCI'] for _, row in plot_data_combined.iterrows()]
    else:
        colors = [COLORS['MI']] * len(plot_data_combined)
    
    # Graficar puntos y líneas de error
    for i, (_, row) in enumerate(plot_data_combined.iterrows()):
        y_pos = y_positions[i]
        
        # Línea para CI
        ax.plot([row['ci_inferior'], row['ci_superior']], [y_pos, y_pos], 
                color=colors[i], linewidth=2.5, alpha=0.8)
        # Punto para efecto
        ax.scatter(row['efecto'], y_pos, color=colors[i], s=120, zorder=3, 
                  edgecolor='black', linewidth=1)
    
    # Línea vertical en cero
    ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, linewidth=1.5)
    
    # Añadir etiquetas y detalles con valores absolutos para Hedges' g
    for i, (_, row) in enumerate(plot_data_combined.iterrows()):
        y_pos = y_positions[i]
        
        # Usar valor absoluto para mostrar el tamaño del efecto
        abs_effect = abs(row['efecto'])
        abs_ci_lower = abs(row['ci_inferior']) if row['ci_inferior'] * row['efecto'] >= 0 else abs(row['ci_inferior'])
        abs_ci_superior = abs(row['ci_superior']) if row['ci_superior'] * row['efecto'] >= 0 else abs(row['ci_superior'])
        
        # Texto con valor de efecto absoluto e intervalo de confianza
        effect_text = f"|g| = {abs_effect:.2f} [CI: {row['ci_inferior']:.2f}, {row['ci_superior']:.2f}]"
        
        # Texto con cambios absolutos y porcentuales (preservar signos)
        if comparacion:
            change_text = (f"Δ = {row['cambio_absoluto']:.2f} vs {row['cambio_absoluto_control']:.2f} " +
                          f"({row['cambio_porcentual']:.1f}% vs {row['cambio_porcentual_control']:.1f}%)")
        else:
            change_text = f"Δ = {row['cambio_absoluto']:.2f} ({row['cambio_porcentual']:.1f}%)"
        
        # Posición del texto según valor del efecto
        if row['efecto'] > 0:
            text_pos = row['ci_superior'] + 0.15
            ha = 'left'
        else:
            text_pos = row['ci_inferior'] - 0.15
            ha = 'right'
        
        # Ajustar posición vertical para evitar solapamiento
        ax.text(text_pos, y_pos + 0.1, effect_text, va='center', ha=ha, 
                fontsize=10, fontweight='bold')
        ax.text(text_pos, y_pos - 0.2, change_text, va='center', ha=ha, 
                fontsize=9, color='#333333')
        
        # Añadir valor p
        p_text = f"p = {row['p_valor']:.4f}"
        if row['p_valor'] < 0.05:
            p_text += "*"
            ax.axhspan(y_pos - 0.35, y_pos + 0.35, alpha=0.1, color='green')
        ax.text(text_pos, y_pos - 0.45, p_text, va='center', ha=ha, 
                fontsize=9, color='#333333')
    
    # Calcular y añadir el efecto combinado
    combined_effect, se_combined, ci_lower, ci_upper, p_value, I_squared = calcular_efecto_combinado(plot_data_combined)
    
    if combined_effect is not None:
        # Posición para el efecto combinado
        combined_y = min(y_positions) - 2
        
        # Añadir línea para el efecto combinado
        ax.axvline(x=combined_effect, color=COLORS['Combinado'], linestyle='-', linewidth=2.5, alpha=0.9)
        
        # Añadir región sombreada para el intervalo de confianza del efecto combinado
        ax.axvspan(ci_lower, ci_upper, alpha=0.2, color=COLORS['Combinado'])
        
        # Etiqueta para el efecto combinado (usar valor absoluto para consistencia)
        abs_combined = abs(combined_effect)
        combined_text = (f"Efecto Combinado: |g| = {abs_combined:.2f} [CI: {ci_lower:.2f}, {ci_upper:.2f}], " +
                        f"p = {p_value:.4f}, I² = {I_squared:.1f}%")
        
        ax.text(combined_effect, combined_y, combined_text, ha='center', va='center', 
                fontsize=11, fontweight='bold', 
                bbox=dict(facecolor='white', alpha=0.8, edgecolor=COLORS['Combinado'], 
                         boxstyle='round,pad=0.5'))
    
    # Configurar ejes
    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlabel('Tamaño del Efecto (g de Hedges)', fontsize=12, fontweight='bold')
    
    # Ajustar título
    param_name = PARAMETROS.get(parametro, parametro)
    ax.set_title(f"{titulo}: {param_name}", fontsize=16, fontweight='bold')
    
    # Añadir interpretación según el tipo de parámetro
    interpretation_y = min(y_positions) - 4
    if parametro in MEJORA_CON_REDUCCION:
        if comparacion:
            ax.text(-0.5, interpretation_y, "Favorece al control", ha='center', fontsize=10, 
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
            ax.text(0.5, interpretation_y, "Favorece al inositol", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        else:
            ax.text(-0.5, interpretation_y, "Empeoramiento", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
            ax.text(0.5, interpretation_y, "Mejoría", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    else:
        # Para parámetros donde el aumento es mejora
        if comparacion:
            ax.text(-0.5, interpretation_y, "Favorece al inositol", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
            ax.text(0.5, interpretation_y, "Favorece al control", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
        else:
            ax.text(-0.5, interpretation_y, "Mejoría", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
            ax.text(0.5, interpretation_y, "Empeoramiento", ha='center', fontsize=10,
                   bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.3'))
    
    # Leyenda para tipo de tratamiento
    if color_by == 'tipo_tratamiento':
        legend_elements = []
        for treatment_type, color in COLORS.items():
            if treatment_type != 'Combinado' and any(row['tipo_tratamiento'] == treatment_type for _, row in plot_data_combined.iterrows()):
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
        plot_data_combined['ci_inferior'].values,
        plot_data_combined['ci_superior'].values,
        [plot_data_combined['efecto'].min() - 0.3, plot_data_combined['efecto'].max() + 0.3]
    ])
    
    if combined_effect is not None:
        all_limits = np.append(all_limits, [ci_lower - 0.2, ci_upper + 0.2])
    
    margin = max(1.0, abs(all_limits.min()) * 0.25, abs(all_limits.max()) * 0.25)
    ax.set_xlim(min(all_limits) - margin, max(all_limits) + margin)
    ax.set_ylim(min(y_positions) - 5, max(y_positions) + 1)
    
    # Añadir cuadrícula sutil
    ax.grid(True, linestyle='--', alpha=0.3)
    
    plt.tight_layout()
    
    if filename:
        plt.savefig(filename, format='pdf', bbox_inches='tight', dpi=300)
    
    return fig

def crear_tabla_resumen_con_abs(data, parametro, vs_control=True):
    """
    Crear tabla resumen con valores absolutos de Hedges' g pero preservando signos en cambios
    """
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
                
                # Usar valores absolutos para el tamaño del efecto
                weighted_effect_abs = np.average(np.abs(obes_data['efecto']), weights=weights)
                weighted_change_abs = np.average(obes_data['cambio_absoluto'], weights=weights)  # Preservar signo
                weighted_change_pct = np.average(obes_data['cambio_porcentual'], weights=weights)  # Preservar signo
                
                # Contar estudios significativos
                sig_studies = sum(obes_data['p_valor'] < 0.05)
                pct_sig = (sig_studies / len(obes_data)) * 100
                
                summary_data.append({
                    'Parámetro': PARAMETROS.get(parametro, parametro),
                    'Tipo Tratamiento': tipo,
                    'Obesidad': obesidad,
                    'Estudios': len(obes_data),
                    'Participantes': sum(weights),
                    'Tamaño Efecto (|g|)': round(weighted_effect_abs, 2),  # Valor absoluto
                    'Cambio Absoluto': round(weighted_change_abs, 2),  # Preservar signo
                    'Cambio %': round(weighted_change_pct, 1),  # Preservar signo
                    'Estudios Significativos': sig_studies,
                    '% Significativos': round(pct_sig, 1)
                })
    
    # Añadir resumen general por obesidad
    for obesidad in ['Sí', 'No']:
        obes_data = param_data[param_data['obesidad'] == obesidad]
        
        if len(obes_data) > 0:
            weights = obes_data['n']
            weighted_effect_abs = np.average(np.abs(obes_data['efecto']), weights=weights)
            weighted_change_abs = np.average(obes_data['cambio_absoluto'], weights=weights)
            weighted_change_pct = np.average(obes_data['cambio_porcentual'], weights=weights)
            
            sig_studies = sum(obes_data['p_valor'] < 0.05)
            pct_sig = (sig_studies / len(obes_data)) * 100
            
            summary_data.append({
                'Parámetro': PARAMETROS.get(parametro, parametro),
                'Tipo Tratamiento': 'TODOS',
                'Obesidad': obesidad,
                'Estudios': len(obes_data),
                'Participantes': sum(weights),
                'Tamaño Efecto (|g|)': round(weighted_effect_abs, 2),
                'Cambio Absoluto': round(weighted_change_abs, 2),
                'Cambio %': round(weighted_change_pct, 1),
                'Estudios Significativos': sig_studies,
                '% Significativos': round(pct_sig, 1)
            })
    
    return pd.DataFrame(summary_data)