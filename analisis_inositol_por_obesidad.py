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
    'axes.titlesize': 16,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 12,
    'figure.figsize': (16, 10)
})

# Cargar los datos
df = pd.read_csv('data/example_data.csv')

# Diccionario para traducir parámetros metabólicos
PARAMETROS = {
    'bmi_kg_m2': 'IMC (kg/m²)',
    'glucose_mg_dl': 'Glucosa (mg/dL)',
    'insulin_uU_ml': 'Insulina (μU/mL)',
    'homa_ir': 'HOMA-IR',
    'testosterone_ng_ml': 'Testosterona (ng/mL)',
    'ovarian_volume_ml': 'Volumen ovárico (mL)'
}

# Filtrar solo los parámetros que nos interesan
parametros_objetivo = ['bmi_kg_m2', 'glucose_mg_dl', 'insulin_uU_ml', 'homa_ir', 'testosterone_ng_ml', 'ovarian_volume_ml']

# Filtrar solo grupos de intervención con inositol
def contiene_inositol(tratamiento):
    if pd.isna(tratamiento):
        return False
    tratamiento = tratamiento.lower()
    return 'inositol' in tratamiento

# Función para calcular g de Hedges (pre-post dentro del mismo grupo)
def calcular_hedges_g(pre_mean, post_mean, pre_sd, post_sd, n):
    # Diferencia de medias (pre - post, valores negativos indican reducción/mejora)
    mean_diff = pre_mean - post_mean
    
    # Desviación estándar pooled para medidas repetidas
    pooled_sd = np.sqrt((pre_sd**2 + post_sd**2) / 2)
    
    if pooled_sd < 0.0001:
        pooled_sd = 0.0001
    
    # d de Cohen
    d = mean_diff / pooled_sd
    
    # Corrección para muestras pequeñas (g de Hedges)
    correction = 1 - (3 / (4 * n - 1))
    g = d * correction
    
    # Error estándar
    se = np.sqrt((2/n) + ((g**2) / (2*n)))
    
    # Intervalo de confianza del 95%
    ci_lower = g - 1.96 * se
    ci_upper = g + 1.96 * se
    
    # Valor p
    z = g / se
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return g, se, ci_lower, ci_upper, p_value, mean_diff

# Función para meta-análisis de efectos aleatorios
def calcular_efecto_combinado(efectos, errores):
    if len(efectos) < 2:
        return None, None, None, None, None
    
    # Calcular varianzas
    variances = errores ** 2
    
    # Pesos iniciales (inverso de la varianza)
    weights = 1 / variances
    
    # Calcular Q (heterogeneidad)
    Q = np.sum(weights * (efectos - np.average(efectos, weights=weights))**2)
    df = len(efectos) - 1
    
    # Calcular tau² (varianza entre estudios)
    if Q > df:
        tau_squared = (Q - df) / (np.sum(weights) - np.sum(weights**2) / np.sum(weights))
    else:
        tau_squared = 0
    
    # Nuevos pesos con tau²
    weights_random = 1 / (variances + tau_squared)
    
    # Efecto combinado
    combined_effect = np.average(efectos, weights=weights_random)
    se_combined = np.sqrt(1 / np.sum(weights_random))
    
    # Intervalo de confianza
    ci_lower = combined_effect - 1.96 * se_combined
    ci_upper = combined_effect + 1.96 * se_combined
    
    # Valor p
    z = combined_effect / se_combined
    p_value = 2 * (1 - stats.norm.cdf(abs(z)))
    
    return combined_effect, se_combined, ci_lower, ci_upper, p_value

# Preparar datos
datos_analisis = []

# Filtrar datos: solo intervenciones con inositol
df_inositol = df[(df['group_type'] == 'intervention') & 
                 (df['treatment_formulation'].apply(contiene_inositol)) &
                 (df['metabolic_parameter'].isin(parametros_objetivo))]

# Procesar cada fila
for _, row in df_inositol.iterrows():
    g, se, ci_lower, ci_upper, p_value, mean_diff = calcular_hedges_g(
        row['baseline_mean'], row['post_mean'],
        row['baseline_sd'], row['post_sd'],
        row['sample_size']
    )
    
    # Calcular cambio porcentual
    pct_change = ((row['post_mean'] - row['baseline_mean']) / row['baseline_mean']) * 100
    
    datos_analisis.append({
        'estudio': row['study_id'],
        'parametro': row['metabolic_parameter'],
        'parametro_es': PARAMETROS[row['metabolic_parameter']],
        'obesidad': 'Con sobrepeso/obesidad' if row['has_obesity'] == 'yes' else 'Peso normal',
        'n': row['sample_size'],
        'baseline_mean': row['baseline_mean'],
        'post_mean': row['post_mean'],
        'hedges_g': g,
        'se': se,
        'ci_lower': ci_lower,
        'ci_upper': ci_upper,
        'p_value': p_value,
        'mean_diff': mean_diff,
        'pct_change': pct_change
    })

# Convertir a DataFrame
df_analisis = pd.DataFrame(datos_analisis)

# Crear forest plot comparativo
def crear_forest_plot_comparativo():
    # Configurar la figura
    fig, axes = plt.subplots(1, len(parametros_objetivo), figsize=(20, 12), sharey=False)
    fig.suptitle('Efectos del Inositol: Comparación entre Pacientes con Peso Normal vs Sobrepeso/Obesidad', 
                 fontsize=18, fontweight='bold', y=0.95)
    
    colores = {'Peso normal': '#2E86AB', 'Con sobrepeso/obesidad': '#A23B72'}
    
    for idx, parametro in enumerate(parametros_objetivo):
        ax = axes[idx]
        
        # Filtrar datos para este parámetro
        param_data = df_analisis[df_analisis['parametro'] == parametro]
        
        if len(param_data) == 0:
            ax.text(0.5, 0.5, 'Sin datos', ha='center', va='center', transform=ax.transAxes)
            ax.set_title(PARAMETROS[parametro], fontweight='bold', fontsize=12)
            continue
        
        # Separar por grupo de obesidad
        normal_data = param_data[param_data['obesidad'] == 'Peso normal']
        obeso_data = param_data[param_data['obesidad'] == 'Con sobrepeso/obesidad']
        
        y_pos = 0
        all_effects = []
        all_se = []
        
        # Plotear datos de peso normal
        if len(normal_data) > 0:
            for _, row in normal_data.iterrows():
                # Línea de confianza (más gruesa)
                ax.plot([row['ci_lower'], row['ci_upper']], [y_pos, y_pos], 
                       color=colores['Peso normal'], linewidth=4, alpha=0.8)
                # Punto central (más grande)
                ax.scatter(row['hedges_g'], y_pos, color=colores['Peso normal'], 
                          s=150, zorder=3, edgecolor='white', linewidth=2)
                
                # Etiquetas con valores (fuente más pequeña)
                effect_text = f"g = {row['hedges_g']:.2f}"
                change_text = f"Δ = {row['mean_diff']:.1f} ({row['pct_change']:.1f}%)"
                
                ax.text(row['ci_upper'] + 0.05, y_pos + 0.1, effect_text, 
                       fontsize=10, va='center', ha='left', fontweight='bold')
                ax.text(row['ci_upper'] + 0.05, y_pos - 0.1, change_text, 
                       fontsize=10, va='center', ha='left', color='#555555')
                
                y_pos += 1
                all_effects.append(row['hedges_g'])
                all_se.append(row['se'])
            
            # Efecto combinado para peso normal
            if len(normal_data) > 1:
                combined_effect, se_combined, ci_lower_comb, ci_upper_comb, p_comb = calcular_efecto_combinado(
                    normal_data['hedges_g'].values, normal_data['se'].values
                )
                
                if combined_effect is not None:
                    # Línea combinada (más gruesa y distintiva)
                    ax.plot([ci_lower_comb, ci_upper_comb], [y_pos, y_pos], 
                           color=colores['Peso normal'], linewidth=6, alpha=1.0)
                    # Punto combinado (rombo)
                    ax.scatter(combined_effect, y_pos, color=colores['Peso normal'], 
                              s=200, zorder=4, marker='D', edgecolor='white', linewidth=2)
                    
                    # Etiqueta combinada
                    ax.text(ci_upper_comb + 0.05, y_pos, f"Combinado: g = {combined_effect:.2f}", 
                           fontsize=11, va='center', ha='left', fontweight='bold',
                           bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
                    
                    y_pos += 1.5
        
        # Plotear datos con sobrepeso/obesidad
        if len(obeso_data) > 0:
            for _, row in obeso_data.iterrows():
                # Línea de confianza (más gruesa)
                ax.plot([row['ci_lower'], row['ci_upper']], [y_pos, y_pos], 
                       color=colores['Con sobrepeso/obesidad'], linewidth=4, alpha=0.8)
                # Punto central (más grande)
                ax.scatter(row['hedges_g'], y_pos, color=colores['Con sobrepeso/obesidad'], 
                          s=150, zorder=3, edgecolor='white', linewidth=2)
                
                # Etiquetas con valores (fuente más pequeña)
                effect_text = f"g = {row['hedges_g']:.2f}"
                change_text = f"Δ = {row['mean_diff']:.1f} ({row['pct_change']:.1f}%)"
                
                ax.text(row['ci_upper'] + 0.05, y_pos + 0.1, effect_text, 
                       fontsize=10, va='center', ha='left', fontweight='bold')
                ax.text(row['ci_upper'] + 0.05, y_pos - 0.1, change_text, 
                       fontsize=10, va='center', ha='left', color='#555555')
                
                y_pos += 1
                all_effects.append(row['hedges_g'])
                all_se.append(row['se'])
            
            # Efecto combinado para obesidad
            if len(obeso_data) > 1:
                combined_effect, se_combined, ci_lower_comb, ci_upper_comb, p_comb = calcular_efecto_combinado(
                    obeso_data['hedges_g'].values, obeso_data['se'].values
                )
                
                if combined_effect is not None:
                    # Línea combinada (más gruesa y distintiva)
                    ax.plot([ci_lower_comb, ci_upper_comb], [y_pos, y_pos], 
                           color=colores['Con sobrepeso/obesidad'], linewidth=6, alpha=1.0)
                    # Punto combinado (rombo)
                    ax.scatter(combined_effect, y_pos, color=colores['Con sobrepeso/obesidad'], 
                              s=200, zorder=4, marker='D', edgecolor='white', linewidth=2)
                    
                    # Etiqueta combinada
                    ax.text(ci_upper_comb + 0.05, y_pos, f"Combinado: g = {combined_effect:.2f}", 
                           fontsize=11, va='center', ha='left', fontweight='bold',
                           bbox=dict(facecolor='white', alpha=0.8, boxstyle='round,pad=0.3'))
                    
                    y_pos += 1.5
        
        # Efecto global combinado
        if len(all_effects) > 1:
            global_effect, global_se, global_ci_lower, global_ci_upper, global_p = calcular_efecto_combinado(
                np.array(all_effects), np.array(all_se)
            )
            
            if global_effect is not None:
                # Línea global (negra, muy gruesa)
                ax.plot([global_ci_lower, global_ci_upper], [y_pos, y_pos], 
                       color='black', linewidth=8, alpha=0.9)
                # Punto global (estrella)
                ax.scatter(global_effect, y_pos, color='black', 
                          s=250, zorder=5, marker='*', edgecolor='white', linewidth=2)
                
                # Etiqueta global
                ax.text(global_ci_upper + 0.05, y_pos, f"Global: g = {global_effect:.2f}", 
                       fontsize=12, va='center', ha='left', fontweight='bold',
                       bbox=dict(facecolor='yellow', alpha=0.8, boxstyle='round,pad=0.3'))
        
        # Línea vertical en cero
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, linewidth=2)
        
        # Configurar ejes
        ax.set_title(PARAMETROS[parametro], fontweight='bold', fontsize=14)
        ax.set_xlabel('g de Hedges', fontsize=12, fontweight='bold')
        
        # Ajustar límites
        if len(param_data) > 0:
            all_ci = np.concatenate([param_data['ci_lower'].values, param_data['ci_upper'].values])
            margin = max(0.5, (max(all_ci) - min(all_ci)) * 0.2)
            ax.set_xlim(min(all_ci) - margin, max(all_ci) + margin)
            ax.set_ylim(-0.5, y_pos + 0.5)
        
        # Añadir interpretación
        ax.text(0, -0.3, "Sin efecto", ha='center', va='top', fontsize=10, 
               bbox=dict(facecolor='white', alpha=0.7, boxstyle='round,pad=0.2'))
        
        # Grilla sutil
        ax.grid(True, linestyle=':', alpha=0.3)
        ax.set_yticks([])
    
    # Leyenda global
    legend_elements = [
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colores['Peso normal'],
                   markersize=12, label='Peso normal', markeredgecolor='white', linewidth=2),
        plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=colores['Con sobrepeso/obesidad'],
                   markersize=12, label='Con sobrepeso/obesidad', markeredgecolor='white', linewidth=2),
        plt.Line2D([0], [0], marker='D', color='w', markerfacecolor='gray',
                   markersize=12, label='Efecto combinado por grupo', markeredgecolor='white', linewidth=2),
        plt.Line2D([0], [0], marker='*', color='w', markerfacecolor='black',
                   markersize=15, label='Efecto global', markeredgecolor='white', linewidth=2)
    ]
    
    fig.legend(handles=legend_elements, loc='lower center', ncol=4, fontsize=12,
               bbox_to_anchor=(0.5, 0.02), frameon=True, fancybox=True, shadow=True)
    
    # Nota explicativa
    fig.text(0.5, 0.08, 'Valores negativos en Δ y g indican reducción (mejoría). IC 95% mostrado.',
             ha='center', fontsize=11, style='italic')
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.15, top=0.9)
    
    return fig

# Crear el gráfico
fig = crear_forest_plot_comparativo()

# Guardar
os.makedirs('resultados', exist_ok=True)
plt.savefig('resultados/forest_plot_comparativo_obesidad.pdf', format='pdf', dpi=300, bbox_inches='tight')
plt.savefig('resultados/forest_plot_comparativo_obesidad.png', format='png', dpi=300, bbox_inches='tight')

# Crear tabla resumen
tabla_resumen = []

for parametro in parametros_objetivo:
    param_data = df_analisis[df_analisis['parametro'] == parametro]
    
    if len(param_data) > 0:
        # Por grupo de obesidad
        for grupo in ['Peso normal', 'Con sobrepeso/obesidad']:
            grupo_data = param_data[param_data['obesidad'] == grupo]
            
            if len(grupo_data) > 0:
                # Efectos ponderados por tamaño de muestra
                weights = grupo_data['n']
                efecto_ponderado = np.average(grupo_data['hedges_g'], weights=weights)
                cambio_ponderado = np.average(grupo_data['mean_diff'], weights=weights)
                pct_ponderado = np.average(grupo_data['pct_change'], weights=weights)
                
                # Efecto combinado
                if len(grupo_data) > 1:
                    combined_effect, _, ci_lower, ci_upper, p_value = calcular_efecto_combinado(
                        grupo_data['hedges_g'].values, grupo_data['se'].values
                    )
                else:
                    combined_effect = grupo_data['hedges_g'].iloc[0]
                    ci_lower = grupo_data['ci_lower'].iloc[0]
                    ci_upper = grupo_data['ci_upper'].iloc[0]
                    p_value = grupo_data['p_value'].iloc[0]
                
                tabla_resumen.append({
                    'Parámetro': PARAMETROS[parametro],
                    'Grupo': grupo,
                    'Estudios': len(grupo_data),
                    'Participantes': int(grupo_data['n'].sum()),
                    'g de Hedges': f"{combined_effect:.2f}",
                    'IC 95%': f"[{ci_lower:.2f}, {ci_upper:.2f}]",
                    'Cambio medio': f"{cambio_ponderado:.1f}",
                    'Cambio %': f"{pct_ponderado:.1f}%",
                    'Valor p': f"{p_value:.3f}",
                    'Significativo': 'Sí' if p_value < 0.05 else 'No'
                })

# Convertir a DataFrame y guardar
df_tabla = pd.DataFrame(tabla_resumen)
df_tabla.to_csv('resultados/tabla_resumen_comparativo.csv', index=False)

# Mostrar tabla
print("TABLA RESUMEN - EFECTOS DEL INOSITOL POR GRUPO DE PESO")
print("=" * 80)
print(df_tabla.to_string(index=False))

plt.show()

print(f"\nAnálisis completado. Archivos guardados en 'resultados/':")
print("- forest_plot_comparativo_obesidad.pdf")
print("- forest_plot_comparativo_obesidad.png") 
print("- tabla_resumen_comparativo.csv")
