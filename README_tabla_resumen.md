# Análisis de Efectos del Inositol - Tabla Resumen

## Descripción

El script `tabla_resumen_inositol.py` analiza los efectos del inositol en diferentes variables metabólicas, enfocándose únicamente en el análisis pre-post de los grupos de intervención que recibieron inositol.

## Características principales

- **Análisis pre-post**: Solo evalúa el cambio dentro de los grupos de intervención con inositol
- **No comparación con control**: Se enfoca en medir el efecto del inositol, no la diferencia respecto a grupos control
- **Tamaño del efecto**: Utiliza g de Hedges para calcular el tamaño del efecto
- **Variables incluidas**: Analiza todas las variables metabólicas disponibles, incluyendo volumen ovárico

## Uso

### Ejecución básica
```bash
python tabla_resumen_inositol.py
```

### Importar como módulo
```python
from tabla_resumen_inositol import analizar_datos_inositol, generar_tabla_resumen

# Analizar datos
df_resultados = analizar_datos_inositol('data/example_data.csv')

# Generar tabla
tabla_final, df_completo = generar_tabla_resumen(df_resultados)
```

## Salidas

El script genera:

1. **Tabla en pantalla**: Resultados formateados mostrados en consola
2. **CSV completo**: `output/tabla_resumen_efectos_inositol.csv` (datos originales)
3. **CSV formateado**: `output/tabla_resumen_efectos_inositol_formatted.csv` (tabla para presentación)

## Estructura de la tabla resumen

| Columna | Descripción |
|---------|-------------|
| Variable Analizada | Nombre de la variable metabólica |
| Nº Estudios | Número de estudios que evalúan esa variable |
| Tamaño Muestra | Tamaño de muestra total del grupo de intervención |
| Diferencia Medias Pre-Post | Diferencia de medias pre-post en el grupo de intervención |
| Tamaño Efecto (g Hedges) | Tamaño del efecto calculado con g de Hedges |
| Interpretación Efecto | Interpretación (Grande, Moderado, Pequeño, Sin efecto) |
| Intervalo Confianza 95% | Intervalo de confianza del 95% |

## Interpretación del tamaño del efecto

Según la convención de Cohen:
- **Sin efecto**: |g| < 0.2
- **Pequeño**: 0.2 ≤ |g| < 0.5
- **Moderado**: 0.5 ≤ |g| < 0.8
- **Grande**: |g| ≥ 0.8

## Tipos de tratamiento con inositol incluidos

- **MI**: Myo-inositol solo
- **DCI**: D-chiro-inositol solo
- **MI+DCI**: Combinación de myo-inositol y d-chiro-inositol
- **MI+MET**: Myo-inositol combinado con metformina

## Requisitos

- Python 3.7+
- pandas ≥ 2.0.0
- numpy ≥ 1.24.0
- scipy ≥ 1.10.0

## Notas importantes

- Los resultados están ordenados por tamaño del efecto (valor absoluto) de mayor a menor
- Solo se incluyen grupos de intervención, no grupos control
- El análisis se enfoca en el cambio pre-post dentro de cada grupo de intervención
- Las variables donde la reducción indica mejora (como IMC, glucosa, etc.) tienen el signo ajustado apropiadamente