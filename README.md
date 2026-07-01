# TFE_Carla-Gonzalez-de-la-Cruz

Este repositorio contiene el flujo de trabajo bioinformático desarrollado para el trabajo titulado: 
**“Relevancia clínica poblacional de variantes del gen UGT1A1 asociadas a toxicidad por irinotecán”**

El objetivo principal es caracterizar la variabilidad genética del gen *UGT1A1* en distintas poblaciones humanas, estudiar los patrones de desequilibrio de ligamiento (LD) entre variantes farmacogenéticas relevantes y evaluar su posible impacto en la respuesta clínica al tratamiento con irinotecán.Se utilizan los datos genómicos procedentes de bases de datos públicas: 1000genomes y genomAD v2.1.1. 

## Fundamentación biológica

El gen *UGT1A1* codifica una enzima clave en la glucuronidación del metabolito activo del irinotecán (SN-38). Variantes genéticas en este gen pueden modificar la actividad enzimática y aumentar el riesgo de desarrollar efectos adversos graves, especialmente neutropenia severa y toxicidad gastrointestinal.
La frecuencia de estas variantes y su patrón de herencia conjunta (desequilibrio de ligamiento) difieren entre poblaciones, lo que puede tener implicaciones relevantes en farmacogenética y medicina personalizada.

## Objetivos del estudio

Este análisis pretende:

1)	Consultar bases de datos públicas para procesamiento de muestras de secuenciación del gen UGT1A1.
2)	Caracterizar las variantes encontradas en los datos de secuenciación del gen UGT1A1 alineados con el cromosoma de referencia hg19. 
3)	Analizar el desequilibrio de ligamiento de las diferentes variantes genéticas en las diferentes poblaciones. 
4)	Describir frecuencias poblacionales de variantes del gen UGT1A1 asociadas a toxicidad por Irinotecan y cómo se distribuyen en las diferentes poblaciones. 


## La metodología empleada es la siguiente: 

### 1. Obtención y filtrado de variantes

Se descargan los archivos VCF correspondientes al cromosoma 2 de las dos bases públicas utilizadas. Posteriormente se filtra la región genómica correspondiente al gen *UGT1A1* utilizando coordenadas de referencia GRCh37. Para ello, las herramientas utilizadas han sido: 
* bcftools
* samtools
* tabix

### 2. Análisis de frecuencias poblacionales

Se calculan las frecuencias alélicas de cada variante para las principales poblaciones representadas en 1000 Genomes:

* AFR (África)
* EUR (Europa)
* EAS (Asia Oriental)
* SAS (Asia del Sur)
* AMR (América)


### 3. Asignación de rsID y normalización de variantes

Las variantes obtenidas se normalizan utilizando el genoma de referencia hg19/GRCh37 y se cruzan con la base de datos dbSNP para incorporar sus correspondientes identificadores rsID.

### 4. Análisis de desequilibrio de ligamiento (LD)

El análisis de LD se realiza utilizando los genotipos individuales y faseados de 1000 Genomes. Se seleccionan variantes comunes mediante filtros de:

* Frecuencia alélica menor (MAF) > 0.25
* Frecuencia alélica menor (MAF) > 0.1
* Frecuencia de 4 variantes relevantes en clínica: rs8175347 (UGT1A1*28), rs4148323 (UGT1A1*6), rs887829, rs6717546
Las matrices de LD (r²) se calculan mediante el paquete SNPRelate en R. El objetivo es identificar bloques de variantes heredadas conjuntamente y estudiar las diferencias poblacionales en la estructura genética de UGT1A1.
Además, se hacen los Test de Mantel para ver patrones de genética de poblaciones diferenciadores. 

### 5. Estudio de diversidad genética utilizando gnomAD

La base de datos gnomAD se utiliza para:

* Analizar frecuencias poblacionales ampliadas.
* Identificar variantes con elevada variabilidad entre poblaciones.
* Detectar variantes enriquecidas en grupos específicos.
* Visualizar patrones poblacionales mediante mapas de calor.

Entre los análisis realizados se incluyen:

* Selección de las variantes con mayor variabilidad poblacional.
* Identificación de variantes potencialmente específicas de población.
* Evaluación de variantes farmacogenéticas relevantes.

### 6. Anotación funcional de variantes

Las variantes más relevantes se anotan mediante la API de Ensembl VEP.

La información recuperada incluye:

* Posición genómica.
* Consecuencia molecular.
* Nomenclatura HGVS.
* Transcritos afectados.
* Cambios proteicos.
* Identificadores UniProt.
* Clasificación funcional.

Los resultados se exportan a archivos Excel para facilitar su interpretación y revisión. Las herramientas Linux utilizadas son: bcftools, samtools, tabix, awk, join, wget

### Paquetes de R

```r
SNPRelate
gdsfmt
corrplot
vegan
readxl
openxlsx
gplots
pheatmap
dplyr
jsonlite
httr2
```
## Resultados principales

* Tablas de frecuencias alélicas por población.
* Variantes anotadas con rsID.
* Matrices de desequilibrio de ligamiento.
* Heatmap comparativos entre poblaciones con las variantes con mayor variabilidad interétnica.
* Anotaciones funcionales de variantes relevantes.
* Interpretación de resultados.
---
## Relevancia clínica
Este estudio proporciona una caracterización detallada de la variabilidad genética y de la estructura de desequilibrio de ligamiento del gen UGT1A1 en diferentes poblaciones humanas, poniendo especial interés en las variantes relevantes en clínica. 
Los resultados contribuyen a mejorar la interpretación de biomarcadores farmacogenéticos asociados a la toxicidad por irinotecán y aportan información útil para futuras estrategias de medicina personalizada y farmacogenómica de precisión.



