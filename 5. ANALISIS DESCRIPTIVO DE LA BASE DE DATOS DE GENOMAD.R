##5. ANALISIS DESCRIPTIVO DE LA BASE DE DATOS DE GENOMAD
```r
library(readxl)
UGT1A1_gnomADfinal <- read_excel("Master en Bioinformática/AÑO 2 (2025-2026)/TFM/ANALISIS VCF_CHR37/obtencion de rs ID/metodo dsSNP/UGT1A1_gnomADfinal.xlsx")
View(UGT1A1_gnomADfinal)

# 5.1 Elegimos qué variantes son las que tiene mayor diversidad de frecuencia entre las poblaciones
#Seleccionar, por ejemplo, las 50 variantes con mayor variabilidad (desviación estándar)
varianza <- apply(datos_genomAD_freq, 1, sd)
top_variantes50 <- datos_genomAD_freq[order(varianza, decreasing = TRUE)[1:50], ]

#Guardar lista de datos de top_variantes50:

library(openxlsx)

write.xlsx(top_variantes50,
           file = "top_variantes50.xlsx",
           rowNames = FALSE)

#Convertir a matriz: 

matriz_heatmap <- as.matrix(top_variantes50)

#Poner los rsID como nombres de fila para que salgan en el gráfico:

rownames(matriz_heatmap) <- datos_genomAD_freq$rsID[order(varianza, decreasing = TRUE)[1:50]]

#Hacer el heatmap:

#Definir paleta de colores (de blanco a azul oscuro):

library(gplots)
colores <- colorRampPalette(c("white", "yellow", "orange", "red"))(100)

heatmap.2(matriz_heatmap,
          main = "Frecuencias de UGT1A1 por Población",
          trace = "none",          
          col = colores,           
          scale = "column",        
          margins = c(10, 8),      
          key = TRUE,              
          dendrogram = "both",     
          cexRow = 0.6,            
          cexCol = 0.8)            

#Vamos a coger las poblaciones específicas sin tener en cuenta la AF_global. Para coger solamente poblaciones especificas, sin tener en cuenta las globales:

poblaciones_especificas <- c(
  "AF_Europeos_No_Fin",
  "AS_Europa_Sur",
  "AF_Latino_Amr",
  "AF_Africano",
  "AS_Adiatico_Este",
  "AS_Judio_Ashkenazi",
  "AF_Finlandes"
)
media_poblacionesespecificas <- rowMeans(freq_num[, poblaciones_especificas], na.rm = TRUE)

#Seleccionar las poblaciones especificas para hacer heatmap

matriz_especifica <- as.matrix(freq_num[, poblaciones_especificas])
matriz_clean <- na.omit(matriz_especifica) #eliminamos valores NA

library(pheatmap)

#Utilizamos log 
mat_log <- log10(matriz_especifica + 1e-6)
mat_log <-  na.omit(mat_log)
pheatmap(mat_log,
         scale = "none",
         color = colorRampPalette(c("blue","white","red"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Frecuencias UGT1A1 (gnomAD, log10)")

#Como hay 3077 variantes, vamos a ver la variabilidad por poblacion. Nos quedamos con las 50 con más varianza. 

rango <- apply(mat_log, 1, function(x) max(x, na.rm=TRUE) - min(x, na.rm=TRUE))

top_var <- mat_log[order(-rango), ][1:50, ]

library(pheatmap)

pheatmap(top_var,
         scale = "row",   # CLAVE
         color = colorRampPalette(c("blue","white","red"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "50 variantes con mayor varianza UGT1A1 (gnomAD)")
```
# 5.2 Nos vamos a centrar en las variantes clínicas relevantes, vamos a ver qué variabilidad de las frecuencias hay entre poblaciones
#Vamos a construir la base de datos que queremos calcular las frecuencias, necesitamos filtarar por polaciones (sin las globales) y los rs que queremos
```r
vars_clinicas <- c(
  "rs887829",  
  "rs4148323", 
  "rs6742078", 
  "rs4124874", 
  "rs34983651")

poblaciones_especificas <- c(
  "AF_Europeos_No_Fin",
  "AS_Europa_Sur",
  "AF_Latino_Amr",
  "AF_Africano",
  "AS_Adiatico_Este",
  "AS_Judio_Ashkenazi",
  "AF_Finlandes"
)

#Sacar los rs de la base de datos junto a REF, y ALT
tabla_base <- UGT1A1_gnomADfinal[, c("rs", "REF", "ALT", poblaciones_especificas)]

#Construir tabla con estos rs y las poblaciones especificas solamente:

library(dplyr)

tabla_final_clinica <- UGT1A1_gnomADfinal %>%
  select(-AF_Global, -AF_PopMax, -AF_Otros, -AF_Homnres, -AF_Mujeres)

#Filtrar por rs
tabla_clinica_completa <- tabla_base[tabla_base$rs %in% vars_clinicas, ]

#A partir de esta base de datos, calculamos las frecuencias y comparaciones. Vamos a crear una base de datos de UGT1A1_genomAD pero quitandole las columnas de frecuencias globales
UGT1A1_filtred <- UGT1A1_gnomADfinal %>%
  select(-AF_Global, -AF_PopMax, -AF_Otros, -AF_Homnres, -AF_Mujeres)
```
# 5.3 Analisis de qué variantes son únicas en cada población 
```r
UGT1A1_filtred$diferencia_poblacional <- rs_diferencialestotales

#Vamos a filtrar los resultados por las variantes que tienen un 5% mas variables entre poblaciones
threshold <- quantile(UGT1A1_filtred$diferencia_poblacional, 0.95, na.rm = TRUE)

variantes_top <- UGT1A1_filtred[UGT1A1_filtred$diferencia_poblacional >= threshold, ]

library(pheatmap)

matriz_top_variabilidad <- as.matrix(variantes_top[, poblaciones_especificas])
rownames(matriz_top_variabilidad) <- variantes_top$rs

pheatmap(matriz_top_variabilidad,
         scale = "none",
         main = "Variantes con mayor 5% variabilidad interpoblacional")
```
#AHORA VAMOS A HACER LO MISMO PERO QUEDANDONOS SOLAMENTE CON 15-20 VARIANTES 
```r
variantes15 <- head(UGT1A1_filtred[order(-UGT1A1_filtred$diferencia_poblacional), ], 15)

#Graficamos y normalizamos:

UGT1A1_filtred[, poblaciones_especificas] <- 
  UGT1A1_filtred[, poblaciones_especificas] / 100

matriz_15 <- as.matrix(variantes15[, poblaciones_especificas])
rownames(matriz_15) <- variantes15$rs

pheatmap(matriz_15,
         scale = "row",
         color = colorRampPalette(c("white","white", "red"))(100),
         main = "Variantes de UGT1A1 con mayor variabilidad interpoblacional")
```