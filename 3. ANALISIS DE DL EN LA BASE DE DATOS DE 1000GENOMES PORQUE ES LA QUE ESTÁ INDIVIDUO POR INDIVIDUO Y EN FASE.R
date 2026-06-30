# 3. ANALISIS DE DL EN LA BASE DE DATOS DE 1000GENOMES PORQUE ES LA QUE ESTÁ INDIVIDUO POR INDIVIDUO Y EN FASE
# 3.1 Quedarnos solamente con variantes relevantes 
#en bash, archivo aparte.

#Una vez tenemos los archivos de las variantes relevantes por población, hacemos el análisis de desequilibrio de ligamiento.
#Instalacion de paquetes necesarios:
```r
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "gdsfmt"))

#carga de librerias 
library(SNPRelate)
library(gdsfmt)
setwd(utils::choose.dir())
```
# 3.1.1 Hacemos el análisis de DL para variantes conn MAF > 0.25%
#Abrir el archivo
```r
genofile <- snpgdsOpen("EUR.gds")

#SELECCIONAR SNPs con filtro estricto MAF > 0.25

snp_seleccionados <- snpgdsSelectSNP(genofile, maf=0.25, missing.rate=0.05)

E#XTRAER LOS NOMBRES REALES (rsID)

#Primero sacamos la lista completa de rsIDs que hay en el archivo, y luego sacamos la lista de todos los IDs numéricos
lista_completa_rsids <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
lista_completa_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))

#Filtramos los rsIDs para que coincidan EXACTAMENTE con los SNPs que pasaron el filtro 0.25
nombres_para_grafico <- lista_completa_rsids[lista_completa_ids %in% snp_seleccionados]
```

#CALCULAR MATRIZ DE DL:
```r
ld_obj <- snpgdsLDMat(genofile, snp.id=snp_seleccionados, method="r", slide=0)
matrix_final <- ld_obj$LD^2

#ASIGNAR NOMBRES
colnames(matrix_final) <- nombres_para_grafico
rownames(matrix_final) <- nombres_para_grafico

#CERRAR
snpgdsClose(genofile)
#graficar
corrplot(matrix_final, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.cex = 0.7,      # Tamaño de letra de los rsID
         addgrid.col = "gray", 
         main = "DL en variantes del gen UGT1A1 en población EUROPEA con MAF > 0.25",
         mar = c(0,0,2,0))

#Ver cuales son las variantes con MAF >0.25 en Europa y comparamos con el resto de poblaciones
genofile <- snpgdsOpen("EUR.gds")
```

#Filtro MAF > 0.25
```r
snp_sel <- snpgdsSelectSNP(genofile, maf=0.25, missing.rate=0.05)

#EXTRAER DATOS (Forzamos la lectura de rsID y de Posición)
ids_numericos <- read.gdsn(index.gdsn(genofile, "snp.id"))
rsids_vcf <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
posiciones <- read.gdsn(index.gdsn(genofile, "snp.position"))

#CREAMOS ETIQUETAS ÚNICAS

#Si el rsID es vacío o ".", usamos la posición genómica (Chr2:pos)
etiquetas <- ifelse(rsids_vcf == "" | rsids_vcf == ".", 
                    paste0("2:", posiciones), 
                    rsids_vcf)

#FILTRAR ETIQUETAS PARA EL GRÁFICO
nombres_finales <- etiquetas[ids_numericos %in% snp_sel]

#CALCULAR DL
ld_obj <- snpgdsLDMat(genofile, snp.id=snp_sel, method="r", slide=0)
matrix_eur <- ld_obj$LD^2
colnames(matrix_eur) <- nombres_finales
rownames(matrix_eur) <- nombres_finales

snpgdsClose(genofile)
```
#Graficamos
```r
corrplot(matrix_eur, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.cex = 0.7,      # Tamaño de letra de los rsID
         addgrid.col = "gray", 
         main = "DL en variantes del gen UGT1A1 en población EUROPEA con MAF > 0.25",
         mar = c(0,0,2,0))

```
#Hacemos lo mismo en poblacion Africana
```r
snpgdsVCF2GDS("UGT1A1_1000G_filtered_AFR.vcf.gz", "AFR.gds", method="biallelic.only")

genofile_afr <- snpgdsOpen("AFR.gds")

ld_obj_afr <- snpgdsLDMat(genofile_afr, snp.id=snp_sel, method="r", slide=0)
matrix_afr <- ld_obj_afr$LD^2

#Ponemos los mismos nombres para que la comparación sea visualmente idéntica:

colnames(matrix_afr) <- nombres_finales
rownames(matrix_afr) <- nombres_finales

snpgdsClose(genofile_afr)

#VISUALIZAR AMBOS GRAFICOS 
par(mfrow=c(1,2)) # Divide el gráfico en 2 columnas

library(corrplot)

#Gráfico EUR
corrplot(matrix_eur, method="color", type="upper", tl.cex=0.6, main="EUR (MAF > 0.2)")

#Gráfico AFR
corrplot(matrix_afr, method="color", type="upper", tl.cex=0.6, main="AFR (Mismos SNPs)")

#Convertimos las 3 poblaciones restantes:
snpgdsVCF2GDS("UGT1A1_1000G_filtered_EAS.vcf.gz", "EAS.gds", method="biallelic.only")
snpgdsVCF2GDS("UGT1A1_1000G_filtered_AMR.vcf.gz", "AMR.gds", method="biallelic.only")
snpgdsVCF2GDS("UGT1A1_1000G_filtered_SAS.vcf.gz", "SAS.gds", method="biallelic.only")

#Ponemos un codigo comparaito de todas las poblaciones

obtener_matriz <- function(archivo_gds, snps, etiquetas) {
  g <- snpgdsOpen(archivo_gds)
  ld <- snpgdsLDMat(g, snp.id=snps, method="r", slide=0)
  m <- ld$LD^2
  colnames(m) <- etiquetas
  rownames(m) <- etiquetas
  snpgdsClose(g)
  return(m)
}

#Calculamos las matrices restantes

matrix_afr <- obtener_matriz("AFR.gds", snp_sel, nombres_finales)
matrix_eas <- obtener_matriz("EAS.gds", snp_sel, nombres_finales)
matrix_amr <- obtener_matriz("AMR.gds", snp_sel, nombres_finales)
matrix_sas <- obtener_matriz("SAS.gds", snp_sel, nombres_finales)

#Configuramos el panel (2 filas, 3 columnas)

par(mfrow = c(2, 3), mar = c(2, 2, 4, 2))

#Función para graficar mostrando los nombres de las variantes

dibujar_con_nombres <- function(mat, titulo) {
  corrplot(mat, 
           method = "color", 
           type = "upper", 
           tl.col = "black",   # Color del texto
           tl.cex = 0.5,      # Tamaño de la letra (ajústalo si son muchos SNPs)
           tl.srt = 45,       # Rotación de las etiquetas a 45 grados
           tl.pos = "lt",     # 'lt' coloca las etiquetas arriba y a la izquierda
           main = titulo, 
           cl.lim = c(0, 1))
}
```
#Dibujamos las 5 poblaciones
```r
dibujar_con_nombres(matrix_eur, "EUR (MAF > 0.25)")
dibujar_con_nombres(matrix_afr, "AFR")
dibujar_con_nombres(matrix_eas, "EAS")
dibujar_con_nombres(matrix_amr, "AMR")
dibujar_con_nombres(matrix_sas, "SAS")

#Volver a 1 solo gráfico

par(mfrow = c(1, 1))
```
# 3.1.2 Creamos un PDF de los gráficos para que tenga una mejor calidad y no se corten los ejes:
```r
#1)Definimos el nombre del archivo y el tamaño (grande para que quepa todo)
pdf("Comparativa_LD_UGT1A1_Completa.pdf", width = 15, height = 10)

#2)Configuramos el panel (2 filas, 3 columnas) con márgenes más amplios
par(mfrow = c(2, 3), mar = c(4, 4, 6, 2))

#3)Función optimizada para que no se corten los rsID
dibujar_tfm <- function(mat, titulo) {
  corrplot(mat, 
           method = "color", 
           type = "upper", 
           tl.col = "black", 
           tl.cex = 0.8,      # Texto más legible
           tl.srt = 90,       # Rotación vertical para que ocupen menos ancho
           tl.pos = "lt",     # Etiquetas arriba y a la izquierda
           diag = FALSE,      # Quitar la diagonal para limpiar el gráfico
           main = titulo, 
           mar = c(0, 0, 4, 0), # Espacio para el título
           cl.lim = c(0, 1))
}

#4)Dibujamos cada una

dibujar_tfm(matrix_eur, "EUROPE (EUR)")
dibujar_tfm(matrix_afr, "AFRICA (AFR)")
dibujar_tfm(matrix_eas, "EAST ASIA (EAS)")
dibujar_tfm(matrix_amr, "ADMIXED AMERICAS (AMR)")
dibujar_tfm(matrix_sas, "SOUTH ASIA (SAS)")

#5)Cerramos el PDF

dev.off()
```