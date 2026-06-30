## 4. ANALISIS DE DESEQUILIBRIO DE LIGAMIENTO Y FRECUENCIAS DE VARIANTES INMPORTANTES EN CLINICA: rs8175347, rs4148323, rs887829 y rs6717546. 
#Filtrar nel archivo original vcf de 1000Genomes para cada poblacion, sin tener en cuenta un MAF: 
```r
bcftools view -S EUR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.EURcomplet.vcf.gz
bcftools view -S AFR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.AFRcomplet.vcf.gz
bcftools view -S AMR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.AMRcomplet.vcf.gz
bcftools view -S EAS.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.EAScomplet.vcf.gz
bcftools view -S SAS.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.SAScomplet.vcf.gz

#Indexar:

for pop in EUR AFR AMR EAS SAS
do
bcftools index UGT1A1.${pop}complet.vcf.gz
done
```
# 4.1 Calcular el DL de cada poblacion
#Convertir todos los VCF a GDS para que los trabaje R. 
```r
library(SNPRelate)

snpgdsVCF2GDS("UGT1A1.EURcomplet.vcf.gz", "EUR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.AFRcomplet.vcf.gz", "AFR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.AMRcomplet.vcf.gz", "AMR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.EAScomplet.vcf.gz", "EAS.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.SAScomplet.vcf.gz", "SAS.gds", method="copy.num.of.ref")

#Empezamos con poblacion EUR:

gen_eur <- snpgdsOpen("EUR.gds")

#Ver que las posiciones estan dentro del archivo

pos <- read.gdsn(index.gdsn(gen_eur, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_eur, "snp.id"))

#Ver si estan las posiciones que me interesan

mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel <- ids[pos %in% mis_posiciones]
ids_sel
length(ids_sel)

#Calculo de LD
ld <- snpgdsLDMat(gen_eur,
                  snp.id = ids_sel,
                  method = "r",
                  slide = 0)

mat_ld <- ld$LD^2

#Añadir etiquetas

labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld) <- rownames(mat_ld) <- labels

#Heatmap

library(corrplot)

corrplot(mat_ld,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Europea",
         mar = c(0,0,3,0))

snpgdsClose(gen_eur)
```
#AHORA VAMOS A HACERLO PARA EL RESTO DE POBLACIONES, Y PARA NO REPETIR CÓDIGO, AJUSTAMOS: 
```r
gen_afr <- snpgdsOpen("AFR.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_afr, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_afr, "snp.id"))

#Ver si estan las posiciones que me interesan:

mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel_afr <- ids[pos %in% mis_posiciones]
ids_sel_afr
length(ids_sel)

#Calculo de LD

ld_afr <- snpgdsLDMat(gen_afr,
                      snp.id = ids_sel_afr,
                      method = "r",
                      slide = 0)

mat_ld_afr <- ld_afr$LD^2

#Añadir etiquetas

labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_afr) <- rownames(mat_ld_afr) <- labels

#Heatmap

library(corrplot)

corrplot(mat_ld_afr,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Africana",
         mar = c(0,0,3,0))

snpgdsClose(gen_afr)
```
#POBLACION AMERICANA
```r
gen_amr <- snpgdsOpen("AMR.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_amr, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_amr, "snp.id"))

#Ver si estan las posiciones que me interesan
mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel_amr <- ids[pos %in% mis_posiciones]
ids_sel_amr
length(ids_sel_amr)
#Calculo de LD
ld_amr <- snpgdsLDMat(gen_amr,
                      snp.id = ids_sel_amr,
                      method = "r",
                      slide = 0)

mat_ld_amr <- ld_amr$LD^2

#Añadir etiquetas

labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_amr) <- rownames(mat_ld_amr) <- labels

#Heatmap

library(corrplot)

corrplot(mat_ld_amr,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Americana",
         mar = c(0,0,3,0))

snpgdsClose(gen_amr)
```
#POBLACION SAS
```r
gen_sas <- snpgdsOpen("SAS.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_sas, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_sas, "snp.id"))

#Ver si estan las posiciones que me interesan

mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel_sas <- ids[pos %in% mis_posiciones]
ids_sel_sas
length(ids_sel_sas)
#Calculo de LD
ld_sas <- snpgdsLDMat(gen_sas,
                      snp.id = ids_sel_sas,
                      method = "r",
                      slide = 0)

mat_ld_sas <- ld_sas$LD^2

#Añadir etiquetas

labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_sas) <- rownames(mat_ld_sas) <- labels

#Heatmap

library(corrplot)

corrplot(mat_ld_sas,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Sur Asiática",
         mar = c(0,0,3,0))

snpgdsClose(gen_sas)
```
#POBLACION ESTE DE ASIA- EAS
```r
gen_eas <- snpgdsOpen("EAS.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_eas, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_eas, "snp.id"))

#Ver si estan las posiciones que me interesan

mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel_eas <- ids[pos %in% mis_posiciones]
ids_sel_eas
length(ids_sel_eas)
#Calculo de LD
ld_eas <- snpgdsLDMat(gen_eas,
                      snp.id = ids_sel_eas,
                      method = "r",
                      slide = 0)

mat_ld_eas <- ld_eas$LD^2

#Añadir etiquetas
labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_eas) <- rownames(mat_ld_eas) <- labels

#Heatmap

library(corrplot)

corrplot(mat_ld_eas,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Este Asiática",
         mar = c(0,0,3,0))

snpgdsClose(gen_eas)
```
#Vamos a guardar en un PDF todas las figuras que hemos cread. Para eso, abrimos unn PDF que vamos a darle un título y luego volvemos a cargar todo el código anterior y finalmente cerramos con devoff():
```r
pdf("LD_UGT1A1_VariantesRelevantesClinica_xPoblaciones.pdf",
    width = 8,
    height = 6)

dev.off()
```

# 4.2 Hacemos el Test de Mantel con los 22 rs que tienen mayor variabilidad en todas las poblaciones (MAF>0.25)
```r
gen_eur <- snpgdsOpen("EUR.gds")
#En lugar de filtrar por posición manual, dejamos que R filtre por calidad (MAF): 
ids_filtrados_eur <- snpgdsSelectSNP(gen_eur, 
                                     maf = 0.25,       # Tu filtro de MAF > 0.25
                                     missing.rate = 0.1) 

#Vemos cuántas variantes ha rescatado en esta población y calculamos la matriz de DL
message("Variantes rescatadas en EUR: ", length(ids_filtrados_eur))
ld_grande_eur <- snpgdsLDMat(gen_eur,
                             snp.id = ids_filtrados_eur,
                             method = "r",
                             slide = 0)

mat_ld_22_eur <- ld_grande_eur$LD^2
mat_ld_22_eur

#Guardamos las posiciones exactas de ESTOS SNPs para poder sincronizar las matrices después
pos_22_eur <- read.gdsn(index.gdsn(gen_eur, "snp.position"))[ids_filtrados_eur]

snpgdsClose(gen_eur)
```
#POBLACIÓN AFRICANA
```r
gen_afr <- snpgdsOpen("AFR.gds")
ids_filt_afr <- snpgdsSelectSNP(gen_afr, maf = 0.25, missing.rate = 0.1)
pos_22_afr   <- read.gdsn(index.gdsn(gen_afr, "snp.position"))[ids_filt_afr]
snpgdsClose(gen_afr)
```
#POBLACIÓN AMERICANA
```r
gen_amr <- snpgdsOpen("AMR.gds")
ids_filt_amr <- snpgdsSelectSNP(gen_amr, maf = 0.25, missing.rate = 0.1)
pos_22_amr   <- read.gdsn(index.gdsn(gen_amr, "snp.position"))[ids_filt_amr]
snpgdsClose(gen_amr)
```
#POBLACIÓN ESTE DE ASIA (EAS)
```r
gen_eas <- snpgdsOpen("EAS.gds")
ids_filt_eas <- snpgdsSelectSNP(gen_eas, maf = 0.25, missing.rate = 0.1)
pos_22_eas   <- read.gdsn(index.gdsn(gen_eas, "snp.position"))[ids_filt_eas]
snpgdsClose(gen_eas)
```
#POBLACIÓN SUR ASIÁTICA (SAS)
```r
gen_sas <- snpgdsOpen("SAS.gds")
ids_filt_sas <- snpgdsSelectSNP(gen_sas, maf = 0.25, missing.rate = 0.1)
pos_22_sas   <- read.gdsn(index.gdsn(gen_sas, "snp.position"))[ids_filt_sas]
snpgdsClose(gen_sas)
```
#Buscamos los SNPs que pasaron el filtro en TODAS las poblaciones
```r
posiciones_comunes <- intersect(pos_22_eur, pos_22_afr) %>%
  intersect(pos_22_amr) %>%
  intersect(pos_22_eas) %>%
  intersect(pos_22_sas)

message("¡Perfecto! SNPs hipervariables compartidos para Mantel: ", length(posiciones_comunes))
#Vemos que comparten 14 snps con maf> 0.25
#Volvemos a generar las matrices de DL. Como vamos  usar las mismas posiciones, no se necesita renombrar filas ni columnas a mano:

#EUR
gen_eur <- snpgdsOpen("EUR.gds")
ids_comunes_eur <- read.gdsn(index.gdsn(gen_eur, "snp.id"))[read.gdsn(index.gdsn(gen_eur, "snp.position")) %in% posiciones_comunes]
mat_ld_eur <- (snpgdsLDMat(gen_eur, snp.id = ids_comunes_eur, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_eur)

#AFR

gen_afr <- snpgdsOpen("AFR.gds")
ids_comunes_afr <- read.gdsn(index.gdsn(gen_afr, "snp.id"))[read.gdsn(index.gdsn(gen_afr, "snp.position")) %in% posiciones_comunes]
mat_ld_afr <- (snpgdsLDMat(gen_afr, snp.id = ids_comunes_afr, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_afr)

#AMR 

gen_amr <- snpgdsOpen("AMR.gds")
ids_comunes_amr <- read.gdsn(index.gdsn(gen_amr, "snp.id"))[read.gdsn(index.gdsn(gen_amr, "snp.position")) %in% posiciones_comunes]
mat_ld_amr <- (snpgdsLDMat(gen_amr, snp.id = ids_comunes_amr, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_amr)

#EAS

gen_eas <- snpgdsOpen("EAS.gds")
ids_comunes_eas <- read.gdsn(index.gdsn(gen_eas, "snp.id"))[read.gdsn(index.gdsn(gen_eas, "snp.position")) %in% posiciones_comunes]
mat_ld_eas <- (snpgdsLDMat(gen_eas, snp.id = ids_comunes_eas, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_eas)

#Recálculo SAS

gen_sas <- snpgdsOpen("SAS.gds")
ids_comunes_sas <- read.gdsn(index.gdsn(gen_sas, "snp.id"))[read.gdsn(index.gdsn(gen_sas, "snp.position")) %in% posiciones_comunes]
mat_ld_sas <- (snpgdsLDMat(gen_sas, snp.id = ids_comunes_sas, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_sas)
```

#########Para saber cuales son los 14 rs comunes donde se aplica el test de mantel
```r
#Abrimos cualquiera de los archivos GDS
gen_eur <- snpgdsOpen("EUR.gds")

#Extraemos todos los rsIDs y todas las posiciones globales del archivo
todos_los_rsids     <- read.gdsn(index.gdsn(gen_eur, "snp.rs.id"))
todas_las_posiciones <- read.gdsn(index.gdsn(gen_eur, "snp.position"))

#Filtramos los rsIDs que coincidan exactamente con tus 14 posiciones comunes
rsids_14_comunes <- todos_los_rsids[todas_las_posiciones %in% posiciones_comunes]

#Cerramos el archivo para evitar bloqueos de memoria
snpgdsClose(gen_eur)

#Ver exactamente cuales son las 14 coordenadas genómicas en el genoma de referencia GRCh37
print(posiciones_comunes)

#TEST DE MANTEL 

library(vegan)

message("--- RESULTADOS DEL TEST DE MANTEL (SET EXPANDIDO) ---")
mantel(mat_ld_eur, mat_ld_afr, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_amr, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_eas, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_sas, method = "pearson", permutations = 9999)
```
#HACER EL MISMO ANALISIS DE LD Y TEST DE MANTEL CON VARIANTES MAF>0.1#
```r
gen_eur <- snpgdsOpen("EUR.gds")

ids_filtrados_eur <- snpgdsSelectSNP(gen_eur, 
                                     maf = 0.1,       # Tu filtro de MAF > 0.25
                                     missing.rate = 0.1) 

#Ver cuántas variantes ha rescatado en esta población
message("Variantes rescatadas en EUR: ", length(ids_filtrados_eur))

#Calculamos la matriz de DL con el set grande
ld_grande_eur <- snpgdsLDMat(gen_eur,
                             snp.id = ids_filtrados_eur,
                             method = "r",
                             slide = 0)

mat_ld_22_eur <- ld_grande_eur$LD^2
mat_ld_22_eur

#Guardamos las posiciones exactas de ESTOS SNPs para poder sincronizar las matrices después
pos_22_eur <- read.gdsn(index.gdsn(gen_eur, "snp.position"))[ids_filtrados_eur]

snpgdsClose(gen_eur)
```
#POBLACIÓN AFRICANA
```r
gen_afr <- snpgdsOpen("AFR.gds")
ids_filt_afr <- snpgdsSelectSNP(gen_afr, maf = 0.1, missing.rate = 0.1)
pos_22_afr   <- read.gdsn(index.gdsn(gen_afr, "snp.position"))[ids_filt_afr]
snpgdsClose(gen_afr)
```
#POBLACIÓN AMERICANA
```r
gen_amr <- snpgdsOpen("AMR.gds")
ids_filt_amr <- snpgdsSelectSNP(gen_amr, maf = 0.1, missing.rate = 0.1)
pos_22_amr   <- read.gdsn(index.gdsn(gen_amr, "snp.position"))[ids_filt_amr]
snpgdsClose(gen_amr)
```
#POBLACIÓN ESTE DE ASIA (EAS)
```r
gen_eas <- snpgdsOpen("EAS.gds")
ids_filt_eas <- snpgdsSelectSNP(gen_eas, maf = 0.1, missing.rate = 0.1)
pos_22_eas   <- read.gdsn(index.gdsn(gen_eas, "snp.position"))[ids_filt_eas]
snpgdsClose(gen_eas)
```
#POBLACIÓN SUR ASIÁTICA (SAS)
```r
gen_sas <- snpgdsOpen("SAS.gds")
ids_filt_sas <- snpgdsSelectSNP(gen_sas, maf = 0.1, missing.rate = 0.1)
pos_22_sas   <- read.gdsn(index.gdsn(gen_sas, "snp.position"))[ids_filt_sas]
snpgdsClose(gen_sas)
```
#SNPs
```r
posiciones_comunes <- intersect(pos_22_eur, pos_22_afr) %>%
  intersect(pos_22_amr) %>%
  intersect(pos_22_eas) %>%
  intersect(pos_22_sas)

message("¡Perfecto! SNPs hipervariables compartidos para Mantel: ", length(posiciones_comunes))

#Recálculo EUR

gen_eur <- snpgdsOpen("EUR.gds")
ids_comunes_eur <- read.gdsn(index.gdsn(gen_eur, "snp.id"))[read.gdsn(index.gdsn(gen_eur, "snp.position")) %in% posiciones_comunes]
mat_ld_eur <- (snpgdsLDMat(gen_eur, snp.id = ids_comunes_eur, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_eur)

#AFR 

gen_afr <- snpgdsOpen("AFR.gds")
ids_comunes_afr <- read.gdsn(index.gdsn(gen_afr, "snp.id"))[read.gdsn(index.gdsn(gen_afr, "snp.position")) %in% posiciones_comunes]
mat_ld_afr <- (snpgdsLDMat(gen_afr, snp.id = ids_comunes_afr, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_afr)

#AMR

gen_amr <- snpgdsOpen("AMR.gds")
ids_comunes_amr <- read.gdsn(index.gdsn(gen_amr, "snp.id"))[read.gdsn(index.gdsn(gen_amr, "snp.position")) %in% posiciones_comunes]
mat_ld_amr <- (snpgdsLDMat(gen_amr, snp.id = ids_comunes_amr, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_amr)

#EAS

gen_eas <- snpgdsOpen("EAS.gds")
ids_comunes_eas <- read.gdsn(index.gdsn(gen_eas, "snp.id"))[read.gdsn(index.gdsn(gen_eas, "snp.position")) %in% posiciones_comunes]
mat_ld_eas <- (snpgdsLDMat(gen_eas, snp.id = ids_comunes_eas, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_eas)

#SAS

gen_sas <- snpgdsOpen("SAS.gds")
ids_comunes_sas <- read.gdsn(index.gdsn(gen_sas, "snp.id"))[read.gdsn(index.gdsn(gen_sas, "snp.position")) %in% posiciones_comunes]
mat_ld_sas <- (snpgdsLDMat(gen_sas, snp.id = ids_comunes_sas, method = "r", slide = 0)$LD)^2
snpgdsClose(gen_sas)

```
#Para saber cuales son los rs (MAF>0.1) comunes donde se aplica el test de mantel
```r
gen_eur <- snpgdsOpen("EUR.gds")

todos_los_rsids     <- read.gdsn(index.gdsn(gen_eur, "snp.rs.id"))
todas_las_posiciones <- read.gdsn(index.gdsn(gen_eur, "snp.position"))

rsids_comunes <- todos_los_rsids[todas_las_posiciones %in% posiciones_comunes]
snpgdsClose(gen_eur)

print(posiciones_comunes)

#EJECUCIÓN DEL TEST DE MANTEL ROBUSTO(MAF>0.1)

library(vegan)

message("--- RESULTADOS DEL TEST DE MANTEL (SET EXPANDIDO) ---")
mantel(mat_ld_eur, mat_ld_afr, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_amr, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_eas, method = "pearson", permutations = 9999)
mantel(mat_ld_eur, mat_ld_sas, method = "pearson", permutations = 9999)

```