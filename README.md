# TFE_Carla-Gonzalez-de-la-Cruz: Relevancia clínica poblacional de variantes del gen UGT1A1 asociadas a toxicidad por Irinotecan

Este repositorio contiene el flujo de trabajo bioinformático desarrollado para el Trabajo Fin de Máster:

**“Relevancia clínica poblacional de variantes del gen UGT1A1 asociadas a toxicidad por irinotecán”**

El objetivo principal es caracterizar la variabilidad genética del gen *UGT1A1* en distintas poblaciones humanas, estudiar los patrones de desequilibrio de ligamiento (LD) entre variantes farmacogenéticas relevantes y evaluar su posible impacto en la respuesta clínica al tratamiento con irinotecán.Se utilizan los datos genómicos procedentes de bases de datos públicas: 1000genomes y genomAD v2.1.1. 

## Fundamentación biológica

El gen *UGT1A1* codifica una enzima clave en la glucuronidación del metabolito activo del irinotecán (SN-38). Variantes genéticas en este gen pueden modificar la actividad enzimática y aumentar el riesgo de desarrollar efectos adversos graves, especialmente neutropenia severa, toxicidad gastrointestinal y alteraciones en el metabolismo del fármaco. 
La frecuencia de estas variantes y su patrón de herencia conjunta (desequilibrio de ligamiento) difieren entre poblaciones, lo que puede tener implicaciones relevantes en farmacogenética y medicina personalizada.

## Objetivos del estudio

Este análisis pretende:

1. Obtener todas las variantes localizadas en la región genómica de UGT1A1.
2. Comparar las frecuencias alélicas entre poblaciones mundiales.
3. Anotar las variantes mediante identificadores rsID.
4. Identificar variantes farmacogenéticas de interés clínico.
5. Analizar la estructura de desequilibrio de ligamiento dentro del gen.
6. Comparar los patrones de LD entre poblaciones.
7. Estudiar la diversidad genética poblacional utilizando gnomAD.
8. Realizar la anotación funcional de variantes relevantes.

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

Posteriormente se integran todas las frecuencias en una única tabla comparativa. Este análisis permite identificar variantes. 

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
* Catálogos de variantes anotadas con rsID.
* Matrices de desequilibrio de ligamiento.
* Mapas de calor comparativos entre poblaciones.
* Informes de diversidad genética.
* Anotaciones funcionales de variantes relevantes.
* Archivos PDF y Excel para interpretación de resultados.
---
## Relevancia clínica
Este estudio proporciona una caracterización detallada de la variabilidad genética y de la estructura de desequilibrio de ligamiento del gen UGT1A1 en diferentes poblaciones humanas.
Los resultados contribuyen a mejorar la interpretación de biomarcadores farmacogenéticos asociados a la toxicidad por irinotecán y aportan información útil para futuras estrategias de medicina personalizada y farmacogenómica de precisión.


## 1. OBTENCIÓN DE LOS ARCHIVOS VCF DE SECUENCIACION DE LAS BASES DE DATOS PÚBLICAS DEL GEN UGT1A (CROMOSOMA DE REFERENCIA CHr37)

# 1.1 Base de datos 1000Genomes
#Se descarga el cromosoma 2 y el índice
```r
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
```
#Filtrado del gen UGT1A1
```r
bcftools view -r 2:234668000-234700000 \ ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \ -Oz -o UGT1A1_1000G.vcf.gz
```
#Se descarga el archivo de las poblaciones incluidas en la base de datos y se comprueba que se ha obtenido correctamente
```r
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
head integrated_call_samples_v3.20130502.ALL.panel.1
```
#Creacion de la lista de poblaciones
```r
awk '$3=="AFR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > AFR.txt
awk '$3=="EUR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > EUR.txt
awk '$3=="EAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > EAS.txt
awk '$3=="SAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > SAS.txt
awk '$3=="AMR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > AMR.txt

wc -l AFR.txt EUR.txt EAS.txt SAS.txt AMR.txt #comprobar que las listas no están vacías
```
#Calculo de las frecuencias de las diferentes variantes en cada poblacion sin rsID
```r
bcftools view -S AFR.txt UGT1A1_1000G.vcf.gz | \ bcftools +fill-tags -- -t AF | \ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > AFR_freq.txt
bcftools view -S EUR.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > EUR_freq.txt
bcftools view -S EAS.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > EAS_freq.txt
bcftools view -S SAS.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > SAS_freq.txt
bcftools view -S AMR.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > AMR_freq.txt
```
#Creacion de una clave unica para cada archivo de poblacion para que nos de las posiciones correctas y las frecuencias como en el ejemplo: 2   234668879   TA   T   0.12 PASE A: 2:234668879:TA:T    0.12. 
```r
awk '{print $1":"$2":"$3":"$4"\t"$5}' AFR_freq.txt > AFR_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' EUR_freq.txt > EUR_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' EAS_freq.txt > EAS_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' SAS_freq.txt > SAS_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' AMR_freq.txt > AMR_key.txt
```
#Ordenamos las columnas de cada arhicvo y las juntamos en uno único con el comando JOIN
```r
sort -k1,1 AFR_key.txt -o AFR_key.txt
sort -k1,1 EUR_key.txt -o EUR_key.txt
sort -k1,1 EAS_key.txt -o EAS_key.txt
sort -k1,1 SAS_key.txt -o SAS_key.txt
sort -k1,1 AMR_key.txt -o AMR_key.txt

join AFR_key.txt EUR_key.txt > tmp1.txt
join tmp1.txt EAS_key.txt > tmp2.txt
join tmp2.txt SAS_key.txt > tmp3.txt
join tmp3.txt AMR_key.txt > tabla_final.txt

sed '1iVARIANT\tAFR\tEUR\tEAS\tSAS\tAMR' tabla_final.txt > tabla_final_con_header.txt #se crea un unico archivio y se añade la cabezera
```
#El archivo tabla_final_con_header.txt obtiene todas las variantes descritas con posiciones del gen UGTA1 y la frecuencia en cada poblacion. 
#Exportamos  a excel ambos archivos y descargamos en nuestro usuario
```r
tr '\t' ',' < tabla_final_con_header.txt > tabla_final_con_header.txt.csv
cp tabla_final_con_header.txt /mnt/c/Users/Usuario/Downloads/
```
# 1.2 Base de datos genomAD
Seguimos los mismos pasos, pero consultando los datos de la base de datos genomAD. 
```r
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.vcf.bgz.tbi
```
#Filtrado del gen UGT1A1
```r
bcftools view -r 2:234668000-234700000 \ gnomad.genomes.r2.1.1.sites.2.vcf.bgz \ -Oz -o UGT1A1_gnomad.vcf.gz
```
#Extraccion del archivo de frecuencias
```r
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%INFO/AF_afr\t%INFO/AF_nfe\t%INFO/AF_eas\t%INFO/AF_amr\n' \
UGT1A1_gnomad.vcf.gz > gnomad_raw.txt
```
#Creacion del identificador de las variantes y añadir la cabecera
```r
$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' gnomad_raw.txt > gnomad_clean.txt
sed '1iVARIANT\tAF\tAFR\tEUR\tEAS\tSAS\tAMR' gnomad_clean.txt > gnomad_final.txt
```
#Exportamos a excel y guardamos resultado
```r
tr '\t' ',' < gnomad_final.txt > gnomad_final.csv
cp gnomad_final.csv /mnt/c/Users/Usuario/Downloads/
```
#Filtramos por variantes relevantes con frecuencia > 0.1 en cualquier poblacion
```r
awk '$3 > 0.1 || $4 > 0.1 || $5 > 0.1 || $6 > 0.1' gnomad_final.txt > variantes_relevantes_0.1.txt 
```
#Exportamos este archivo a excel y guardamos
```r
tr '\t' ',' < variantes_relevantes_0.01.txt > variantes_relevantescompletas_0.01.csv
cp variantes_relevantescompletas_0.01.csv /mnt/c/Users/Usuario/Downloads/
```
## 2. CRUZAR LAS VARIANTES OBTENIDAS, NUESTRAS COORDENADAS CON LOS RSID QUE SABEMOS. COMO NOS DA ERROR CRUZAR CON ENSEMBL, PROBAMOS A CRUZAR AMBAS BASES DE DATOS CON sdSNP (CROMOSOMA DE REFERENCIA CHr37)
#Descargamos dbSNP (CHr37)
```r
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
```
#Indexar VCF
```r
bcftools index UGT1A1_gnomad.vcf.gz
bcftools index UGT1A1_1000G.vcf.gz
```
#Normalizar las variantes
#Descargamos primero Genoma Chr37(hg19)
```r
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
```
#Lo descomprimimos quedandonos la secuencia fasta 
```r
gunzip human_g1k_v37.fasta.gz
human_g1k_v37.fasta
```
#Indexar
```r
samtools faidx human_g1k_v37.fasta
```
#Este comando crea el archivo fasta descomprimido e indexado human_g1k_v37.fasta.fai

#Normalización para alinear las variantes correctamente 
```r
bcftools norm -f human_g1k_v37.fasta \
UGT1A1_gnomad.vcf.gz -Oz -o gnomad.norm.vcf.gz

bcftools norm -f human_g1k_v37.fasta \
UGT1A1_1000G.vcf.gz -Oz -o 1000G.norm.vcf.gz
```
#Indexamos de nuevo ambos archivos 
```r
bcftools index -t gnomad.norm.vcf.gz
bcftools index -t 1000G.norm.vcf.gz
```
#AÑADIR rs ID a genomAD

#Es necesario ver como se llaman los cromosomas en los dos archivos antes: 
```r
bcftools view -H gnomad.norm.vcf.gz | head -n 5 | cut -f1
bcftools view -H GCF_000001405.25.gz | head -n 5 | cut -f1
bcftools view -H 1000G.norm.vcf.gz | head -n 5 | cut -f1
```
#Como nos da dos archivos que se llaman 2 y el archivo de dbSNP nos devuelve la posicion NC_000001.10, hay que reenombrarlos para que se crucen correctamente:
```r
cat <<EOF > mapping_chrs.txt
NC_000001.10 1
NC_000002.11 2
NC_000003.11 3
NC_000004.11 4
NC_000005.9  5
NC_000006.11 6
NC_000007.13 7
NC_000008.10 8
NC_000009.11 9
NC_000010.10 10
NC_000011.9  11
NC_000012.11 12
NC_000013.10 13
NC_000014.8  14
NC_000015.9  15
NC_000016.9  16
NC_000017.10 17
NC_000018.9  18
NC_000019.9  19
NC_000020.10 20
NC_000021.8  21
NC_000022.10 22
NC_000023.10 X
NC_000024.9  Y
NC_012920.1  MT
EOF
```
#Renombramos para que la posicion del cromosoma 2 la nombre como 1. Ahora generamos una nueva versión del archivo de referencia (el que tiene los rsID) pero con los nombres de cromosomas simplificados:
```r
bcftools annotate --rename-chrs mapping_chrs.txt GCF_000001405.25.gz -Oz -o ref_renamed.vcf.gz

#Indexamos este archivo renombrado

tabix -p vcf ref_renamed.vcf.gz

bcftools annotate \
  -a ref_renamed.vcf.gz \
  -c ID \
  gnomad.norm.vcf.gz \
  -Oz -o gnomad.with_rsID.vcf.gz
```
#Comprobamos que se ha ejecutaado correctamente:
```r
bcftools view -H gnomad.with_rsID.vcf.gz | head -n 20 | awk '{print $1, $2, $3}'
```
#Añadir rs ID a 1000Genomes
```r
bcftools annotate \
-a ref_renamed.vcf.gz \
-c ID \
1000G.norm.vcf.gz \
-Oz -o 1000G.with_rsID.vcf.gz

#Comprobamos que se ha ejecutaado correctamente

bcftools view -H 1000G.with_rsID.vcf.gz | head -n 20 | awk '{print $1, $2, $3}'

#Indexar nuevamente ambos archivos 

tabix -p vcf gnomad.with_rsID.vcf.gz
tabix -p vcf 1000G.with_rsID.vcf.gz
```

#Extraer los datos y convertirlos a formato CSV (cambiando tabuladores por comas)

#genomAD:
```r
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_popmax\t%INFO/AF_nfe\t%INFO/AF_nfe_seu\t%INFO/AF_amr\t%INFO/AF_afr\t%INFO/AF_eas\t%INFO/AF_asj\t%INFO/AF_fin\t%INFO/AF_oth\t%INFO/AF_male\t%INFO/AF_female\n' gnomad.with_rsID.vcf.gz > datos_poblaciones_total.tsv

#Conversion de archivo tsv a csv y creacion de la cabecera (títulos de las columnas)
echo "CHR,POS,rsID,REF,ALT,AF_Global,AF_PopMax,AF_Europeos,AF_Europa_Sur,AF_Latino,AF_Africano,AF_Asiatico,AF_Judio,AF_Fin,AF_Otros,AF_Hombres,AF_Mujeres" > UGT1A1_final.csv

#Convertir los tabuladores del archivo .tsv en comas y añadirlo al archivo:

tr '\t' ',' < datos_poblaciones_total.tsv >> UGT1A1_gnomADfinal.csv
```

#1000Genomes:
```r
bcftools annotate \
  -a ref_renamed.vcf.gz \
  -c ID \
  1000G.norm.vcf.gz \
  -Oz -o 1000G.with_rsID.vcf.gz

#Indexar el nuevo archivo:

tabix -p vcf 1000G.with_rsID.vcf.gz

#Crear la cabecera

echo "CHR,POS,rsID,REF,ALT,AF_Global,AF_Europeos,AF_Africanos,AF_Americanos,AF_Asiaticos_Este,AF_Asiaticos_Sur" > Final_1000G_con_RS.csv

#Extraer los datos del archivo NUEVO (1000G.with_rsID.vcf.gz):

bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/EUR_AF\t%INFO/AFR_AF\t%INFO/AMR_AF\t%INFO/EAS_AF\t%INFO/SAS_AF\n' 1000G.with_rsID.vcf.gz | tr '\t' ',' >> Final_1000G_con_RS.csv
```

#Exportar al ordenador:
```r
explorer.exe 
```
# 3. ANALISIS DE DL EN LA BASE DE DATOS DE 1000GENOMES PORQUE ES LA QUE ESTÁ INDIVIDUO POR INDIVIDUO Y EN FASE
# 3.1 Quedarnos solamente con variantes relevantes 
#Primero se hace un filtrado global y luego poblacion por poblacion
```r
 bcftools view -v snps UGT1A1_1000G.vcf.gz | \
bcftools filter -i 'MAF>0.05 && F_MISSING<0.1' \
-Oz -o UGT1A1_1000G_filtered.vcf.gz

#Indexar
bcftools index

#Filtro por poblacion en EUR:

bcftools view -S EUR.txt UGT1A1_1000G_filtered.vcf.gz -Oz -o UGT1A1_1000G_filtered_EUR.vcf.gz
bcftools index UGT1A1_1000G_filtered_EUR.vcf.gz

#Filtro por poblacion  EAS

bcftools view -S EAS.txt UGT1A1_1000G_filtered.vcf.gz -Oz -o UGT1A1_1000G_filtered_EAS.vcf.gz
bcftools index UGT1A1_1000G_filtered_EAS.vcf.gz

#Filtro por poblacion AFR_AF

bcftools view -S AFR.txt UGT1A1_1000G_filtered.vcf.gz -Oz -o UGT1A1_1000G_filtered_AFR.vcf.gz
bcftools index UGT1A1_1000G_filtered_AFR.vcf.gz

#Filtro por poblacion AMR:

bcftools view -S AMR.txt UGT1A1_1000G_filtered.vcf.gz -Oz -o UGT1A1_1000G_filtered_AMR.vcf.gz
bcftools index UGT1A1_1000G_filtered_AMR.vcf.gz

#Filtro por poblacion SAS:

bcftools view -S SAS.txt UGT1A1_1000G_filtered.vcf.gz -Oz -o UGT1A1_1000G_filtered_SAS.vcf.gz
bcftools index UGT1A1_1000G_filtered_SAS.vcf.gz
```

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
# 4.2 Hacer test estadistico para comparar el LD de las diferentes poblaciones
```r
library(vegan)

#Coger las matrices calculadas para cada poblacion y poner el mismo orden.
#Primero aplicamos el nombre de las columnas que tenemos para la matriz de europa(mat_ld) al siguiente:

colnames(mat_ld_afr) <- labels
rownames(mat_ld_afr) <- labels

colnames(mat_ld_amr) <- labels
rownames(mat_ld_amr) <- labels

colnames(mat_ld_eas) <- labels
rownames(mat_ld_eas) <- labels

colnames(mat_ld_eas) <- labels
rownames(mat_ld_eas) <- labels

#Aplicar Test de Mantel cogiendo como indicativo la poblacion europea

mantel(mat_ld, mat_ld_amr, method = "pearson", permutations = 999999)
mantel(mat_ld, mat_ld_afr, method = "pearson", permutations = 999999)
mantel(mat_ld, mat_ld_eas, method = "pearson", permutations = 999999)
mantel(mat_ld, mat_ld_sas, method = "pearson", permutations = 999999)
```
# 4.3 Hacemos el Test de Mantel con los 22 rs que tienen mayor variabilidad en todas las poblaciones (MAF>0.25)
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
## 5. ANALISIS DESCRIPTIVO DE LA BASE DE DATOS DE GENOMAD
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
## 6. ANALISIS DE QUE TIPO DE VARIANTES SON LAS OBTENIDAS EN LAS 22 RS CON MAF > 0.25 VÍA API DE ENSEMBL
```r

library(httr2)
library(jsonlite)
library(openxlsx)

#22 rsIDs originales
rsids <- c("rs759174", "rs28946889", "rs4148326", "rs4663971", "rs4148328", 
           "rs4148329", "rs6717546", "rs6719561", "rs2003363", "rs56133230", 
           "rs10203853", "rs10209214", "rs6728940", "rs6746002", "rs1500475", 
           "rs6728520", "rs6735414", "rs10202032", "rs10199512", "rs4663336", 
           "rs6723936", "rs62648950")

url_vep <- "https://grch37.rest.ensembl.org/vep/human/id?hgvs=1"
body_data <- list(ids = rsids)

#Consulta a la API
message("1. Consultando Ensembl VEP (Extrayendo Dataset Máximo)...")
req <- request(url_vep) %>%
  req_headers(`Content-Type` = "application/json", `Accept` = "application/json") %>%
  req_body_raw(jsonlite::toJSON(body_data, auto_unbox = TRUE), "application/json") %>%
  req_timeout(20) # Evita bloqueos en el servidor

res <- req_perform(req)

#Procesamiento de los campos
if (resp_status(res) == 200) {
  raw_data <- resp_body_string(res) %>% fromJSON(flatten = TRUE)
  
  message("2. Desempaquetando variables globales y consecuencias de transcritos...")
  lista_total <- list()
  contador <- 1
  
  for (i in 1:nrow(raw_data)) {
    variante <- raw_data[i, ]
    
    # --- Datos Nivel Variante ---
    rsid       <- variante$id
    cromosoma  <- if("seq_region_name" %in% names(variante)) variante$seq_region_name else NA
    pos_inicio <- if("start" %in% names(variante)) variante$start else NA
    pos_fin    <- if("end" %in% names(variante)) variante$end else NA
    hebra      <- if("strand" %in% names(variante)) variante$strand else NA
    localizacion_global <- paste0("chr", cromosoma, ":", pos_inicio, "-", pos_fin, " (", ifelse(hebra == 1, "+", "-"), ")")
    
    # Bloque de transcritos
    consecuencias <- variante$transcript_consequences[[1]]
    
    if (!is.null(consecuencias) && nrow(consecuencias) > 0) {
      for (j in 1:nrow(consecuencias)) {
        tx <- consecuencias[j, ]
        
        # Extracción y mapeo de los campos que has solicitado 
        tx_id        <- if("transcript_id" %in% names(tx)) tx$transcript_id else NA
        coding_impact<- if("consequence_terms" %in% names(tx)) paste(unlist(tx$consequence_terms), collapse = ", ") else NA
        gen          <- if("gene_symbol" %in% names(tx)) tx$gene_symbol else NA
        
        # Nomenclaturas HGVS puras de Ensembl
        hgvs_coding  <- if("hgvsc" %in% names(tx)) tx$hgvsc else "No disponible/No codificante"
        hgvs_protein <- if("hgvsp" %in% names(tx)) tx$hgvsp else "No modifica proteína"
        
        # Datos moleculares del transcrito
        splice_dist  <- if("splice_distance" %in% names(tx)) tx$splice_distance else "N/A"
        tsl          <- if("tsl" %in% names(tx)) tx$tsl else NA
        mane_select  <- if("mane_select" %in% names(tx)) tx$mane_select else "N/A"
        appris       <- if("appris" %in% names(tx)) tx$appris else "N/A"
        
        # Accesiones a bases de datos de proteínas externas
        uniprot_acc  <- if("swissprot" %in% names(tx)) {
          paste(unlist(tx$swissprot), collapse = ", ")
        } else if ("uniprot_isoform" %in% names(tx)) {
          paste(unlist(tx$uniprot_isoform), collapse = ", ")
        } else {
          "N/A"
        }
        
        biotype      <- if("biotype" %in% names(tx)) tx$biotype else NA
        
        # Guardar todo estructurado en la fila
        lista_total[[contador]] <- data.frame(
          RSID = rsid,
          Location_hg19 = localizacion_global,
          Chrom = cromosoma,
          Pos_Start = pos_inicio,
          Transcript = tx_id,
          Gene = gen,
          Coding_Impact = coding_impact,
          HGVS_Coding = hgvs_coding,
          HGVS_Protein = hgvs_protein,
          Splice_Distance = splice_dist,
          TSL = tsl,
          MANE_Select = mane_select,
          APPRIS = appris,
          UniProt_Accession = uniprot_acc,
          Biotype = biotype,
          stringsAsFactors = FALSE
        )
        contador <- contador + 1
      }
    }
  }
  
  tabla_maxima <- do.call(rbind, lista_total)
  
  #Exporotar a Excel
  message("3. Escribiendo reporte masivo en Excel...")
  wb <- createWorkbook()
  addWorksheet(wb, "Anotación Completa VEP")
  
  #Estilo
  estilo_top <- createStyle(fontName = "Calibri", fontColour = "#FFFFFF", 
                            textDecoration = "bold", bgFill = "#1F4E79", halign = "center")
  
  writeData(wb, "Anotación Completa VEP", tabla_maxima, headerStyle = estilo_top)
  setColWidths(wb, "Anotación Completa VEP", cols = 1:ncol(tabla_maxima), widths = "auto")
  
  saveWorkbook(wb, "Anotacion_Completa_Ensembl_UGT1A.xlsx", overwrite = TRUE)
  
  message("¡PROCESO COMPLETADO!")
  message("Revisa tu directorio de trabajo. Se ha generado: 'Anotacion_Completa_Ensembl_UGT1A.xlsx'")
  
} else {
  stop("El servidor de Ensembl no ha respondido de forma correcta. Código: ", resp_status(res))
}
```


