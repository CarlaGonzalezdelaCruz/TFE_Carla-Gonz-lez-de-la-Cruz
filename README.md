# TFE_Carla-Gonzalez-de-la-Cruz
Relevancia clínica poblacional de variantes del gen UGT1A1 asociadas a toxicidad por Irinotecan

En el siguiente repositorio se adjunta los códigos empleados en la metodología del Trabajo Fin de Estudios. 

## 1. OBTENCIÓN DE LOS ARCHIVOS VCF DE SECUENCIACION DE LAS BASES DE DATOS PÚBLICAS DEL GEN UGT1A (CROMOSOMA DE REFERENCIA CHr37)
# 1.1 Base de datos 1000Genomes
#Se descarga el cromosoma 2 y el índice
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi

# Filtrado del gen UGT1A1
bcftools view -r 2:234668000-234700000 \ ALL.chr2.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz \ -Oz -o UGT1A1_1000G.vcf.gz
#Se descarga el archivo de las poblaciones incluidas en la base de datos y se comprueba que se ha obtenido correctamente
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
head integrated_call_samples_v3.20130502.ALL.panel.1

# Creacion de la lista de poblaciones
awk '$3=="AFR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > AFR.txt
awk '$3=="EUR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > EUR.txt
awk '$3=="EAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > EAS.txt
awk '$3=="SAS"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > SAS.txt
awk '$3=="AMR"{print $1}' integrated_call_samples_v3.20130502.ALL.panel > AMR.txt

wc -l AFR.txt EUR.txt EAS.txt SAS.txt AMR.txt #comprobar que las listas no están vacías

# Calculo de las frecuencias de las diferentes variantes en cada poblacion sin rsID

bcftools view -S AFR.txt UGT1A1_1000G.vcf.gz | \ bcftools +fill-tags -- -t AF | \ bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > AFR_freq.txt
bcftools view -S EUR.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > EUR_freq.txt
bcftools view -S EAS.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > EAS_freq.txt
bcftools view -S SAS.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > SAS_freq.txt
bcftools view -S AMR.txt UGT1A1_1000G.vcf.gz | bcftools +fill-tags -- -t AF | bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\n' > AMR_freq.txt

# Creacion de una clave unica para cada archivo de poblacion para que nos de las posiciones correctas y las frecuencias como en el ejemplo: 2   234668879   TA   T   0.12 PASE A: 2:234668879:TA:T    0.12. 
awk '{print $1":"$2":"$3":"$4"\t"$5}' AFR_freq.txt > AFR_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' EUR_freq.txt > EUR_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' EAS_freq.txt > EAS_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' SAS_freq.txt > SAS_key.txt
awk '{print $1":"$2":"$3":"$4"\t"$5}' AMR_freq.txt > AMR_key.txt

# Ordenamos las columnas de cada arhicvo y las juntamos en uno único con el comando JOIN
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

# El archivo tabla_final_con_header.txt obtiene todas las variantes descritas con posiciones del gen UGTA1 y la frecuencia en cada poblacion. 
# Filtramos este archivo quedándonos solo con las variantes con una frecuencia mayor a 0.1 en algunas de los grupos poblacionales
awk '$2 > 0.1 || $3 > 0.1 || $4 > 0.1 || $5 > 0.1 || $6 > 0.1' tabla_final_con_header.txt > variantes_relevantes_1000G.txt

# Exportamos  a excel ambos archivos y descargamos en nuestro usuario
tr '\t' ',' < tabla_final_con_header.txt > tabla_final_con_header.txt.csv
cp tabla_final_con_header.txt /mnt/c/Users/Usuario/Downloads/

tr '\t' ',' < variantes_relevantes_1000G.txt > variantes_relevantes_1000G.csv
cp variantes_relevantes_1000G.csv /mnt/c/Users/Usuario/Downloads/

# 1.2 Base de datos genomAD
Seguimos los mismos pasos, pero consultando los datos de la base de datos genomAD. 
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.vcf.bgz
wget https://storage.googleapis.com/gcp-public-data--gnomad/release/2.1.1/vcf/genomes/gnomad.genomes.r2.1.1.sites.2.vcf.bgz.tbi

# Filtrado del gen UGT1A1
bcftools view -r 2:234668000-234700000 \ gnomad.genomes.r2.1.1.sites.2.vcf.bgz \ -Oz -o UGT1A1_gnomad.vcf.gz
#Extraccion del archivo de frecuencias
bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\t%AF\t%INFO/AF_afr\t%INFO/AF_nfe\t%INFO/AF_eas\t%INFO/AF_amr\n' \
UGT1A1_gnomad.vcf.gz > gnomad_raw.txt
#Creacion del identificador de las variantes y añadir la cabecera
$1":"$2":"$3":"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' gnomad_raw.txt > gnomad_clean.txt
sed '1iVARIANT\tAF\tAFR\tEUR\tEAS\tSAS\tAMR' gnomad_clean.txt > gnomad_final.txt

# Exportamos a excel y guardamos resultado
tr '\t' ',' < gnomad_final.txt > gnomad_final.csv
cp gnomad_final.csv /mnt/c/Users/Usuario/Downloads/

#Filtramos por variantes relevantes con frecuencia > 0.1 en cualquier poblacion
awk '$3 > 0.1 || $4 > 0.1 || $5 > 0.1 || $6 > 0.1' gnomad_final.txt > variantes_relevantes_0.1.txt 
#Exportamos este archivo a excel y guardamos
tr '\t' ',' < variantes_relevantes_0.01.txt > variantes_relevantescompletas_0.01.csv
cp variantes_relevantescompletas_0.01.csv /mnt/c/Users/Usuario/Downloads/

## 2. CRUZAR LAS VARIANTES OBTENIDAS, NUESTRAS COORDENADAS CON LOS RSID QUE SABEMOS. COMO NOS DA ERROR CRUZAR CON ENSEMBL, PROBAMOS A CRUZAR AMBAS BASES DE DATOS CON sdSNP (CROMOSOMA DE REFERENCIA CHr37)
# Descargamos dbSNP (CHr37)
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz
wget https://ftp.ncbi.nih.gov/snp/latest_release/VCF/GCF_000001405.25.gz.tbi
# Indexar VCF
bcftools index UGT1A1_gnomad.vcf.gz
bcftools index UGT1A1_1000G.vcf.gz
# Normalizar las variantes
#Descargamos primero Genoma Chr37(hg19)
wget https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz
#Lo descomprimimos quedandonos la secuencia fasta 
gunzip human_g1k_v37.fasta.gz
human_g1k_v37.fasta
# Indexar
samtools faidx human_g1k_v37.fasta
#Este comando crea el archivo fasta descomprimido e indexado human_g1k_v37.fasta.fai
# Normalización para alinear las variantes correctamente 
bcftools norm -f human_g1k_v37.fasta \
UGT1A1_gnomad.vcf.gz -Oz -o gnomad.norm.vcf.gz

bcftools norm -f human_g1k_v37.fasta \
UGT1A1_1000G.vcf.gz -Oz -o 1000G.norm.vcf.gz
# Indexamos de nuevo ambos archivos 
bcftools index -t gnomad.norm.vcf.gz
bcftools index -t 1000G.norm.vcf.gz
# AÑADIR rs ID a genomAD
#Es necesatio ver como se llaman los cromosomas en los dos archivos antes: 
bcftools view -H gnomad.norm.vcf.gz | head -n 5 | cut -f1
bcftools view -H GCF_000001405.25.gz | head -n 5 | cut -f1
bcftools view -H 1000G.norm.vcf.gz | head -n 5 | cut -f1
#Como nos da dos archivos que se llaman 2 y el archivo de dbSNP nos devuelve la posicion NC_000001.10, hay que reenombrarlos para que se crucen correctamente:

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

#Renombramos para que la posicion del cromosoma 2 la nombre como 1. 
Ahora generamos una nueva versión del archivo de referencia (el que tiene los rsID) pero con los nombres de cromosomas simplificados:
bcftools annotate --rename-chrs mapping_chrs.txt GCF_000001405.25.gz -Oz -o ref_renamed.vcf.gz

#indexamos este archivo renombrado
tabix -p vcf ref_renamed.vcf.gz



bcftools annotate \
  -a ref_renamed.vcf.gz \
  -c ID \
  gnomad.norm.vcf.gz \
  -Oz -o gnomad.with_rsID.vcf.gz

#Comprobamos que se ha ejecutaado correctamente
bcftools view -H gnomad.with_rsID.vcf.gz | head -n 20 | awk '{print $1, $2, $3}'

# Añadir rs ID a 1000Genomes
bcftools annotate \
-a ref_renamed.vcf.gz \
-c ID \
1000G.norm.vcf.gz \
-Oz -o 1000G.with_rsID.vcf.gz

#Comprobamos que se ha ejecutaado correctamente
bcftools view -H 1000G.with_rsID.vcf.gz | head -n 20 | awk '{print $1, $2, $3}'

# Indexar nuevamente ambos archivos 
tabix -p vcf gnomad.with_rsID.vcf.gz
tabix -p vcf 1000G.with_rsID.vcf.gz

# Extraer los datos y convertirlos a formato CSV (cambiando tabuladores por comas)
#genomAD
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_popmax\t%INFO/AF_nfe\t%INFO/AF_nfe_seu\t%INFO/AF_amr\t%INFO/AF_afr\t%INFO/AF_eas\t%INFO/AF_asj\t%INFO/AF_fin\t%INFO/AF_oth\t%INFO/AF_male\t%INFO/AF_female\n' gnomad.with_rsID.vcf.gz > datos_poblaciones_total.tsv

#Conversion de archivo tsv a csv
#Crear la cabecera (títulos de las columnas)
echo "CHR,POS,rsID,REF,ALT,AF_Global,AF_PopMax,AF_Europeos,AF_Europa_Sur,AF_Latino,AF_Africano,AF_Asiatico,AF_Judio,AF_Fin,AF_Otros,AF_Hombres,AF_Mujeres" > UGT1A1_final.csv

#Convertir los tabuladores del archivo .tsv en comas y añadirlo al archivo
tr '\t' ',' < datos_poblaciones_total.tsv >> UGT1A1_gnomADfinal.csv

#1000Genomes
bcftools annotate \
  -a ref_renamed.vcf.gz \
  -c ID \
  1000G.norm.vcf.gz \
  -Oz -o 1000G.with_rsID.vcf.gz

#Indexar el nuevo archivo (siempre recomendable)
tabix -p vcf 1000G.with_rsID.vcf.gz
#Crear la cabecera
echo "CHR,POS,rsID,REF,ALT,AF_Global,AF_Europeos,AF_Africanos,AF_Americanos,AF_Asiaticos_Este,AF_Asiaticos_Sur" > Final_1000G_con_RS.csv

#Extraer los datos del archivo NUEVO (1000G.with_rsID.vcf.gz)
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/AF\t%INFO/EUR_AF\t%INFO/AFR_AF\t%INFO/AMR_AF\t%INFO/EAS_AF\t%INFO/SAS_AF\n' 1000G.with_rsID.vcf.gz | tr '\t' ',' >> Final_1000G_con_RS.csv


# Exportar al ordenador
explorer.exe 





# 3. ANALISIS DE DL EN LA BASE DE DATOS DE 1000GENOMES PORQUE ES LA QUE ESTÁ INDIVIDUO POR INDIVIDUO Y EN FASE
# 3.1 Quedarnos solamente con variantes relevantes 
#Primero se hace un filtrado global y luego poblacion por poblacion

 bcftools view -v snps UGT1A1_1000G.vcf.gz | \
bcftools filter -i 'MAF>0.05 && F_MISSING<0.1' \
-Oz -o UGT1A1_1000G_filtered.vcf.gz

#indexar
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


#Una vez tenemos los archivos de las variantes relevantes por población, hacemos el análisis de desequilibrio de ligamiento.
# Instalacion de paquetes necesarios
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("SNPRelate", "gdsfmt"))
# carga de librerias 
library(SNPRelate)
library(gdsfmt)
setwd(utils::choose.dir())
# Hacemos el análisis para variantes conn MAF > 0.25%
# 1. Abrir el archivo
genofile <- snpgdsOpen("EUR.gds")

# 2. SELECCIONAR SNPs con filtro estricto MAF > 0.25
# Esto dejará solo las variantes MUY comunes
snp_seleccionados <- snpgdsSelectSNP(genofile, maf=0.25, missing.rate=0.05)

# 3. EXTRAER LOS NOMBRES REALES (rsID)
# Primero sacamos la lista completa de rsIDs que hay en el archivo
lista_completa_rsids <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))

# Ahora sacamos la lista de TODOS los IDs numéricos
lista_completa_ids <- read.gdsn(index.gdsn(genofile, "snp.id"))

# Filtramos los rsIDs para que coincidan EXACTAMENTE con los SNPs que pasaron el filtro 0.25
nombres_para_grafico <- lista_completa_rsids[lista_completa_ids %in% snp_seleccionados]

# 4. CALCULAR MATRIZ DE LD
ld_obj <- snpgdsLDMat(genofile, snp.id=snp_seleccionados, method="r", slide=0)
matrix_final <- ld_obj$LD^2

# 5. ASIGNAR NOMBRES
colnames(matrix_final) <- nombres_para_grafico
rownames(matrix_final) <- nombres_para_grafico

# 6. CERRAR
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
# A continuacion vemos cuales son las variantes con MAF >0.25 en Europa y comparamos con el resto de poblaciones
genofile <- snpgdsOpen("EUR.gds")

# 2. Filtro MAF > 0.25
snp_sel <- snpgdsSelectSNP(genofile, maf=0.25, missing.rate=0.05)

# 3. EXTRAER DATOS (Forzamos la lectura de rsID y de Posición)
ids_numericos <- read.gdsn(index.gdsn(genofile, "snp.id"))
rsids_vcf <- read.gdsn(index.gdsn(genofile, "snp.rs.id"))
posiciones <- read.gdsn(index.gdsn(genofile, "snp.position"))

# 4. CREAR ETIQUETAS ÚNICAS
# Si el rsID es vacío o ".", usamos la posición genómica (Chr2:pos)
etiquetas <- ifelse(rsids_vcf == "" | rsids_vcf == ".", 
                    paste0("2:", posiciones), 
                    rsids_vcf)

# 5. FILTRAR ETIQUETAS PARA EL GRÁFICO
nombres_finales <- etiquetas[ids_numericos %in% snp_sel]

# 6. CALCULAR LD
ld_obj <- snpgdsLDMat(genofile, snp.id=snp_sel, method="r", slide=0)
matrix_eur <- ld_obj$LD^2
colnames(matrix_eur) <- nombres_finales
rownames(matrix_eur) <- nombres_finales

snpgdsClose(genofile)
# 7. Graficar
corrplot(matrix_eur, 
         method = "color", 
         type = "upper", 
         tl.col = "black", 
         tl.cex = 0.7,      # Tamaño de letra de los rsID
         addgrid.col = "gray", 
         main = "DL en variantes del gen UGT1A1 en población EUROPEA con MAF > 0.25",
         mar = c(0,0,2,0))
# Hacemos lo mismo en poblacion Africana
# Sustituye el nombre entre comillas por el nombre REAL de tu archivo VCF africano
snpgdsVCF2GDS("UGT1A1_1000G_filtered_AFR.vcf.gz", "AFR.gds", method="biallelic.only")
# 1. Abrir AFR
genofile_afr <- snpgdsOpen("AFR.gds")

# 2. USAR LOS MISMOS SNPs QUE EN EUR (Usamos snp_sel que guardamos antes)
ld_obj_afr <- snpgdsLDMat(genofile_afr, snp.id=snp_sel, method="r", slide=0)
matrix_afr <- ld_obj_afr$LD^2

# 3. Ponemos los mismos nombres para que la comparación sea visualmente idéntica
colnames(matrix_afr) <- nombres_finales
rownames(matrix_afr) <- nombres_finales

snpgdsClose(genofile_afr)

#VISUALIZAR AMBOS GRAFICOS 
par(mfrow=c(1,2)) # Divide el gráfico en 2 columnas

library(corrplot)
# Gráfico EUR
corrplot(matrix_eur, method="color", type="upper", tl.cex=0.6, main="EUR (MAF > 0.2)")

# Gráfico AFR
corrplot(matrix_afr, method="color", type="upper", tl.cex=0.6, main="AFR (Mismos SNPs)")

#CONVERTIR LAS OTRAS TRES POBLACIONES
# Convertimos las 3 poblaciones restantes
snpgdsVCF2GDS("UGT1A1_1000G_filtered_EAS.vcf.gz", "EAS.gds", method="biallelic.only")
snpgdsVCF2GDS("UGT1A1_1000G_filtered_AMR.vcf.gz", "AMR.gds", method="biallelic.only")
snpgdsVCF2GDS("UGT1A1_1000G_filtered_SAS.vcf.gz", "SAS.gds", method="biallelic.only")

# codigo comparaito de todas las poblaciones
obtener_matriz <- function(archivo_gds, snps, etiquetas) {
  g <- snpgdsOpen(archivo_gds)
  ld <- snpgdsLDMat(g, snp.id=snps, method="r", slide=0)
  m <- ld$LD^2
  colnames(m) <- etiquetas
  rownames(m) <- etiquetas
  snpgdsClose(g)
  return(m)
}

# Calculamos las matrices restantes
matrix_afr <- obtener_matriz("AFR.gds", snp_sel, nombres_finales)
matrix_eas <- obtener_matriz("EAS.gds", snp_sel, nombres_finales)
matrix_amr <- obtener_matriz("AMR.gds", snp_sel, nombres_finales)
matrix_sas <- obtener_matriz("SAS.gds", snp_sel, nombres_finales)

# Configuramos el panel (2 filas, 3 columnas)
par(mfrow = c(2, 3), mar = c(2, 2, 4, 2))

# Función para graficar mostrando los nombres de las variantes
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

# Dibujamos las 5 poblaciones
dibujar_con_nombres(matrix_eur, "EUR (MAF > 0.25)")
dibujar_con_nombres(matrix_afr, "AFR")
dibujar_con_nombres(matrix_eas, "EAS")
dibujar_con_nombres(matrix_amr, "AMR")
dibujar_con_nombres(matrix_sas, "SAS")

# Volver a 1 solo gráfico
par(mfrow = c(1, 1))

# Creamos un PDF de los gráficos para que tenga una mejor calidad y no se corten los ejes:
CREACION DE UN PDF DONDE SE MUESTREN TODAS LAS GRAFICAS SIN QUE SE CORTEN #######


# 1. Definimos el nombre del archivo y el tamaño (grande para que quepa todo)
pdf("Comparativa_LD_UGT1A1_Completa.pdf", width = 15, height = 10)

# 2. Configuramos el panel (2 filas, 3 columnas) con márgenes más amplios
par(mfrow = c(2, 3), mar = c(4, 4, 6, 2))

# Función optimizada para que no se corten los rsID
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

# 3. Dibujamos cada una
dibujar_tfm(matrix_eur, "EUROPE (EUR)")
dibujar_tfm(matrix_afr, "AFRICA (AFR)")
dibujar_tfm(matrix_eas, "EAST ASIA (EAS)")
dibujar_tfm(matrix_amr, "ADMIXED AMERICAS (AMR)")
dibujar_tfm(matrix_sas, "SOUTH ASIA (SAS)")

# 4. Cerramos el PDF
dev.off()

## 4. ANALISIS DE DESEQUILIBRIO DE LIGAMIENTO Y FRECUENCIAS DE VARIANTES INMPORTANTES EN CLINICA: rs8175347, rs4148323, rs887829 y rs6717546. 
# 1. Filtrar nel archivo original vcf de 1000Genomes para cada poblacion, sin tener en cuenta un MAF. 
 bcftools view -S EUR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.EURcomplet.vcf.gz
bcftools view -S AFR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.AFRcomplet.vcf.gz
bcftools view -S AMR.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.AMRcomplet.vcf.gz
bcftools view -S EAS.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.EAScomplet.vcf.gz
bcftools view -S SAS.txt UGT1A1_1000G.vcf.gz -Oz -o UGT1A1.SAScomplet.vcf.gz

# Indexar
for pop in EUR AFR AMR EAS SAS
do
  bcftools index UGT1A1.${pop}complet.vcf.gz
done

# 2. Calcular el LD de cada poblacion
#Convertir todos los VCF a GDS para que los trabaje R. 
library(SNPRelate)

snpgdsVCF2GDS("UGT1A1.EURcomplet.vcf.gz", "EUR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.AFRcomplet.vcf.gz", "AFR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.AMRcomplet.vcf.gz", "AMR.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.EAScomplet.vcf.gz", "EAS.gds", method="copy.num.of.ref")
snpgdsVCF2GDS("UGT1A1.SAScomplet.vcf.gz", "SAS.gds", method="copy.num.of.ref")


# empezar con poblacion EUR
gen_eur <- snpgdsOpen("EUR.gds")
# ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_eur, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_eur, "snp.id"))
# ver si estan las posiciones que me interesan
mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel <- ids[pos %in% mis_posiciones]
ids_sel
length(ids_sel)
# Calculo de LD
ld <- snpgdsLDMat(gen_eur,
                  snp.id = ids_sel,
                  method = "r",
                  slide = 0)

mat_ld <- ld$LD^2

# Añadir etiquetas
labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld) <- rownames(mat_ld) <- labels
# Heatmap
library(corrplot)

corrplot(mat_ld,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Europea",
         mar = c(0,0,3,0))

snpgdsClose(gen_eur)
# AHORA VAMOS A HACERLO PARA EL RESTO DE POBLACIONES, Y PARA NO REPETIR CÓDIGO, AJUSTAMOS: 
gen_afr <- snpgdsOpen("AFR.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_afr, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_afr, "snp.id"))
# ver si estan las posiciones que me interesan
mis_posiciones <- c(
  234668879,
  234668570, 
  234669144,
  234682119
)
ids_sel_afr <- ids[pos %in% mis_posiciones]
ids_sel_afr
length(ids_sel)
# Calculo de LD
ld_afr <- snpgdsLDMat(gen_afr,
                  snp.id = ids_sel_afr,
                  method = "r",
                  slide = 0)

mat_ld_afr <- ld_afr$LD^2

# Añadir etiquetas
labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_afr) <- rownames(mat_ld_afr) <- labels
# Heatmap
library(corrplot)

corrplot(mat_ld_afr,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Africana",
         mar = c(0,0,3,0))

snpgdsClose(gen_afr)

# POBLACION AMERICANA
gen_amr <- snpgdsOpen("AMR.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_amr, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_amr, "snp.id"))
# ver si estan las posiciones que me interesan
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

# Añadir etiquetas
labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_amr) <- rownames(mat_ld_amr) <- labels
# Heatmap
library(corrplot)

corrplot(mat_ld_amr,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Americana",
         mar = c(0,0,3,0))

snpgdsClose(gen_amr)

# POBLACION SAS
gen_sas <- snpgdsOpen("SAS.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_sas, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_sas, "snp.id"))
# ver si estan las posiciones que me interesan
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

# Añadir etiquetas
labels <- paste0("chr2:", pos[pos %in% mis_posiciones])

colnames(mat_ld_sas) <- rownames(mat_ld_sas) <- labels
# Heatmap
library(corrplot)

corrplot(mat_ld_sas,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.cex = 0.8,
         main = "DL de Variantes relevantes en clínica en población Sur Asiática",
         mar = c(0,0,3,0))

snpgdsClose(gen_sas)

# POBLACION ESTE DE ASIA- EAS
gen_eas <- snpgdsOpen("EAS.gds")
#ver que las posiciones eraan dentro del archivo
pos <- read.gdsn(index.gdsn(gen_eas, "snp.position"))
range(pos) #comprobar que tiene toda la region de UGT1A1
ids <- read.gdsn(index.gdsn(gen_eas, "snp.id"))
# ver si estan las posiciones que me interesan
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

# Vamos a guardar en un PDF todas las figuras que hemos cread. Para eso, abrimos unn PDF que vamos a darle un título y luego volvemos a cargar todo el código anterior y finalmente cerramos con devoff()
pdf("LD_UGT1A1_VariantesRelevantesClinica_xPoblaciones.pdf",
    width = 8,
    height = 6)
    
dev.off()



