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


