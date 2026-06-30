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