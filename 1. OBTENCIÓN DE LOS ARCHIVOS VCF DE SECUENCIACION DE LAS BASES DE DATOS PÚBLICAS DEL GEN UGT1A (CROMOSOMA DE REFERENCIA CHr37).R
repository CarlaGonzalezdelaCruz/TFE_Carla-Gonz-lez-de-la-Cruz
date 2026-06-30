
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