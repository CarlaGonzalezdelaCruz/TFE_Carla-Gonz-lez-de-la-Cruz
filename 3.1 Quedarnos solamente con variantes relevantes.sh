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
