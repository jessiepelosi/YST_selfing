# Plotting Linkage Disequalibrium Decay Across Genome

To plot the decay of linkage disequilibirum across the genome, first a **.vcf** files for each sample population is needed. We used vcftools to create the neccesary files from a single vcf containing
all sample populations with population maps containing indiviudals in each population. 
```
for file in *_samples; do vcftools --vcf california.maf0.01.biallelic.25mis.recode.vcf --keep $file --recode --out "$file".california.maf0.01.biallelic.25mis ; done 
```
We used plink to calculate linkage disequilibrium. In order for plink to preform the calculatiions it needs **.bed, .bim, and .fam files**, which can also be created with **plink**.
```
# Make bed, bim, fam files: 
plink --vcf samples.california.maf0.01.biallelic.25mis.recode.vcf --make-bed --double-id --allow-extra-chr --out samples   
```
With **.bed, .bim, and .fam** files now created, plink is able to calculate linkage disequilibrium.
```
# Calculate LD
plink --bfile samples --double-id --allow-extra-chr --set-missing-var-ids @:# --r2 gz --ld-window 100 --ld-window-kb 1000 --ld-window-r2 0 --out samples_ld 
```
To make the linkage decay data able to be easily plotted, we implemented additional python script sourced from **https://speciationgenomics.github.io/ld_decay/** which was modified. This script generates
files ending in **.ld_decay_bins**.
```
python ld_decay_calc.py -i samples_ld.ld.gz -o samples_Scaffold_1na
```
To plot this data in order to view the linkage decay curve, we read our files ending in **.ld_decay_bins** into **Rstudio**. Packages used include **ggplot2** and **readr**. Finally, we used **ggplot2** to generate a plot with linkage disequilibrium decay. We displayed ld deacy for each chromosome individually using **facet_wrap**. 
