# Generating Population Structure of Samples

Before we can plot PCAs and ADMIXTURE plots to determine the population structure of samples in R, we need to process samples in HPC. 

<b> Filtering samples based on quality and specified thresholds </b>

Before we created any files, we made sure to filter out unique samples (no siblings) with low number of loci (n>2500) and the following criteria: 1) minor allele frequency (MAF) > 0.01 2) biallelic 3) missing less than 25% of data.
We then wrote the samples following the criteria into a seperate .vcf files. The "sample_list" represents a list of samples with no siblings that have greater than 2500 loci.

```
module load vcftools
vcftools --vcf samples.vcf --keep sample_list --maf 0.01 --min-allele 2 --max-allele 2 --max-missing 0.75 --out samples.maf0.01.biallelic.missing25.vcf
```

<b> Using PLINK to generate output files for PCA</b>

We used PLINK to generate the .eigenvector and .eigenvalue files necessary for PCA (why???????). This first step converts our filtered .vcf files into a .vcf file that is recognizable by PLINK itself and generates
.bed, .bim, .fam, .log, and .nosex files.

```
module load plink
plink --vcf samples.maf0.01.biallelic.missing25.vcf --make-bed --double-id --allow-extra-chr
```

This next step identifies linkage sites to perform linkage pruning, and notably generates a .prune.in and .prune.out file.

```
plink --bfile samples.maf0.01.biallelic.missing25.vcf --double-id --allow-extra-chr  --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out samples.maf0.01.biallelic.missing25.pca
```
Finally, we generate necessary files for PCA using this final step.

```
plink --bfile samples.maf0.01.biallelic.missing25.vcf --samples.maf0.01.biallelic.missing25.pca.prune.in --double-id --allow-extra-chr --set-missing-var-ids @:# --make-bed --pca --out samples.maf0.01.biallelic.missing25.pca.structure
```
```
