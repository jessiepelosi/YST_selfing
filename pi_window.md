# Nucleotide diversity (π)
We used **vcftools** calculate nucleotide diversity or π across set windows with a length of **1Mb** with a step size of **50Kb**.

```
vcftools --vcf samples.california.maf0.01.biallelic.25mis.recode.vcf --out 1mb.cali.samples.maf0.01.biallelic.25mis.recode.vcf --window-pi 1000000 --window-pi-step 50000
```
