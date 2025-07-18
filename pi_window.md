#Nucleotide diversity (π)
We used **vcftools** calculate nucleotide diversity or π across set windows with a

```
vcftools --vcf des_samples.california.maf0.01.biallelic.25mis.recode.vcf --out 1mb.cali.des.maf0.01.biallelic.25mis.recode.vcf --window-pi 1000000 --window-pi-step 50000
```
