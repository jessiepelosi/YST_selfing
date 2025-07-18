# Population differentiation (F<sub>ST</sub>)
We used **vcftools** to calculate F<sub>ST</sub> in order to compare relatedness of our sample populations. 

```
vcftools
--vcf california.maf0.01.biallelic.25mis.recode.vcf
	--fst-window-size 1000000
	--fst-window-step 50000
	--weir-fst-pop keepsamples.noneeorcdeswat.txt
	--weir-fst-pop des_samples_list.txt
	--keep keepsamples.noneeorcdeswat.txt
	--keep des_samples_list.txt
	--out des_fst_out
```
