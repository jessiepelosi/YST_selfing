# Site Frequency Spectra

We generated a site frequency spectrum for each of the 17 California populations to visualize distribution of unique alleles. I NEED MORE EXPLANATION I DON'T KNOW WHY WERE DOING THIS!!

<b> Installing and Running easySFS.py </b>

This easySFS.py python script was developed by (CITE THE SOURCE HERE). Instructions for installation come directly from this developer's page.

```
conda create -n easySFS
conda activate easySFS
conda install -c conda-forge numpy pandas scipy -y
git clone https://github.com/isaacovercast/easySFS.git
cd easySFS
chmod 777 easySFS.py
./easySFS.py
```

<b> Projecting Coordinates </b>

This program finds the number of alleles mapped to the number of sites for each population. Coordinates picked for each population aim to maximize the number of alleles (x) while having high number of sites.
Coordinates for each number of alleles per population can be viewed in a .preview file generated for each population.

```
./easySFS.py -i population1.maf0.01.biallelic.25mis.vcf -p samples_population1 --preview --window-bp 50000 > population1.preview
```
<b> Optimizing Coordinates </b>

After choosing the optimal coordinate for each population, the following code will be run. The --proj flag represents the number of alleles (x) from the optimized coordinate of the chosen population.
In this case, 14 will represent our optimized x-coordinate.
```
./easySFS.py -i population1.maf0.01.biallelic.25mis.recode.vcf -p samples_population1 --window-bp 50000 --proj 14 -o population1.output
```

After this process has been repeated for each population, output files can be downloaded to be visualized in R.

