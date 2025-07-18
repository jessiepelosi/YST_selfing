# Intersect Significant Windows (low pi, high F<sub>ST</sub>) with Protein IDs

We ran intersect on HPC to examine which proteins were found in significant windows of high F<sub>ST</sub> and low pi values (windows found in R, script in repo)
and the protein IDs of proteins found in the yellowstar thistle. After running this command, we were able to identify which proteins were present in our significsant
windows.

The first step of this program requires a .bed file - in our case, we needed to convert a .gff file into a .bed file.

```
module load bedops
gff2bed < GCA_number_genomic.gff > GCA_number_genomic.bed
```

After generating our .bed file, we were able to run intersect on HPC. It is important to note that this step requires a .bed file for the significant windows
as well, so we created a separate file from our filtered dataframe in R that had the 3 necessary columns of a .bed file (chromosome, start position, end position).

```
module load bedtools
bedtools intersect -a GCA_number_genomic.bed -b significant.windows.bed > genes.in.significantwindow
```

