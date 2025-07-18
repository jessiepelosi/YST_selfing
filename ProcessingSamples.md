# Processing Samples 

Put in some fancy intro about processing our samples.
<b> We will use placeholder objects to represent our data. </b> 

<b> process_radtags </b>

The data used for genomic analyses were RADseq data, which includes barcodes in sequences. Our data needed to be formatted <i> in a way STACKS</i> (Catchen et al., 2013) would accept. Inset reasoning here I don't know it

```
awk '{print $2, $3, $1}' barcode_data > barcode_data_editedforStacks
```

We used process_radtags from the module STACKS to demultiplex DNA strands. This can only be done <i> after verifying barcodes and RAD cutsites</i>.

```
module load stacks
process_radtags -p barcode_data -P -b barcode_data_editedforStacks -o ./ -c -q -r --renz-1 pstI --renz-2 mseI --inline-inline 
```

<b> fastqc </b>

This step is performed to check the quality of each read. 

```
module load fastqc
for file in *.fq.gz; do fastqc $file;done
```

<b> Aligning reads to reference genome from NCBI </b>

To align reads to a reference genome, a genome needs to be downloaded onto a local computer from NCBI. This command will take the GCA accesion code(?) of the specific genome.

```
datasets download genome accession gca_number --include gff3,rtein,genome,seq-report
```

After locally obtaining the genome information, the reference genome needs to be indexed using BWA. This command will take the .fna file generated form the prior step.

```
module load bwa
bwa index GCA_number_genomic.fna
```

After this step, reads are able to be mapped to the reference genome. Our read data contained both forward and reverse reads, denoted as .1 and .2 in this example. Our 'sample' object represents a list of samples
we ran this command on to maximize efficieny. This step outputs the aligned reads as a .sam file.
```

module load bwa
for file in $(cat samples); do bwa mem GCA_number_genomic.fna "$file".1.fq.gz "$file".2.fq.gz > "$file".sam;done
```

<b> Generating .sam statistics </b>

After generating .sam files from the previous step, we wanted to obtain basic statistics about the alignment percentages for each read. 
```

module load samtools
for file in *.sam; do samtools flagstat "$file" > "$file".flagstat;done
```

<b> Converting .sam file to .bam file </b>

This step converts .sam files to .bam files, which are needed to the next step, as well as to improve efficiency when further processing data.

```
module load samtools
for file in $(cat samples); do samtools sort "$file".sam -o "$file".bam -O bam;done
```

<b> STACKS ref_map.pl </b>

This step requires .bam files of each read aligned to a reference genome, as well as a population map for all samples (column 1 sample name, column 2 population name, tab-delimited). In our case,
we wanted each population to be specified as a .vcf file. In this case, the tilde represents the absolute path to the file locations.

ASK JESSIE!!!

```
ref_map.pl --samples ~/bam.files --popmap population_map_of_samples -o ~/ref_pl.output -X "populations:--vcf"
```
