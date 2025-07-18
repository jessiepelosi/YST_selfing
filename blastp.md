# BLASTp YST against Asteraceae and Arabidopsis

We utilized BLAST (Basic Local Alignment Search Tool) to identify potenital areas of the genome containing proteins believed to be associated with the
s-locus. 

First we made a database containing all the protein sequences from the Asteraceae and Arabadopsis reference genomes we selected.
```
makeblastdb
  -in proteins.faa
  -out protein_db
  -dbtype prot
```
Then we queried against the database we created with a selection of protein IDs obtained from our F<sub>ST</sub> and π intersect with a statistiaclly
significant p-value.

```
blastp
  -query protein_seqs.fq
  -db db_aster.arabid.fasta
  -evalue 1e-5
  -outfmt 6
  -num_threads 8
  -out blast_results
```
