################################################
#   MIK2 blastx results analysis               #
#   Taylor Curry, Abby Pearse, Jessie Pelosi   #
#   August 14, 2025                            #
################################################

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

#Read in blastx results (remove hits on JARXY in a txt editor first)
MIK2_blastx <- read_tsv("MIK2_ortho_blast_results")

#Sort and filter for best hits for each query
best_MIK2_blastx <- MIK2_blastx %>%
  group_by(qseqid) %>%
  arrange(evalue) %>%
  filter(row_number()==1) %>% 
  filter(evalue<1e-10)


#Clean data
clean_best_MIK2_blastx <- separate(
  data = best_MIK2_blastx,
  col = qseqid,
  into = c("chromosome", "source", "qseqid", "delete"),
  sep = "_"
)

#Hide useless columns
clean_best_MIK2_blastx$source <- NULL
clean_best_MIK2_blastx$delete <- NULL

#Read in key for protein to gene ID
gene_key <- read_tsv("ID_protein_gene_key.txt", col_names=F)
colnames(gene_key) <- c("qseqid", "geneID")

#Join gene key and clean blastx results
geneID_MIK2_blastx <- inner_join(clean_best_MIK2_blastx,gene_key)

#Save column with geneID to obtain positions from .gff with HPC
geneID_hpc<- as_data_frame(geneID_MIK2_blastx$geneID)
write.table(geneID_hpc,
            row.names = F,
            col.names = F,
            file = "GeneID_hpc.txt",
            quote = F
)

#Read in and clean GeneIDs with postions from .gff
gene_pos <- read_tsv("GeneID_positions", col_names = F)
gene_pos <- separate(
  data = gene_pos,
  col = X4,
  into = c("id"),
  sep = ";"
)

gene_pos <- separate(
  data = gene_pos,
  col =id,
  into = c("delete", "geneID"),
  sep = "-"
)

colnames(gene_pos) <- c("chromosome", "start", "end", "delete","geneID")
gene_pos$delete <- NULL
gene_pos$chromosome <- NULL

#Join df containing blastx results with geneID positions
fin_MIK2_blastx <- inner_join(geneID_MIK2_blastx,gene_pos, by = "geneID")

#Read in intersect from bedtools
MIK2_blast_intersect <- read_tsv("MIK2_bed_intersect", col_names = F)

#Remove negative values
MIK2_blast_intersect <- filter(MIK2_blast_intersect,X13 >= 0)

#Read in file from bedtools intersect DES fst pi windows and MIK2 blast
DES_MIK2_fstpi_int <- read_tsv("MIK2_fstpiwindows_bed_intersect", col_names = F)

#Separate out gene id
DES_MIK2_fstpi_int <- separate(
  data = DES_MIK2_fstpi_int,
  col = X9,
  into = c("delete", "GeneID"),
  sep = "-"
)
DES_MIK2_fstpi_int$delete <- NULL
DES_MIK2_fstpi_int <- separate(
  data = DES_MIK2_fstpi_int,
  col = GeneID,
  into = c("GeneID", "MISC"),
  sep = ";"
)

#Distinct genes in significant windows DES
DES_MIK2_fstpi_int <- distinct(DES_MIK2_fstpi_int,GeneID)

#Read in file from bedtools intersect WAT fst pi windows and MIK2 blast
wat_MIK2_fstpi_int <- read_tsv("wat_MIK2_fstpiwindows_bed_intersect", col_names = F)

#Separate out gene id
wat_MIK2_fstpi_int <- separate(
  data = wat_MIK2_fstpi_int,
  col = X9,
  into = c("delete", "GeneID"),
  sep = "-"
)
wat_MIK2_fstpi_int$delete <- NULL
wat_MIK2_fstpi_int <- separate(
  data = wat_MIK2_fstpi_int,
  col = GeneID,
  into = c("GeneID", "MISC"),
  sep = ";"
)

#Distinct genes in significant windows WAT
wat_MIK2_fstpi_int <- distinct(wat_MIK2_fstpi_int,GeneID)

