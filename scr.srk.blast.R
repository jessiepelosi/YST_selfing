###############################

# finding genes corresponding with SRK/SCR hits in DES and WAT
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated 7/28/2025

##############################

library(ggplot2)
library(dplyr)

######## finding hit with highest e-value per protein ID

blastp_srk <- read.delim("protein_srk_blasted.txt", header=FALSE)
colnames(blastp_srk) <-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                     'qend', 'sstart', 'send', 'evalue', 'bitscore')

best_hit_srk <- blastp_srk %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number()==1)

blastp_scr <- read.delim("scr_protein.blasted.txt", header=FALSE)
colnames(blastp_scr) <- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart',
                          'qend', 'sstart', 'send', 'evalue', 'bitscore')

best_hit_scr <- blastp_scr %>% 
  group_by(qseqid) %>% 
  arrange(evalue) %>% 
  filter(row_number()==1)

####### reading in .gff file of YST, genes with corresponding protein ID

cleaned.genes <- read.csv("cleaned.genomicgff.csv", header=FALSE)
cleaned.genes <- cleaned.genes %>% select(V1, V4, V5, V9)
colnames(cleaned.genes) <- c("chrom", "startpos", "endpos", "geneid")
cleaned.genes$geneid <- gsub("ID=gene-", "", cleaned.genes$geneid)

protein.names <- read.delim("protein_names_cleaned.txt", sep="", header=FALSE)
protein.names <- protein.names %>% select(V1, V4)
colnames(protein.names) <- c("protein_id", "geneid")

genes_pos <- inner_join(cleaned.genes, protein.names, by="geneid")
genes_pos$protein_id <- gsub(">", "", genes_pos$protein_id)
colnames(genes_pos) <- c("chrom", "startpos", "endpos", "geneid", "qseqid")


#filtering SRK hits on each chromosome

hits_srk_positions <- inner_join(best_hit_srk, genes_pos)
srk_locs <- hits_srk_positions %>% select("qseqid", "startpos", "endpos", "chrom")
srk_locs$gene <- "SRK"

srk.chrom1 <- srk_locs %>% filter(chrom=="CM058040.1")
srk.chrom2 <- srk_locs %>% filter(chrom=="CM058041.1")
srk.chrom3 <- srk_locs %>% filter(chrom=="CM058042.1")
srk.chrom4 <- srk_locs %>% filter(chrom=="CM058043.1")
srk.chrom5 <- srk_locs %>% filter(chrom=="CM058044.1")
srk.chrom6 <- srk_locs %>% filter(chrom=="CM058045.1")
srk.chrom7 <- srk_locs %>% filter(chrom=="CM058046.1")
srk.chrom8 <- srk_locs %>% filter(chrom=="CM058047.1")
 

## filtering out SCR hits located on each chromosome

hits_scr_positions <- inner_join(best_hit_scr, genes_pos, by="qseqid")
scr_locs <- hits_scr_positions %>% select("qseqid", "startpos", "endpos", "chrom")
scr_locs$gene <- "SCR"

scr.chrom1 <- scr_locs %>% filter(chrom=="CM058040.1")
scr.chrom2 <- scr_locs %>% filter(chrom=="CM058041.1")
scr.chrom3 <- scr_locs %>% filter(chrom=="CM058042.1")
scr.chrom4 <- scr_locs %>% filter(chrom=="CM058043.1")
#no hits on chrom 5
scr.chrom6 <- scr_locs %>% filter(chrom=="CM058045.1")
scr.chrom7 <- scr_locs %>% filter(chrom=="CM058046.1")
scr.chrom8 <- scr_locs %>% filter(chrom=="CM058047.1")

##### combining SRK and SCR hits into one dataframe per chromosome
chrom1 <- bind_rows(srk.chrom1, scr.chrom1)
chrom1 <- arrange(chrom1, startpos)

chrom2 <- bind_rows(srk.chrom2, scr.chrom2)
chrom2 <- arrange(chrom2, startpos)

chrom3 <- bind_rows(srk.chrom3, scr.chrom3)
chrom3 <- arrange(chrom3, startpos)

chrom4 <- bind_rows(srk.chrom4, scr.chrom4)
chrom4 <- arrange(chrom4, startpos)

chrom6 <- bind_rows(srk.chrom6, scr.chrom6)
chrom6 <- arrange(chrom6, startpos)

chrom7 <- bind_rows(srk.chrom7, scr.chrom7)
chrom7 <- arrange(chrom7, startpos)

chrom8 <- bind_rows(srk.chrom8, scr.chrom8)
chrom8 <- arrange(chrom8, startpos)

####### finding distances between start positions of srk/scr genes

distances <- function(srk.df, scr.df){
  clean_scr_df <- na.omit(scr.df)
  clean_srk_df <- na.omit(srk.df)
  scr_qseqids <- character()
  srk_qseqids <- character()
  differences <- numeric()
  for (i in 1:nrow(scr.df)){
    for (j in 1:nrow(srk.df)){
      scr_current_qseqid <- clean_scr_df$qseqid[i] 
      srk_current_qseqid <- clean_srk_df$qseqid[j] 
      current_difference <- clean_scr_df$startpos[i] - clean_srk_df$startpos[j]
      scr_qseqids <- c(scr_qseqids, scr_current_qseqid)
      srk_qseqids <- c(srk_qseqids, srk_current_qseqid)
      differences <- c(differences, current_difference)
    }
  }
  result_df <- data.frame(scr_qseqids, srk_qseqids, differences)
  colnames(result_df) <- c('scr_qseqids', 'protein_id', 'differences')
  result_df$abs_dis=abs(result_df$differences)
  result_df <- arrange(result_df, abs_dis)
  result_df <- result_df %>% filter(abs_dis != 0)
  return(result_df)
}

### finding genes within 100,000 bp of one another, write out file with genes/protein ID within 100kb

# chromosome 1
chrom1_distances <- distances(srk.chrom1, scr.chrom1)
chrom1_genesin100kb <- chrom1_distances %>% filter (abs_dis < 100000)
chrom1_gene_protein <- inner_join(chrom1_genesin100kb, protein.names, by="protein_id")
chrom1_gene_protein <- chrom1_gene_protein %>% select(protein_id, differences, 
                                                      abs_dis, geneid, scr_qseqids)
write.table(chrom1_gene_protein, row.names = FALSE, "chrom1.bed", sep="\t", quote=FALSE)

#chromosome 2
chrom2_distances <- distances(srk.chrom2, scr.chrom2)
chrom2_genesin100kb <- chrom2_distances %>% filter (abs_dis < 100000)
chrom2_gene_protein <- inner_join(chrom2_genesin100kb, proteinnames.df, by="protein_id")
chrom2_gene_protein <- chrom2_gene_protein %>% select(protein_id, differences, 
                                                     abs_dis, geneid, scr_qseqids)
write.table(chrom2_gene_protein, row.names = FALSE, "chrom2.bed", sep="\t", quote=FALSE)


#chromosome 3 - no SRK/SCR genes within 100kb of one another
chrom3_distances <- distances(srk.chrom3, scr.chrom3)
chrom3_genesin100kb <- chrom3_distances %>% filter (abs_dis < 100000)

#chromosome 4
chrom4_distances <- distances(srk.chrom4, scr.chrom4)
chrom4_genesin100kb <- chrom4_distances %>% filter (abs_dis < 100000)
chrom4_gene_protein <- inner_join(chrom4_genesin100kb, proteinnames.df)
chrom4_gene_protein <- chrom4_gene_protein %>% select(protein_id, differences, 
                                                     abs_dis, geneid, scr_qseqids)
write.table(chrom4_gene_protein, row.names = FALSE, "chrom4.bed", sep="\t", quote=FALSE)

#chromosome 6
chrom6_distances <- distances(srk.chrom6, scr.chrom6)
chrom6_genesin100kb <- chrom6_distances %>% filter (abs_dis < 100000)
chrom6_gene_protein <- inner_join(chrom6_genesin100kb, proteinnames.df, by="protein_id")
chrom6_gene_protein <- chrom6_gene_protein %>% select(protein_id, differences, 
                                                     abs_dis, geneid, scr_qseqids)
write.table(chrom6_gene_protein, row.names = FALSE, "chrom6.bed", sep="\t", quote=FALSE)

#chromosome 7 - no SRK/SCR genes within 100kb of one another
chrom7_distances <- distances(srk.chrom7, scr.chrom7)
chrom7_genesin100kb <- chrom7_distances %>% filter (abs_dis < 100000)

#chromosome 8 - no SRK/SCR genes within 100kb of one another
chrom8_distances <- distances(srk.chrom8, scr.chrom8)
chrom8_genesin100kb <- chrom8_distances %>% filter (abs_dis < 100000)










