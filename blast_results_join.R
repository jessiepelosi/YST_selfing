library(dplyr)

#des
des_blast <- read.delim("des_blast_results",header=F)
colnames(des_blast)<- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart','qend', 'sstart', 'send', 'evalue', 'bitscore')

#wat
wat_blast <- read.delim("wat_blast_results",header=F)
colnames(wat_blast)<- c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart','qend', 'sstart', 'send', 'evalue', 'bitscore')

#protein function db
protID <- read.delim(sep = "\t", "proteinID_function_db_clean", header=F)
colnames(protID)<- c('sseqid', 'function')

#Create functional annotations for DES
des_blast_with_func <- inner_join(protID,des_blast)
write.table(des_blast_with_func,"des_blast_with_func", row.names = F, quote = F, col.names = F, sep = "\t")

#Create functional annotations for WAT
wat_blast_with_func <- inner_join(protID,wat_blast)
write.table(wat_blast_with_func,"wat_blast_with_func", row.names = F, quote = F, col.names = F, sep = "\t")
