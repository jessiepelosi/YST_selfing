#############################
#FST PI windows for California Invasion (1mb)
#Taylor Curry, Abby Pearse, Jessie Pelosi
#Dlugosch lab
#July 31, 2025
############################

library(dplyr)
library(ggplot2)
library(gridExtra)

###Fst###
#read in DES and process Fst data
fst_des_outcrossing <- read.delim("des_fst_out.windowed.weir.fst")
fst_des_outcrossing <- mutate(fst_des_outcrossing,BIN_MID=(BIN_END-BIN_START)/2 + BIN_START)
fst_des_outcrossing$Z <- (fst_des_outcrossing$WEIGHTED_FST - mean(fst_des_outcrossing$WEIGHTED_FST))/sd(fst_des_outcrossing$WEIGHTED_FST)
fst_des_outcrossing$pval_fst <- pnorm(fst_des_outcrossing$Z, lower.tail = F)
fst_des_outcrossing_sig <- fst_des_outcrossing %>% filter(pval_fst<0.05)

#read in WAT and process Fst data
fst_wat_outcrossing <- read.delim("wat_fst_out.windowed.weir.fst")
fst_wat_outcrossing <- mutate(fst_wat_outcrossing,BIN_MID=(BIN_END-BIN_START)/2 + BIN_START)
fst_wat_outcrossing$Z <- (fst_wat_outcrossing$WEIGHTED_FST - mean(fst_wat_outcrossing$WEIGHTED_FST))/sd(fst_wat_outcrossing$WEIGHTED_FST)
fst_wat_outcrossing$pval_fst <- pnorm(fst_wat_outcrossing$Z, lower.tail = F)
fst_wat_outcrossing_sig <- fst_wat_outcrossing %>% filter(pval_fst<0.05)

###Pi###
##read in OUTCROSSING pops. and process pi data 
pi_outcrossing <- read.delim("1mb.cali.noneeorcdeswat.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
pi_outcrossing <- mutate(pi_outcrossing,BIN_MID=(BIN_END-BIN_START)/2 + BIN_START)

##read in DES pops. and process pi data 
pi_des <- read.delim("1mb.cali.des.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
pi_des <- mutate(pi_des,BIN_MID=(BIN_END-BIN_START)/2 + BIN_START)

#join DES and OUTCROSSING data and find ratio
pi_data_outcrossing_des <- inner_join(x=pi_des,y=pi_outcrossing,by=c("CHROM","BIN_MID"))
pi_data_outcrossing_des <- mutate(pi_data_outcrossing_des, ratio=PI.x/PI.y)

#Calculate log, zscore, pvalue for the DES and OUTCROSSING ratio
pi_data_outcrossing_des$logpi <- log2(pi_data_outcrossing_des$ratio)
pi_data_outcrossing_des$zscore <- (pi_data_outcrossing_des$logpi-mean(pi_data_outcrossing_des$logpi))/sd(pi_data_outcrossing_des$logpi)
pi_data_outcrossing_des$pval_pi <- pnorm(pi_data_outcrossing_des$zscore)
pi_data_outcrossing_des_sig <- pi_data_outcrossing_des %>% filter(pval_pi<0.05)

##read in WAT pops. and process pi data 
pi_wat <- read.delim("1mb.cali.wat.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
pi_wat <- mutate(pi_wat,BIN_MID=(BIN_END-BIN_START)/2 + BIN_START)

#join WAT and OUTCROSSING data and find ratio
pi_data_outcrossing_wat <- inner_join(x=pi_wat,y=pi_outcrossing,by=c("CHROM","BIN_MID"))
pi_data_outcrossing_wat <- mutate(pi_data_outcrossing_wat, ratio=PI.x/PI.y)

#Calculate log, zscore, pvalue for the WAT and OUTCROSSING ratio
pi_data_outcrossing_wat$logpi <- log2(pi_data_outcrossing_wat$ratio)
pi_data_outcrossing_wat$zscore <- (pi_data_outcrossing_wat$logpi-mean(pi_data_outcrossing_wat$logpi))/sd(pi_data_outcrossing_wat$logpi)
pi_data_outcrossing_wat$pval_pi <- pnorm(pi_data_outcrossing_wat$zscore)
pi_data_outcrossing_wat_sig <- pi_data_outcrossing_wat %>% filter(pval_pi<0.05)

#combine DES Fst and Pi data, select important data for finding sig genes
des_fst_pi_sig <- inner_join(pi_data_outcrossing_des_sig,fst_des_outcrossing_sig,
                             by = c("CHROM",'BIN_MID'))
des_fst_pi_sig_col123 <- des_fst_pi_sig %>% select("CHROM","BIN_START.x","BIN_END.x")
write.table(des_fst_pi_sig_col123,"des_fstpi_sig_windows",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
            )

#combine WAT Fst and Pi data, select important data for finding sig genes
wat_fst_pi_sig <- inner_join(pi_data_outcrossing_wat_sig,fst_wat_outcrossing_sig,
                             by = c("CHROM",'BIN_MID'))

wat_fst_pi_sig_col123 <- wat_fst_pi_sig %>% select("CHROM","BIN_START.x","BIN_END.x")
write.table(wat_fst_pi_sig_col123,"wat_fstpi_sig_windows",
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t"
)

#Calculate window overlaps
overlaps <- inner_join(des_fst_pi_sig_col123,wat_fst_pi_sig_col123)
