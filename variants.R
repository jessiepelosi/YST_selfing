###############
# Classifying variants in each California population of YST
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated 7/25/2025
###############

library(dplyr)
library(ggplot2)

####################### FUNCTIONS, READING AND CLEANING FILES ###################


#the following functions are used to trim dataframes to only include wanted information and join multiple dataframes together


trim_df <- function(variant, df.sum, popname){
  df.variant <- df.sum %>% filter(V5==variant)
  
  df.variant$popname.freq <- df.variant$freq
  
  df.variant.trim <- df.variant %>% select(V1, V2, V3, V4, V5, popname.freq) %>% distinct()
  df.variant.trim <- df.variant.trim[order(df.variant.trim[, 'V2']),]
  df.variant.trim <- df.variant.trim[!duplicated(df.variant.trim$V2),]
  df.variant.trim <- data.frame(df.variant.trim)
  colnames(df.variant.trim) <- c("V1", "V2", "V3", "V4", "V5", popname)
  return(df.variant.trim)
}

join_dfs_full <- function(...) {
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  xxx = full_join(...)
  return(xxx)
}

join_dfs_inner <- function(...) {
  df1 = list(...)[[1]]
  df2 = list(...)[[2]]
  xxx = inner_join(..., by=c('V1', 'V2', 'V3', 'V4', 'V5'))
  return(xxx)
}

############# the following sections cleaned each population dataframe

###### population Col

col <- read.delim("col.cleaned.txt", header = F, sep = "\t", na.strings = ".")
col$V6 <- as.numeric(col$V6)
col$V7 <- as.numeric(col$V7)
col.sum <- cbind(col, rowSums(col[,6:31], na.rm = T))
col.sum <- col.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(col[, 6:31], na.rm = T)`))))
col.sum$freq = col.sum$`rowSums(col[, 6:31], na.rm = T)`/col.sum$num

###### population Des

des <- read.delim("des.cleaned.txt", header = F, sep = "\t", na.strings = ".")
des$V6 <- as.numeric(des$V6)
des$V7 <- as.numeric(des$V7)

des.sum <- cbind(des, rowSums(des[,6:37], na.rm = T))
des.sum <- des.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(des[, 6:37], na.rm = T)`))))
des.sum$freq = des.sum$`rowSums(des[, 6:37], na.rm = T)`/des.sum$num

####### population Dia 

dia <- read.delim("dia.cleaned.txt", header = F, sep = "\t", na.strings = ".")
dia$V6 <- as.numeric(dia$V6)
dia$V7 <- as.numeric(dia$V7)
dia.sum <- cbind(dia, rowSums(dia[,6:43], na.rm = T))
dia.sum <- dia.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(dia[, 6:43], na.rm = T)`))))
dia.sum$freq = dia.sum$`rowSums(dia[, 6:43], na.rm = T)`/dia.sum$num

####### population Gil 

gil <- read.delim("gil.cleaned.txt", header = F, sep = "\t", na.strings = ".")
gil$V6 <- as.numeric(gil$V6)
gil$V7 <- as.numeric(gil$V7)

gil.sum <- cbind(gil, rowSums(gil[,6:45], na.rm = T))
gil.sum <- gil.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(gil[, 6:45], na.rm = T)`))))
gil.sum$freq = gil.sum$`rowSums(gil[, 6:45], na.rm = T)`/gil.sum$num

####### population Gol

gol <- read.delim("gol.cleaned.txt", header = F, sep = "\t", na.strings = ".")
gol$V6 <- as.numeric(gol$V6)
gol$V7 <- as.numeric(gol$V7)
gol.sum <- cbind(gol, rowSums(gol[,6:41], na.rm = T))
gol.sum <- gol.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(gol[, 6:41], na.rm = T)`))))
gol.sum$freq = gol.sum$`rowSums(gol[, 6:41], na.rm = T)`/gol.sum$num

####### population Leb

leb <- read.delim("leb.cleaned.txt", header = F, sep = "\t", na.strings = ".")
leb$V6 <- as.numeric(leb$V6)
leb$V7 <- as.numeric(leb$V7)
leb.sum <- cbind(leb, rowSums(leb[,6:43], na.rm = T))
leb.sum <- leb.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(leb[, 6:43], na.rm = T)`))))
leb.sum$freq = leb.sum$`rowSums(leb[, 6:43], na.rm = T)`/leb.sum$num

####### population Nee

nee<- read.delim("nee.cleaned.txt", header = F, sep = "\t", na.strings = ".")
nee$V6 <- as.numeric(nee$V6)
nee$V7 <- as.numeric(nee$V7)
nee.sum <- cbind(nee, rowSums(nee[,6:17], na.rm = T))
nee.sum <- nee.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(nee[, 6:17], na.rm = T)`))))
nee.sum$freq = nee.sum$`rowSums(nee[, 6:17], na.rm = T)`/nee.sum$num

####### population Orc

orc <- read.delim("orc.cleaned.txt", header = F, sep = "\t", na.strings = ".")
orc$V6 <- as.numeric(orc$V6)
orc$V7 <- as.numeric(orc$V7)
orc.sum <- cbind(orc, rowSums(orc[,6:15], na.rm = T))
orc.sum <- orc.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(orc[, 6:15], na.rm = T)`))))
orc.sum$freq = orc.sum$`rowSums(orc[, 6:15], na.rm = T)`/orc.sum$num

####### population Oro

oro <- read.delim("oro.cleaned.txt", header = F, sep = "\t", na.strings = ".")
oro$V6 <- as.numeric(oro$V6)
oro$V7 <- as.numeric(oro$V7)
oro.sum <- cbind(oro, rowSums(oro[,6:39], na.rm = T))
oro.sum <- oro.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(oro[, 6:39], na.rm = T)`))))
oro.sum$freq = oro.sum$`rowSums(oro[, 6:39], na.rm = T)`/oro.sum$num

####### population Rb

rb <- read.delim("rb.cleaned.txt", header = F, sep = "\t", na.strings = ".")
rb$V6 <- as.numeric(rb$V6)
rb$V7 <- as.numeric(rb$V7)
rb.sum <- cbind(rb, rowSums(rb[,6:49], na.rm = T))
rb.sum <- rb.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                        -`rowSums(rb[, 6:49], na.rm = T)`))))
rb.sum$freq = rb.sum$`rowSums(rb[, 6:49], na.rm = T)`/rb.sum$num

####### population Res

res <- read.delim("res.cleaned.txt", header = F, sep = "\t", na.strings = ".")
res$V6 <- as.numeric(res$V6)
res$V7 <- as.numeric(res$V7)
res.sum <- cbind(res, rowSums(res[,6:23], na.rm = T))
res.sum <- res.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(res[, 6:23], na.rm = T)`))))
res.sum$freq = res.sum$`rowSums(res[, 6:23], na.rm = T)`/res.sum$num

####### population Sie

sie <- read.delim("sie.cleaned.txt", header = F, sep = "\t", na.strings = ".")
sie$V6 <- as.numeric(sie$V6)
sie$V7 <- as.numeric(sie$V7)
sie.sum <- cbind(sie, rowSums(sie[,6:41], na.rm = T))
sie.sum <- sie.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(sie[, 6:41], na.rm = T)`))))
sie.sum$freq = sie.sum$`rowSums(sie[, 6:41], na.rm = T)`/sie.sum$num

####### population Tri

tri <- read.delim("tri.cleaned.txt", header = F, sep = "\t", na.strings = ".")
tri$V6 <- as.numeric(tri$V6)
tri$V7 <- as.numeric(tri$V7)
tri.sum <- cbind(tri, rowSums(tri[,6:37], na.rm = T))
tri.sum <- tri.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(tri[, 6:37], na.rm = T)`))))
tri.sum$freq = tri.sum$`rowSums(tri[, 6:37], na.rm = T)`/tri.sum$num

####### population Uki

uki <- read.delim("uki.cleaned.txt", header = F, sep = "\t", na.strings = ".")
uki$V6 <- as.numeric(uki$V6)
uki$V7 <- as.numeric(uki$V7)
uki.sum <- cbind(uki, rowSums(uki[,6:33], na.rm = T))
uki.sum <- uki.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(uki[, 6:33], na.rm = T)`))))
uki.sum$freq = uki.sum$`rowSums(uki[, 6:33], na.rm = T)`/uki.sum$num

####### population Vet

vet <- read.delim("vet.cleaned.txt", header = F, sep = "\t", na.strings = ".")
vet$V6 <- as.numeric(vet$V6)
vet$V7 <- as.numeric(vet$V7)
vet.sum <- cbind(vet, rowSums(vet[,6:37], na.rm = T))
vet.sum <- vet.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(vet[, 6:37], na.rm = T)`))))
vet.sum$freq = vet.sum$`rowSums(vet[, 6:37], na.rm = T)`/vet.sum$num

####### population Wat

wat <- read.delim("wat.cleaned.txt", header = F, sep = "\t", na.strings = ".")
wat$V6 <- as.numeric(wat$V6)
wat$V7 <- as.numeric(wat$V7)
wat.sum <- cbind(wat, rowSums(wat[,6:17], na.rm = T))
wat.sum <- wat.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(wat[, 6:17], na.rm = T)`))))
wat.sum$freq = wat.sum$`rowSums(wat[, 6:17], na.rm = T)`/wat.sum$num

####### population Yre

yre <- read.delim("yre.cleaned.txt", header = F, sep = "\t", na.strings = ".")
yre$V6 <- as.numeric(yre$V6)
yre$V7 <- as.numeric(yre$V7)
yre.sum <- cbind(yre, rowSums(yre[,6:35], na.rm = T))
yre.sum <- yre.sum %>% mutate(num = rowSums(!is.na(select(., -V1, -V2, -V3, -V4, -V5, 
                                                          -`rowSums(yre[, 6:35], na.rm = T)`))))
yre.sum$freq = yre.sum$`rowSums(yre[, 6:35], na.rm = T)`/yre.sum$num


################################## FILTERING FOR MISSENSE VARIANTS ################################

col.missense.trim <- trim_df("missense_variant", col.sum, "Col")
des.missense.trim <- trim_df("missense_variant", des.sum, "Des")
dia.missense.trim <- trim_df("missense_variant", dia.sum, "Dia")
gil.missense.trim <- trim_df("missense_variant", gil.sum, "Gil")
gol.missense.trim <- trim_df("missense_variant", gol.sum, "Gol")
leb.missense.trim <- trim_df("missense_variant", leb.sum, "Leb")
nee.missense.trim <- trim_df("missense_variant", nee.sum, 'Nee')
orc.missense.trim <- trim_df("missense_variant", orc.sum, "Orc")
oro.missense.trim <- trim_df("missense_variant", oro.sum, "Oro")
rb.missense.trim <- trim_df("missense_variant", rb.sum, "Rb")
res.missense.trim <- trim_df("missense_variant", res.sum, "Res")
sie.missense.trim <- trim_df("missense_variant", sie.sum, "Sie")
tri.missense.trim <- trim_df("missense_variant", tri.sum, "Tri")
uki.missense.trim <- trim_df("missense_variant", uki.sum, "Uki")
vet.missense.trim <- trim_df("missense_variant", vet.sum, "Vet")
wat.missense.trim <- trim_df("missense_variant", wat.sum, "Wat")
yre.missense.trim <- trim_df("missense_variant", yre.sum, "Yre")


missense.popslist <- list(col.missense.trim, des.missense.trim, dia.missense.trim, gil.missense.trim,
                  gol.missense.trim, leb.missense.trim, nee.missense.trim, orc.missense.trim,
                  oro.missense.trim, rb.missense.trim, res.missense.trim, sie.missense.trim, 
                  tri.missense.trim, uki.missense.trim, vet.missense.trim, wat.missense.trim, 
                  yre.missense.trim)


missense_ij <- Reduce(join_dfs_inner, missense.popslist)

pivot_missense_ij <- pivot_longer(missense_ij, cols=6:22)

missense_pivot_summary_ij <- pivot_missense_ij %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

mean(missense_pivot_summary_ij$mean)

ggplot(missense_pivot_summary_ij, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Missense Variants in YST California Populations') +
  xlab("Population") +
  ylab("Mean Frequency")

#####################################   SYNONYMOUS VARIANTS ####################################

col.synonymous.trim <- trim_df("synonymous_variant", col.sum, "Col")
des.synonymous.trim <- trim_df("synonymous_variant", des.sum, "Des")
dia.synonymous.trim <- trim_df("synonymous_variant", dia.sum, "Dia")
gil.synonymous.trim <- trim_df("synonymous_variant", gil.sum, "Gil")
gol.synonymous.trim <- trim_df("synonymous_variant", gol.sum, "Gol")
leb.synonymous.trim <- trim_df("synonymous_variant", leb.sum, "Leb")
nee.synonymous.trim <- trim_df("synonymous_variant", nee.sum, "Nee")
orc.synonymous.trim <- trim_df("synonymous_variant", orc.sum, "Orc")
oro.synonymous.trim <- trim_df("synonymous_variant", oro.sum, "Oro")
rb.synonymous.trim <- trim_df("synonymous_variant", rb.sum, "Rb")
res.synonymous.trim <- trim_df("synonymous_variant", res.sum, "Res")
sie.synonymous.trim <- trim_df("synonymous_variant", sie.sum, "Sie")
tri.synonymous.trim <- trim_df("synonymous_variant", tri.sum, "Tri")
uki.synonymous.trim <- trim_df("synonymous_variant", uki.sum, "Uki")
vet.synonymous.trim <- trim_df("synonymous_variant", vet.sum, "Vet")
wat.synonymous.trim <- trim_df("synonymous_variant", wat.sum, "Wat")
yre.synonymous.trim <- trim_df("synonymous_variant", yre.sum, "Yre")

synonymous.popslist <- list(col.synonymous.trim, des.synonymous.trim, dia.synonymous.trim,
                         gil.synonymous.trim, gol.synonymous.trim, leb.synonymous.trim,
                         nee.synonymous.trim, orc.synonymous.trim, oro.synonymous.trim,
                         rb.synonymous.trim, res.synonymous.trim, sie.synonymous.trim,
                         tri.synonymous.trim, uki.synonymous.trim, vet.synonymous.trim,
                         wat.synonymous.trim, yre.synonymous.trim)

synonymous_compfreq_ij <- Reduce(join_dfs_inner, synonymous.popslist)

pivot_syn_ij <- pivot_longer(synonymous_compfreq_ij, cols=6:22)

syn_pivot_summary_ij <- pivot_syn_ij %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(syn_pivot_summary_ij, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Synonymous Variants in YST California Populations') +
  xlab('Population') +
  ylab('Mean Frequency')

################################# INTERGENIC REGIONS ##################################
col.intergenic.trim <- trim_df("intergenic_region", col.sum, "Col")
des.intergenic.trim <- trim_df("intergenic_region", des.sum, "Des")
dia.intergenic.trim <- trim_df("intergenic_region", dia.sum, "Dia")
gil.intergenic.trim <- trim_df("intergenic_region", gil.sum, "Gil")
gol.intergenic.trim <- trim_df("intergenic_region", gol.sum, "Gol")
leb.intergenic.trim <- trim_df("intergenic_region", leb.sum, "Leb")
nee.intergenic.trim <- trim_df("intergenic_region", nee.sum, "Nee")
orc.intergenic.trim <- trim_df("intergenic_region", orc.sum, "Orc")
oro.intergenic.trim <- trim_df("intergenic_region", oro.sum, "Oro")
rb.intergenic.trim <- trim_df("intergenic_region", rb.sum, "Rb")
res.intergenic.trim <- trim_df("intergenic_region", res.sum, "Res")
sie.intergenic.trim <- trim_df("intergenic_region", sie.sum, "Sie")
tri.intergenic.trim <- trim_df("intergenic_region", tri.sum, "Tri")
uki.intergenic.trim <- trim_df("intergenic_region", uki.sum, "Uki")
vet.intergenic.trim <- trim_df("intergenic_region", vet.sum, "Vet")
wat.intergenic.trim <- trim_df("intergenic_region", wat.sum, "Wat")
yre.intergenic.trim <- trim_df("intergenic_region", yre.sum, "Yre")


intergenic.poplist <- list(col.intergenic.trim, des.intergenic.trim, dia.intergenic.trim,
                           gil.intergenic.trim, gol.intergenic.trim, leb.intergenic.trim,
                           nee.intergenic.trim, orc.intergenic.trim, oro.intergenic.trim,
                           rb.intergenic.trim, res.intergenic.trim, sie.intergenic.trim,
                           tri.intergenic.trim, uki.intergenic.trim, vet.intergenic.trim,
                           wat.intergenic.trim, yre.intergenic.trim)

intergenic_compfreq_ij <- Reduce(join_dfs_inner, intergenic.poplist)

pivot_intergenic_ij <- pivot_longer(intergenic_compfreq_ij, cols=6:22)

intergenic_pivot_summary_ij <- pivot_intergenic_ij %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(syn_pivot_summary_ij, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Intergenic Regions in YST California Populations') +
  xlab('Population') + 
  ylab('Mean Frequency') 

################################# INTRON VARIANTS ##################################

col.introns.trim <- trim_df("intron_variant", col.sum, "Col")
des.introns.trim <- trim_df("intron_variant", des.sum, "Des")
dia.introns.trim <- trim_df("intron_variant", dia.sum, "Dia")
gil.introns.trim <- trim_df("intron_variant", gil.sum, "Gil")
gol.introns.trim <- trim_df("intron_variant", gol.sum, "Gol")
leb.introns.trim <- trim_df("intron_variant", leb.sum, "Leb")
nee.introns.trim <- trim_df("intron_variant", nee.sum, "Nee")
orc.introns.trim <- trim_df("intron_variant", orc.sum, "Orc")
oro.introns.trim <- trim_df("intron_variant", oro.sum, "Oro")
rb.introns.trim <- trim_df("intron_variant", rb.sum, "Rb")
res.introns.trim <- trim_df("intron_variant", res.sum, "Res")
sie.introns.trim <- trim_df("intron_variant", sie.sum, "Sie")
tri.introns.trim <- trim_df("intron_variant", tri.sum, "Tri")
uki.introns.trim <- trim_df("intron_variant", uki.sum, "Uki")
vet.introns.trim <- trim_df("intron_variant", vet.sum, "Vet")
wat.introns.trim <- trim_df("intron_variant", wat.sum, "Wat")
yre.introns.trim <- trim_df("intron_variant", yre.sum, "Yre")

introns.poplist <- list(col.introns.trim, des.introns.trim, dia.introns.trim, gil.introns.trim,
                        gol.introns.trim, leb.introns.trim, nee.introns.trim, orc.introns.trim,
                        oro.introns.trim, rb.introns.trim, res.introns.trim, sie.introns.trim,
                        tri.introns.trim, uki.introns.trim, vet.introns.trim, wat.introns.trim,
                        yre.introns.trim)

intron_compfreq_ij <- Reduce(join_dfs_inner, introns.poplist)

pivot_intron_ij <- pivot_longer(intron_compfreq_ij, cols=6:22)

intron_pivot_summary_ij <- pivot_intron_ij %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(intron_pivot_summary_ij, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Intron Regions in YST California Populations') +
  xlab('Population') + ylab('Mean Frequency') 

##################################### STOP-GAINED VARIANTS ###############################

col.stopgained.trim <- trim_df("stop_gained", col.sum, "Col")
des.stopgained.trim <- trim_df("stop_gained", des.sum, "Des")
dia.stopgained.trim <- trim_df("stop_gained", dia.sum, "Dia")
gil.stopgained.trim <- trim_df("stop_gained", gil.sum, "Gil")
gol.stopgained.trim <- trim_df("stop_gained", gol.sum, "Gol")
leb.stopgained.trim <- trim_df("stop_gained", leb.sum, "Leb")
nee.stopgained.trim <- trim_df("stop_gained", nee.sum, "Nee")
orc.stopgained.trim <- trim_df("stop_gained", orc.sum, "Orc")
oro.stopgained.trim <- trim_df("stop_gained", oro.sum, "Oro")
rb.stopgained.trim <- trim_df("stop_gained", rb.sum, "Rb")
res.stopgained.trim <- trim_df("stop_gained", res.sum, "Res")
sie.stopgained.trim <- trim_df("stop_gained", sie.sum, "Sie")
tri.stopgained.trim <- trim_df("stop_gained", tri.sum, "Tri")
uki.stopgained.trim <- trim_df("stop_gained", uki.sum, "Uki")
vet.stopgained.trim <- trim_df("stop_gained", vet.sum, "Vet")
wat.stopgained.trim <- trim_df("stop_gained", wat.sum, "Wat")
yre.stopgained.trim <- trim_df("stop_gained", yre.sum, "Yre")

stopgained.poplist <- list(col.stopgained.trim, des.stopgained.trim, dia.stopgained.trim, 
                           gil.stopgained.trim, gol.stopgained.trim, leb.stopgained.trim, 
                           nee.stopgained.trim, orc.stopgained.trim, oro.stopgained.trim, 
                           rb.stopgained.trim, res.stopgained.trim, sie.stopgained.trim,
                           tri.stopgained.trim, uki.stopgained.trim, vet.stopgained.trim, 
                           wat.stopgained.trim, yre.stopgained.trim)

stopgained_compfreq_ij <- Reduce(join_dfs_inner, stopgained.poplist)

pivot_stopgained_ij <- pivot_longer(stopgained_compfreq_ij, cols=6:22)

stopgained_pivot_summary_ij <- pivot_stopgained_ij %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(stopgained_pivot_summary_ij, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Stop-gained Instances in YST California Populations') +
  xlab('Population') + ylab('Mean Frequency')

################################### START-LOST VARIANTS ################################

col.startlost <- trim_df("start_lost", col.sum, "Col")
des.startlost <- trim_df("start_lost", des.sum, "Des")
dia.startlost <- trim_df("start_lost", dia.sum, "Dia")
gil.startlost <- trim_df("start_lost", gil.sum, "Gil")
gol.startlost <- trim_df("start_lost", gol.sum, "Gol")
leb.startlost <- trim_df("start_lost", leb.sum, "Leb")
nee.startlost <- trim_df("start_lost", nee.sum, "Nee")
orc.startlost <- trim_df("start_lost", orc.sum, "Orc")
oro.startlost <- trim_df("start_lost", oro.sum, "Oro")
rb.startlost <- trim_df("start_lost", rb.sum, "Rb")
res.startlost <- trim_df("start_lost", res.sum, "Res")
sie.startlost <- trim_df("start_lost", sie.sum, "Sie")
tri.startlost <- trim_df("start_lost", tri.sum, "Tri")
uki.startlost <- trim_df("start_lost", uki.sum, "Uki")
vet.startlost <- trim_df("start_lost", vet.sum, "Vet")
wat.startlost <- trim_df("start_lost", wat.sum, "Wat")
yre.startlost <- trim_df("start_lost", yre.sum, "Yre")

startlost.list <- list(col.startlost, des.startlost, dia.startlost, gil.startlost,
                       gol.startlost, leb.startlost, nee.startlost, orc.startlost,
                       oro.startlost, rb.startlost, res.startlost, sie.startlost,
                       tri.startlost, uki.startlost, vet.startlost, wat.startlost, yre.startlost)

startlost.compfreq <- Reduce(join_dfs_inner, startlost.list)

pivot.startlost.compfreq <- pivot_longer(startlost.compfreq, cols=6:22)

startlost.summary <- pivot.startlost.compfreq %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(startlost.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle(" Start Lost in YST California Populations") +
  xlab('Population') + ylab('Mean Frequency')


#################################### SPLICE REGION VARIANTS #############################

col.srv <- trim_df("splice_region_variant", col.sum, "Col")
des.srv <- trim_df("splice_region_variant", des.sum, "Des")
dia.srv <- trim_df("splice_region_variant", dia.sum, "Dia")
gil.srv <- trim_df("splice_region_variant", gil.sum, "Gil")
gol.srv <- trim_df("splice_region_variant", gol.sum, "Gol")
leb.srv <- trim_df("splice_region_variant", leb.sum, "Leb")
nee.srv <- trim_df("splice_region_variant", nee.sum, "Nee")
orc.srv <- trim_df("splice_region_variant", orc.sum, "Orc")
oro.srv <- trim_df("splice_region_variant", oro.sum, "Oro")
rb.srv <- trim_df("splice_region_variant", rb.sum, "Rb")
res.srv <- trim_df("splice_region_variant", res.sum, "Res")
sie.srv <- trim_df("splice_region_variant", sie.sum, "Sie")
tri.srv <- trim_df("splice_region_variant", tri.sum, "Tri")
uki.srv <- trim_df("splice_region_variant", uki.sum, "Uki")
vet.srv <- trim_df("splice_region_variant", vet.sum, "Vet")
wat.srv <- trim_df("splice_region_variant", wat.sum, "Wat")
yre.srv <- trim_df("splice_region_variant", yre.sum, "Yre")

srv.list <- list(col.srv, des.srv, dia.srv, gil.srv, gol.srv, leb.srv, nee.srv, orc.srv, oro.srv,
                 rb.srv, res.srv, sie.srv, tri.srv, uki.srv, vet.srv, wat.srv, yre.srv)

srv.compfreq <- Reduce(join_dfs_inner, srv.list)

pivot.srv.compfreq <- pivot_longer(srv.compfreq, cols=6:22)

srv.summary <- pivot.srv.compfreq %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(srv.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle("Splice Region Variants in YST California Populations") +
  xlab('Population') + ylab('Mean Frequency')


####################### MISSENSE VARIANTS AND SPLICE REGION VARIANTS #######################

col.mis.splice <- trim_df("missense_variant&splice_region_variant", col.sum, "Col")
des.mis.splice <- trim_df("missense_variant&splice_region_variant", des.sum, "Des")
dia.mis.splice <- trim_df("missense_variant&splice_region_variant", dia.sum, "Dia")
gil.mis.splice <- trim_df("missense_variant&splice_region_variant", gil.sum, "Gil")
gol.mis.splice <- trim_df("missense_variant&splice_region_variant", gol.sum, "Gol")
leb.mis.splice <- trim_df("missense_variant&splice_region_variant", leb.sum, "Leb")
nee.mis.splice <- trim_df("missense_variant&splice_region_variant", nee.sum, "Nee")
orc.mis.splice <- trim_df("missense_variant&splice_region_variant", orc.sum, "Orc")
oro.mis.splice <- trim_df("missense_variant&splice_region_variant", oro.sum, "Oro")
rb.mis.splice <- trim_df("missense_variant&splice_region_variant", rb.sum, "Rb")
res.mis.splice <- trim_df("missense_variant&splice_region_variant", res.sum, "Res")
sie.mis.splice <- trim_df("missense_variant&splice_region_variant", sie.sum, "Sie")
tri.mis.splice <- trim_df("missense_variant&splice_region_variant", tri.sum, "Tri")
uki.mis.splice <- trim_df("missense_variant&splice_region_variant", uki.sum, "Uki")
vet.mis.splice <- trim_df("missense_variant&splice_region_variant", vet.sum, "Vet")
wat.mis.splice <- trim_df("missense_variant&splice_region_variant", wat.sum, "Wat")
yre.mis.splice <- trim_df("missense_variant&splice_region_variant", yre.sum, "Yre")

mis.splice.list <- list(col.mis.splice, des.mis.splice, dia.mis.splice, gil.mis.splice,
                        gol.mis.splice, leb.mis.splice, nee.mis.splice, orc.mis.splice,
                        oro.mis.splice, rb.mis.splice, res.mis.splice, sie.mis.splice,
                        tri.mis.splice, uki.mis.splice, vet.mis.splice, wat.mis.splice,
                        yre.mis.splice)

mis.splice.compfreq <- Reduce(join_dfs_inner, mis.splice.list)

pivot_mis.splice <- pivot_longer(mis.splice.compfreq, cols=6:22)

mis.splice.summary <- pivot_mis.splice %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(mis.splice.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Missense and Splice Region Variants in YST California Populations') +
  xlab('Population') + ylab('Mean Frequency')

###################### SYNONYMOUS VARIANTS AND SPLICE REGION VARIANTS ##################

col.syn.splice <- trim_df("splice_region_variant&synonymous_variant", col.sum, "Col")
des.syn.splice <- trim_df("splice_region_variant&synonymous_variant", des.sum, "Des")
dia.syn.splice <- trim_df("splice_region_variant&synonymous_variant", dia.sum, "Dia")
gil.syn.splice <- trim_df("splice_region_variant&synonymous_variant", gil.sum, "Gil")
gol.syn.splice <- trim_df("splice_region_variant&synonymous_variant", gol.sum, "Gol")
leb.syn.splice <- trim_df("splice_region_variant&synonymous_variant", leb.sum, "Leb")
nee.syn.splice <- trim_df("splice_region_variant&synonymous_variant", nee.sum, "Nee")
orc.syn.splice <- trim_df("splice_region_variant&synonymous_variant", orc.sum, "Orc")
oro.syn.splice <- trim_df("splice_region_variant&synonymous_variant", oro.sum, "Oro")
rb.syn.splice <- trim_df("splice_region_variant&synonymous_variant", rb.sum, "Rb")
res.syn.splice <- trim_df("splice_region_variant&synonymous_variant", res.sum, "Res")
sie.syn.splice <- trim_df("splice_region_variant&synonymous_variant", sie.sum, "Sie")
tri.syn.splice <- trim_df("splice_region_variant&synonymous_variant", tri.sum, "Tri")
uki.syn.splice <- trim_df("splice_region_variant&synonymous_variant", uki.sum, "Uki")
vet.syn.splice <- trim_df("splice_region_variant&synonymous_variant", vet.sum, "Vet")
wat.syn.splice <- trim_df("splice_region_variant&synonymous_variant", wat.sum, "Wat")
yre.syn.splice <- trim_df("splice_region_variant&synonymous_variant", yre.sum, "Yre")

syn.splice.list <- list(col.syn.splice, des.syn.splice, dia.syn.splice, gil.syn.splice,
                        gol.syn.splice, leb.syn.splice, nee.syn.splice, orc.syn.splice,
                        oro.syn.splice, rb.syn.splice, res.syn.splice, sie.syn.splice,
                        tri.syn.splice, uki.syn.splice, vet.syn.splice, wat.syn.splice,
                        yre.syn.splice)

syn.splice.compfreq <- Reduce(join_dfs_inner, syn.splice.list)

pivot_syn.splice <- pivot_longer(syn.splice.compfreq, cols=6:22)

syn.splice.summary <- pivot_syn.splice %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(syn.splice.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Synonymous and Splice Region Variants in YST California Populations') +
  xlab('Population') + ylab('Mean Frequency')

######################### UPSTREAM GENE VARIANTS AND INTERGENIC REGIONS ##################

col.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", col.sum, "Col")
des.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", des.sum, "Des")
dia.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", dia.sum, "Dia")
gil.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", gil.sum, "Gil")
gol.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", gol.sum, "Gol")
leb.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", leb.sum, "Leb")
nee.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", nee.sum, "Nee")
orc.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", orc.sum, "Orc")
oro.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", oro.sum, "Oro")
rb.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", rb.sum, "Rb")
res.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", res.sum, "Res")
sie.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", sie.sum, "Sie")
tri.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", tri.sum, "Tri")
uki.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", uki.sum, "Uki")
vet.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", vet.sum, "Vet")
wat.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", wat.sum, "Wat")
yre.upstream.intergenic <- trim_df("upstream_gene_variant_intergenic_region", yre.sum, "Yre")

upstream.intergenic.list <- list(col.upstream.intergenic, des.upstream.intergenic, 
                                 dia.upstream.intergenic, gil.upstream.intergenic,
                                 gol.upstream.intergenic, leb.upstream.intergenic, 
                                 nee.upstream.intergenic, orc.upstream.intergenic,
                                 oro.upstream.intergenic, rb.upstream.intergenic, 
                                 res.upstream.intergenic, sie.upstream.intergenic,
                                 tri.upstream.intergenic, uki.upstream.intergenic, 
                                 vet.upstream.intergenic, wat.upstream.intergenic,
                                 yre.upstream.intergenic)

upstream.intergenic.compfreq <- Reduce(join_dfs_inner, upstream.intergenic.list)

pivot_upstream.intergenic.compfreq <- pivot_longer(upstream.intergenic.compfreq, cols=6:22)

upstream.intergenic.summary <- pivot_upstream.intergenic.compfreq %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(upstream.intergenic.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle('Upstream Gene Variants and Intergenic Regions in YST California Populations') +
  xlab('Population') + ylab('Mean Frequency')


##################################### 3' UTR VARIANTS ##################################

col.3primeUTR <- trim_df("3_prime_UTR_variant", col.sum, "Col")
des.3primeUTR <- trim_df("3_prime_UTR_variant", des.sum, "Des")
dia.3primeUTR <- trim_df("3_prime_UTR_variant", dia.sum, "Dia")
gil.3primeUTR <- trim_df("3_prime_UTR_variant", gil.sum, "Gil")
gol.3primeUTR <- trim_df("3_prime_UTR_variant", gol.sum, "Gol")
leb.3primeUTR <- trim_df("3_prime_UTR_variant", leb.sum, "Leb")
nee.3primeUTR <- trim_df("3_prime_UTR_variant", nee.sum, "Nee")
orc.3primeUTR <- trim_df("3_prime_UTR_variant", orc.sum, "Orc")
oro.3primeUTR <- trim_df("3_prime_UTR_variant", oro.sum, "Oro")
rb.3primeUTR <- trim_df("3_prime_UTR_variant", rb.sum, "Rb")
res.3primeUTR <- trim_df("3_prime_UTR_variant", res.sum, "Res")
sie.3primeUTR <- trim_df("3_prime_UTR_variant", sie.sum, "Sie")
tri.3primeUTR <- trim_df("3_prime_UTR_variant", tri.sum, "Tri")
uki.3primeUTR <- trim_df("3_prime_UTR_variant", uki.sum, "Uki")
vet.3primeUTR <- trim_df("3_prime_UTR_variant", vet.sum, "Vet")
wat.3primeUTR <- trim_df("3_prime_UTR_variant", wat.sum, "Wat")
yre.3primeUTR <- trim_df("3_prime_UTR_variant", yre.sum, "Yre")

utr3prime.list <- list(col.3primeUTR, des.3primeUTR, dia.3primeUTR, gil.3primeUTR, gol.3primeUTR,
                       leb.3primeUTR, nee.3primeUTR, orc.3primeUTR, oro.3primeUTR, rb.3primeUTR,
                       res.3primeUTR, sie.3primeUTR, tri.3primeUTR, uki.3primeUTR, vet.3primeUTR,
                       wat.3primeUTR, yre.3primeUTR)

utr3prime.compfreq <- Reduce(join_dfs_inner, utr3prime.list)

pivot.utr3prime.compfreq <- pivot_longer(utr3prime.compfreq, cols=6:22)

utr3prime.summary <- pivot.utr3prime.compfreq %>%
  group_by(name) %>%
  summarize(mean=mean(value, na.rm=TRUE), se=sd(value, na.rm=TRUE)/sqrt(length(value)))

ggplot(utr3prime.summary, mapping=aes(x=reorder(name, mean), y=mean)) +
  geom_point() + theme_bw() +
  geom_errorbar(aes(ymin=mean-se, ymax=mean+se)) +
  ggtitle(" 3' UTR Variants in YST California Populations") +
  xlab('Population') + ylab("Mean Frequency")
