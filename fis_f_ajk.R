##########################
# Calculating FIS, F, and AJK relatedness for California YST populations
# Jessie Pelosi, Abby Pearse, Taylor Curry
# Last updated 7/28/2025
##########################

library(adegenet)
library(hierfstat)
library(vcfR)
library(ggplot2)
library(dplyr)

#read in data from VCF
vcf <- vcfR::read.vcfR("california.maf0.01.biallelic.25mis.recode.vcf")

# Convert VCF to genind 
gid <- vcfR2genind(vcf)

# Add population information for genind 
popmap <- read.delim("samples.hz.txt")
pop(gid) <- popmap$Population

# Separate gid into individual population geninds 
obj <- seppop(gid)

# Run basic.stats on each population genind 
get_fis <- function(stat) {
  stat.df <- as.data.frame(stat$overall)
  fis <- stat.df[9,]
  return(fis)
}

Col.stat <- basic.stats(obj$Col)
Col.fis <- get_fis(Col.stat)

Des.stat <- basic.stats(obj$Des)
Des.fis <- get_fis(Des.stat)

Dia.stat <- basic.stats(obj$Dia)
Dia.fis <- get_fis(Dia.stat)

Gil.stat <- basic.stats(obj$Gil)
Gil.fis <- get_fis(Gil.stat)

Gol.stat <- basic.stats(obj$Gol)
Gol.fis <- get_fis(Gol.stat)

Leb.stat <- basic.stats(obj$Leb)
Leb.fis <- get_fis(Leb.stat)

Nee.stat <- basic.stats(obj$Nee)
Nee.fis <- get_fis(Nee.stat)

Orc.stat <- basic.stats(obj$Orc)
Orc.fis <- get_fis(Orc.stat)

Oro.stat <- basic.stats(obj$Oro)
Oro.fis <- get_fis(Oro.stat)

Rb.stat <- basic.stats(obj$Rb)
Rb.fis <- get_fis(Rb.stat)

Res.stat <- basic.stats(obj$Res)
Res.fis <- get_fis(Res.stat)

Sie.stat <- basic.stats(obj$Sie)
Sie.fis <- get_fis(Sie.stat)

Tri.stat <- basic.stats(obj$Tri)
Tri.fis <- get_fis(Tri.stat)

Uki.stat <- basic.stats(obj$Uki)
Uki.fis <- get_fis(Uki.stat)

Vet.stat <- basic.stats(obj$Vet)
Vet.fis <- get_fis(Vet.stat)

Wat.stat <- basic.stats(obj$Wat)
Wat.fis <- get_fis(Wat.stat)

Yre.stat <- basic.stats(obj$Yre)
Yre.fis <- get_fis(Yre.stat)

fis <- c(Col.fis, Des.fis, Dia.fis, Gil.fis, Gol.fis, Leb.fis, Nee.fis, Orc.fis, Rb.fis, Res.fis, Sie.fis, 
         Tri.fis, Uki.fis, Vet.fis, Wat.fis, Yre.fis)

pops <- c("Col", "Des", "Dia", "Gil", "Gol", "Leb", "Nee", "Orc", "Rb", "Res", "Sie", "Tri", "Uki", "Vet", "Wat", "Yre")

fis_df <- as.data.frame(fis, pops)
fis_df$pop <- rownames(fis_df)

ggplot(data = fis_df, mapping = aes(x=reorder(pop, fis), y = fis)) +geom_point() +
  theme_bw() + ggtitle("Inbreeding Coefficient Fis from hierfstat") +
  xlab("Population") + ylab(bquote(F[is]))


##### F

vcftools_het <- read.delim("out.het")
vcftools_het$population <- popmap$Population
ggplot(data = vcftools_het, mapping = aes(x=reorder(population, F), y = F)) +geom_boxplot() +
  theme_bw() + ggtitle("Inbreeding Coefficient F from vcftools") + 
  xlab("Population") + ylab("Median F")


# Let's also look at relatedness as calculated by vcftools
relatedness <- read.delim("out.relatedness")
relatedness_pops <- read.delim("pop_comparisons_relatedness.txt", header = F)
relatedness$INDV1POP <- relatedness_pops$V1
relatedness$INDV2POP <- relatedness_pops$V2

t <- relatedness %>% 
  group_by(INDV1POP, INDV2POP) %>% 
  summarize(median = median(RELATEDNESS_AJK))

selfOnly_rel1 <- t %>% 
  filter(INDV1POP == INDV2POP)

ggplot(selfOnly_rel1, mapping = aes(x=reorder(INDV1POP, median), y = median)) + geom_point() + 
  theme_bw() + xlab("Population") + ylab("Median Within-Population Relatedness (Ajk Method)")

