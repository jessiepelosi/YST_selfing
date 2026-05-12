#########################
#
# Population genetic metrics for YST 
# Abby Pearse, Taylor Curry, Jessie Pelosi
#
#########################

library(dplyr)
library(ggplot2)
library(adegenet)
library(hierfstat)
library(vcfR)

##################
####### Pi #######
##################


#### COL #######
col_pi <- read.delim("data/col.windowed.pi")
col_mean_pi <- mean(col_pi$PI)
col_sterror_pi <- sd(col_pi$PI)/sqrt(length(col_pi$PI))

##### DES #####
des_pi <- read.delim("data/des.windowed.pi")
des_mean_pi <- mean(des_pi$PI)
des_sterror_pi <- sd(des_pi$PI)/sqrt(length(des_pi$PI))


#### DIA #####
dia_pi <- read.delim("data/dia.windowed.pi")
dia_mean_pi <- mean(dia_pi$PI)
dia_sterror_pi <- sd(dia_pi$PI)/sqrt(length(dia_pi$PI))

#### GIL ###
gil_pi <- read.delim("data/gil_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
gil_mean_pi <- mean(gil_pi$PI)
gil_sterror_pi <- sd(gil_pi$PI)/sqrt(length(gil_pi$PI))


#### GOL ####
gol_pi <- read.delim("data/gol_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
gol_mean_pi <- mean(gol_pi$PI)
gol_sterror_pi <- sd(gol_pi$PI)/sqrt(length(gol_pi$PI))

#### LEB ####
leb_pi <- read.delim("data/leb_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
leb_mean_pi <- mean(leb_pi$PI)
leb_sterror_pi <- sd(leb_pi$PI)/sqrt(length(leb_pi$PI))


### NEE ####
nee_pi <- read.delim("data/nee_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
nee_mean_pi <- mean(nee_pi$PI)
nee_sterror_pi <- sd(nee_pi$PI)/sqrt(length(nee_pi$PI))


### ORC ####
orc_pi <- read.delim("data/orc_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
orc_mean_pi <- mean(orc_pi$PI)
orc_sterror_pi <- sd(orc_pi$PI)/sqrt(length(orc_pi$PI))


### ORO ####
oro_pi <- read.delim("data/oro_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
oro_mean_pi <- mean(oro_pi$PI)
oro_sterror_pi <- sd(oro_pi$PI)/sqrt(length(oro_pi$PI))


### RB ##
rb_pi <- read.delim("data/rb_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
rb_mean_pi <- mean(rb_pi$PI)
rb_sterror_pi <- sd(rb_pi$PI)/sqrt(length(rb_pi$PI))


#### res
res_pi <- read.delim("data/res_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
res_mean_pi <- mean(res_pi$PI)
res_sterror_pi <- sd(res_pi$PI)/sqrt(length(res_pi$PI))


### sie
sie_pi <- read.delim("data/sie_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
sie_mean_pi <- mean(sie_pi$PI)
sie_sterror_pi <- sd(sie_pi$PI)/sqrt(length(sie_pi$PI))


## tri
tri_pi <- read.delim("data/tri_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
tri_mean_pi <- mean(tri_pi$PI)
tri_sterror_pi <- sd(tri_pi$PI)/sqrt(length(tri_pi$PI))


### uki
uki_pi <- read.delim("data/uki_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
uki_mean_pi <- mean(uki_pi$PI)
uki_sterror_pi <- sd(uki_pi$PI)/sqrt(length(uki_pi$PI))


## vet
vet_pi <- read.delim("data/vet_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
vet_mean_pi <- mean(vet_pi$PI)
vet_sterror_pi <- sd(vet_pi$PI)/sqrt(length(vet_pi$PI))


### wat
wat_pi <- read.delim("data/wat_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
wat_mean_pi <- mean(wat_pi$PI)
wat_sterror_pi <- sd(wat_pi$PI)/sqrt(length(wat_pi$PI))


### yre
yre_pi <- read.delim("data/yre_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
yre_mean_pi <- mean(yre_pi$PI)
yre_sterror_pi <- sd(yre_pi$PI)/sqrt(length(yre_pi$PI))


means <- c(col_mean_pi, des_mean_pi, dia_mean_pi, gil_mean_pi, gol_mean_pi, leb_mean_pi, nee_mean_pi, orc_mean_pi,
           oro_mean_pi, rb_mean_pi, res_mean_pi, sie_mean_pi, tri_mean_pi, uki_mean_pi, vet_mean_pi, wat_mean_pi,
           yre_mean_pi)

populations <- c("COL", "DES", "DIA", "GIL", "GOL", "LEB", "NEE", "ORC", "ORO", "RB", "RES", "SIE", "TRI", "UKI",
                 "VET", "WAT", "YRE")


sterrors <- c(col_sterror_pi, des_sterror_pi, dia_sterror_pi, gil_sterror_pi, gol_sterror_pi, leb_sterror_pi,
              nee_sterror_pi, orc_sterror_pi, oro_sterror_pi, rb_sterror_pi, res_sterror_pi, sie_sterror_pi,
              tri_sterror_pi, uki_sterror_pi, vet_sterror_pi, wat_sterror_pi, yre_sterror_pi)

pi.dataframe <- as.data.frame(means, populations)
pi.dataframe$StandardError <- sterrors
pi.dataframe$pop <- rownames(pi.dataframe)

ggplot(pi.dataframe, mapping=aes(x=reorder(pop, means), y=means)) + geom_point(size=3) + 
  geom_errorbar(aes(ymin=means-StandardError, ymax=means+StandardError), width=.2) +
  ggtitle("mean Pi") + 
  xlab("Population") + ylab("Mean Nucleotide Diversity (pi)") + theme_bw()


##### Additional metrics used VCF input 


#read in data from VCF
vcf <- vcfR::read.vcfR("data/california.maf0.01.biallelic.25mis.recode.vcf")

# Convert VCF to genind 
gid <- vcfR2genind(vcf)

# Add population information for genind 
popmap <- read.delim("data/samples.hz.txt")
pop(gid) <- popmap$Population

# Separate gid into individual population geninds 
obj <- seppop(gid)

# Run basic.stats on each population genind 
get_fis <- function(stat) {
  stat.df <- as.data.frame(stat$overall)
  fis <- stat.df[9,]
  return(fis)
}

#Heterozygosity is generated in basic.stats. To view both expected and observed heterozygosity, simply call the population.stat argument

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

vcftools_het <- read.delim("data/out.het")
vcftools_het$population <- popmap$Population
ggplot(data = vcftools_het, mapping = aes(x=reorder(population, F), y = F)) +geom_boxplot() +
  theme_bw() + ggtitle("Inbreeding Coefficient F from vcftools") + 
  xlab("Population") + ylab("Median F")


# Let's also look at relatedness as calculated by vcftools
relatedness <- read.delim("data/out.relatedness")
relatedness_pops <- read.delim("data/pop_comparisons_relatedness.txt", header = F)
relatedness$INDV1POP <- relatedness_pops$V1
relatedness$INDV2POP <- relatedness_pops$V2

t <- relatedness %>% 
  group_by(INDV1POP, INDV2POP) %>% 
  summarize(median = median(RELATEDNESS_AJK))

selfOnly_rel1 <- t %>% 
  filter(INDV1POP == INDV2POP)

ggplot(selfOnly_rel1, mapping = aes(x=reorder(INDV1POP, median), y = median)) + geom_point() + 
  theme_bw() + xlab("Population") + ylab("Median Within-Population Relatedness (Ajk Method)")


 
