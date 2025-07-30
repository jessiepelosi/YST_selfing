#########################

# Estimating nucleotide diversity for each California population
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated 7/28/2025

#########################

library(dplyr)
library(ggplot2)


#### COL #######
col_pi <- read.delim("col.windowed.pi")
col_mean_pi <- mean(col_pi$PI)
col_sterror_pi <- sd(col_pi$PI)/sqrt(length(col_pi$PI))

##### DES #####
des_pi <- read.delim("des.windowed.pi")
des_mean_pi <- mean(des_pi$PI)
des_sterror_pi <- sd(des_pi$PI)/sqrt(length(des_pi$PI))


#### DIA #####
dia_pi <- read.delim("dia.windowed.pi")
dia_mean_pi <- mean(dia_pi$PI)
dia_sterror_pi <- sd(dia_pi$PI)/sqrt(length(dia_pi$PI))

#### GIL ###
gil_pi <- read.delim("gil_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
gil_mean_pi <- mean(gil_pi$PI)
gil_sterror_pi <- sd(gil_pi$PI)/sqrt(length(gil_pi$PI))


#### GOL ####
gol_pi <- read.delim("gol_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
gol_mean_pi <- mean(gol_pi$PI)
gol_sterror_pi <- sd(gol_pi$PI)/sqrt(length(gol_pi$PI))

#### LEB ####
leb_pi <- read.delim("leb_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
leb_mean_pi <- mean(leb_pi$PI)
leb_sterror_pi <- sd(leb_pi$PI)/sqrt(length(leb_pi$PI))


### NEE ####
nee_pi <- read.delim("nee_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
nee_mean_pi <- mean(nee_pi$PI)
nee_sterror_pi <- sd(nee_pi$PI)/sqrt(length(nee_pi$PI))


### ORC ####
orc_pi <- read.delim("orc_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
orc_mean_pi <- mean(orc_pi$PI)
orc_sterror_pi <- sd(orc_pi$PI)/sqrt(length(orc_pi$PI))


### ORO ####
oro_pi <- read.delim("oro_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
oro_mean_pi <- mean(oro_pi$PI)
oro_sterror_pi <- sd(oro_pi$PI)/sqrt(length(oro_pi$PI))


### RB ##
rb_pi <- read.delim("rb_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
rb_mean_pi <- mean(rb_pi$PI)
rb_sterror_pi <- sd(rb_pi$PI)/sqrt(length(rb_pi$PI))


#### res
res_pi <- read.delim("res_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
res_mean_pi <- mean(res_pi$PI)
res_sterror_pi <- sd(res_pi$PI)/sqrt(length(res_pi$PI))


### sie
sie_pi <- read.delim("sie_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
sie_mean_pi <- mean(sie_pi$PI)
sie_sterror_pi <- sd(sie_pi$PI)/sqrt(length(sie_pi$PI))


## tri
tri_pi <- read.delim("tri_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
tri_mean_pi <- mean(tri_pi$PI)
tri_sterror_pi <- sd(tri_pi$PI)/sqrt(length(tri_pi$PI))


### uki
uki_pi <- read.delim("uki_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
uki_mean_pi <- mean(uki_pi$PI)
uki_sterror_pi <- sd(uki_pi$PI)/sqrt(length(uki_pi$PI))


## vet
vet_pi <- read.delim("vet_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
vet_mean_pi <- mean(vet_pi$PI)
vet_sterror_pi <- sd(vet_pi$PI)/sqrt(length(vet_pi$PI))


### wat
wat_pi <- read.delim("wat_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
wat_mean_pi <- mean(wat_pi$PI)
wat_sterror_pi <- sd(wat_pi$PI)/sqrt(length(wat_pi$PI))


### yre
yre_pi <- read.delim("yre_samples.cali.maf0.01.biallelic.25mis.recode.vcf.windowed.pi")
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
 
