################
# Site Frequency Spectra for California YST Populations
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated July 25th, 2025
################

library(dplyr)
library(ggplot2)
library(tidyr)
library(reshape2)

# after downloading output file, it is necessary to transpose the table in Excel or similar program.

######## population Col: coordinate (18, 1416)

col_sfs <- read.delim("col_sfs_transposed.txt")
col_sfs <- col_sfs %>% select(x, y)
col_sfs <- filter(col_sfs, y > 0)
col_sfs$Frequency <- (col_sfs$y/(2*1416))
colnames(col_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(col_sfs, aes(num.alleles, frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() + 
  ggtitle("Site Frequency Spectrum: Col") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Des: coordinate (14, 578)

des_sfs <- read.delim("des_transposed.txt")
des_sfs <- des_sfs %>% select(x,y)
des_sfs <- filter(des_sfs, y > 0)
des_sfs$frequency <- (des_sfs$y/(2*578))
colnames(des_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(des_sfs, aes(num.alleles, frequency)) + geom_point() +
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Des") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Dia: coordinate (26, 1562)

dia_sfs <- read.delim("dia_tranposed.txt")
dia_sfs <- dia_sfs %>% select(x,y)
dia_sfs <- filter(dia_sfs, y > 0)
dia_sfs$frequency <- (dia_sfs$y/(2 *1562))
colnames(dia_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(dia_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Dia") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")
  
######## population Gil: coordinate (28, 1611)

gil_sfs <- read.delim("gil_transposed.txt")
gil_sfs <- gil_sfs %>% select(x,y)
gil_sfs <- filter(gil_sfs, y > 0)
gil_sfs$frequency <- (gil_sfs$y/(2* 1611))
colnames(gil_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(gil_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Gil") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Gol: coordinate (24, 1485)

gol_sfs <- read.delim("gol_transposed.txt")
gol_sfs <- gol_sfs %>% select(x,y)
gol_sfs <- filter(gol_sfs, y > 0)
gol_sfs$frequency <- (gol_sfs$y/(2 *1485))
colnames(gol_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(gol_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Gol") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Leb: coordinate (20, 1279)

leb_sfs <- read.delim("leb_transposed.txt")
leb_sfs <- leb_sfs %>% select(x,y)
leb_sfs <- filter(leb_sfs, y > 0)
leb_sfs$frequency <- (leb_sfs$y/(2 *1279))
colnames(leb_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(leb_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Leb") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Nee: coordinate (10, 541)

nee_sfs <- read.delim("nee_transposed.txt")
nee_sfs <- nee_sfs %>% select(x,y)
nee_sfs <- filter(nee_sfs, y>0)
nee_sfs$frequency <- (nee_sfs$y/(2 * 541))
colnames(nee_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(nee_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Nee") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Orc: coordinate (6, 457)

orc_sfs <- read.delim("orc_transposed.txt")
orc_sfs <- orc_sfs %>% select(x,y)
orc_sfs <- filter(orc_sfs, y>0)
orc_sfs$frequency <- (orc_sfs$y/(2 * 541))
colnames(orc_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(orc_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Orc") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Oro: coordinate (24, 1473)

oro_sfs <- read.delim("oro_transposed.txt")
oro_sfs <- oro_sfs %>% select(x,y)
oro_sfs <- filter(oro_sfs, y>0)
oro_sfs$frequency <- (oro_sfs$y/(2 * 541))
colnames(oro_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(oro_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Oro") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Rb: coordinate (30, 1440)

rb_sfs <- read.delim("rb_transposed.txt")
rb_sfs <- rb_sfs %>% select(x,y)
rb_sfs <- filter(rb_sfs, y>0)
rb_sfs$frequency <- (rb_sfs$y/(2 * 1440))
colnames(rb_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(rb_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Rb") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Res: coordinate (12, 1128)

res_sfs <- read.delim("res_transposed.txt")
res_sfs <- res_sfs %>% select(x,y)
res_sfs <- filter(res_sfs, y>0)
res_sfs$frequency <- (res_sfs$y/(2 * 1128))
colnames(res_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(res_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Res") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Sie: coordinate (24, 1318)

sie_sfs <- read.delim("sie_transposed.txt")
sie_sfs <- sie_sfs %>% select(x,y)
sie_sfs <- filter(sie_sfs, y>0)
sie_sfs$frequency <- (sie_sfs$y/(2 * 1128))
colnames(sie_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(sie_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Sie") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Tri: coordinate (20, 1320)

tri_sfs <- read.delim("tri_transposed.txt")
tri_sfs <- tri_sfs %>% select(x,y)
tri_sfs <- filter(tri_sfs, y>0)
tri_sfs$frequency <- (tri_sfs$y/(2 * 1128))
colnames(tri_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(tri_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Tri") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Uki: coordinate (20, 1187)

uki_sfs <- read.delim("uki_transposed.txt")
uki_sfs <- uki_sfs %>% select(x,y)
uki_sfs <- filter(uki_sfs, y>0)
uki_sfs$frequency <- (uki_sfs$y/(2 * 1187))
colnames(uki_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(uki_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Uki") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Vet: coordinate (22, 1287)

vet_sfs <- read.delim("vet_transposed.txt")
vet_sfs <- vet_sfs %>% select(x,y)
vet_sfs <- filter(vet_sfs, y>0)
vet_sfs$frequency <- (vet_sfs$y/(2 * 1287))
colnames(vet_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(vet_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Vet") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Wat: coordinate (6, 388)

wat_sfs <- read.delim("wat_transposed.txt")
wat_sfs <- wat_sfs %>% select(x,y)
wat_sfs <- filter(wat_sfs, y>0)
wat_sfs$frequency <- (wat_sfs$y/(2 * 1287))
colnames(wat_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(wat_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Wat") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

######## population Yre: coordinate (22, 1402)

yre_sfs <- read.delim("yre_transposed.txt")
yre_sfs <- yre_sfs %>% select(x,y)
yre_sfs <- filter(yre_sfs, y>0)
yre_sfs$frequency <- (yre_sfs$y/(2 * 1287))
colnames(yre_sfs) <- c("num.alleles", "num.sites", "frequency")

ggplot(yre_sfs, aes(num.alleles,frequency)) + geom_point() + 
  geom_line() + 
  theme_bw() +
  ggtitle("Site Frequency Spectrum: Yre") +
  xlab("Number of Alleles") +
  ylab("Frequency of Alleles")

#overlapping maps
ggplot() + 
  geom_line(data=col_sfs, aes(num.alleles,frequency), color = "black") +
  geom_line(data=des_sfs, aes(num.alleles,frequency), color='red', linewidth=1) +
  geom_line(data=dia_sfs, aes(num.alleles,frequency), color='black') + 
  geom_line(data=gil_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=gol_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=leb_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=nee_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=orc_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=oro_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=rb_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=res_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=sie_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=tri_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=uki_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=vet_sfs, aes(num.alleles,frequency), color='black') +
  geom_line(data=wat_sfs, aes(num.alleles,frequency), color='royalblue', linewidth=1) +
  geom_line(data=yre_sfs, aes(num.alleles,frequency), color='black') +
  ggtitle("Site Frequency Spectra all YST California Populations") +
  xlab("Number of Alleles") + ylab("Frequency of Alleles") +
  theme_bw()





