#############################
## Population Structure
## Abby Pearse, Taylor Curry 
## Jessie Pelosi
############################

library(ggplot2)
library(dplyr)
library(popkin)

############################
###### ALL SPP + POP #######
############################

pca <- read.delim("data/allpops_allspp_pca_structure.eigenvec", sep="", header = F)
eigenval <- scan(paste("data/allpops_allspp_pca_structure.eigenval", sep =""))

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

spp <- c(rep("CESO",119), rep("CESO-DES",16), rep("CESO",383), rep("CESO-WAT", 6), rep("CESO",15), rep("CEME",3), rep("CENI", 5), rep("CEPA", 5))

pca <- tbl_df(data.frame(pca, spp))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
cumsum(pve$pve)

a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()


ggplot(pca, mapping=aes(x=PC1, y=PC2, color=spp)) + geom_point() +
  theme_bw() +
  ggtitle("Clustering of All California YST Populations")


# Omit all CENI and CEPA

pca_ceso_ceme <- filter(pca, spp == "CESO" | spp == "CEME" | spp == "CESO-DES" | spp == "CESO-WAT")
pca_ceso_ceme$spp <- factor(pca_ceso_ceme$spp, levels = c("CEME","CESO", "CESO-DES", "CESO-WAT"))
pca_ceso_ceme <- pca_ceso_ceme %>% arrange(spp)

ggplot(pca_ceso_ceme, mapping=aes(x=PC1, y=PC2, color=spp)) + geom_point(size = 3) +
  theme_bw() + theme(text = element_text(size = 15)) +
  scale_color_manual(values = c('CEME' = "#f0d25b", 'CESO' = "#8ad0ee", 'CESO-DES' ='#cd327b', 'CESO-WAT'= '#4169e1'), 
                     labels = c('CEME' = 'Centaurea melitensis', "CESO" = "Centaurea solsistialis", 
                                'CESO-DES' = 'Centaurea solsistialis - DES', 'CESO-WAT' = 'Centaurea solsistialis - WAT')) + 
  xlab(paste0("PC1 (", round(pve[1,2],2), "% Variation Expalined)")) + 
  ylab(paste0("PC2 (", round(pve[2,2],2), "% Variation Expalined)")) +
  theme(legend.title = element_blank())

ggsave("ceso_ceme.pdf", height = 6, width = 9)
ggsave("ceso_ceme.png", height = 6, width = 9, dpi = 300)


############################
###### YST CALI ONLY #######
############################

# reading in files from PLINK

pca <- read.delim("data/california_pca_structure.eigenvec",sep="", header=FALSE)
eigenval <- scan(paste("data/california_pca_structure.eigenval",sep=""))

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


loc <- read.delim("data/alphabetized.calisamples.txt", sep="", header=FALSE)
cali <- read.delim("data/california.pop", sep="\t", header=FALSE)
pop <- loc[loc$V1 %in% cali$V1, ]
pop <- pop$V2

pca <- tbl_df(data.frame(pca, pop))
pve <- data.frame(PC = 1:20, pve = eigenval/sum(eigenval)*100)
cumsum(pve$pve)

# plotting percent variance for each PC
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity")
a + ylab("Percentage variance explained") + theme_light()

# plotting PCA clustering for all California populations
ggplot(pca, mapping=aes(x=PC1, y=PC2, color=pop)) + geom_point() +
  theme_bw() +
  ggtitle("Clustering of All California YST Populations")


#plotting PCA clustering for all populations, but only coloring DES and WAT populations
ggplot(pca, mapping=aes(x=PC1, y=PC2, color=pop)) + geom_point(size =3) +
  scale_color_manual(values=c('DES'='#cd327b', 'WAT'='#4169e1', 
                              "COL" = '#d3d2d2', 'DIA' = '#d3d2d2',
                              'GIL' = '#d3d2d2', 'GOL' = '#d3d2d2',
                              'LEB' = '#d3d2d2', 'NEE' = '#d3d2d2',
                              'ORC' = '#d3d2d2', 'ORO' = '#d3d2d2',
                              'RB' = '#d3d2d2', 'RES' = '#d3d2d2',
                              'SIE' = '#d3d2d2', 'TRI' = '#d3d2d2',
                              'UKI' = '#d3d2d2', 'VET' = '#d3d2d2',
                              'YRE' = '#d3d2d2')) +
  xlab(paste0("PC1 (", round(pve[1,2],2), "% Variation Explained)")) +
  ylab(paste0("PC2 (", round(pve[2,2],2), "% Variation Expalined")) +
  theme_bw() + theme(text=element_text(size =15 ))
  

ggsave("ceso_Cali.pdf", height = 6, width = 7)
ggsave("ceso_Cali.png", height = 6, width = 7, dpi = 300)

## ADMIXTURE Results 

# determining optimal K values based on CV error

cverror_cali <- read.delim("data/cverror_cali.out", sep="\t", header=FALSE)

ggplot(cverror_cali, mapping=aes(x=V1, y=V2)) + geom_line() + theme_bw() +
  xlab("K") + ylab("Cross Validation Error") +
  theme(text = element_text(size = 15))


# CV error curve suggested optimal K values are K=2 and K=3

Q2 <- read.delim("data/california_pca_structure.2.Q", sep="", header=FALSE)
Q3 <- read.delim("data/california_pca_structure.3.Q", sep="", header=FALSE)

loc <- read.delim("data/alphabetized.calisamples.txt", sep="", header=FALSE)
cali <- read.delim("data/california.pop", sep="\t", header=FALSE)
pop <- loc[loc$V1 %in% cali$V1, ]
pop <- loc$pop

#plotting K=2
pdf(file = "cali.k2.pdf", width = 12, height = 4)
plot_admix(Q2, col=c('#d3d2d2', '#cd327b'), labs = loc$V2,
           labs_cex = 0.75, ylab_cex = 1, labs_sep=TRUE, labs_lwd = 1.5, labs_col="black", 
           labs_las = c(3), xlab_line=0.5, leg_omit = T, leg_width = 0, xlab = '') 
dev.off()

png(file = "cali.k2.png", width = 12, height = 4, units = "in", res = 300)
plot_admix(Q2, col=c('#d3d2d2', '#cd327b'), labs = loc$V2,
           labs_cex = 0.75, ylab_cex = 1, labs_sep=TRUE, labs_lwd = 1.5, labs_col="black", 
           labs_las = c(3), xlab_line=0.5, leg_omit = T, leg_width = 0, xlab = '') 
dev.off()


#plotting K=3
pdf(file = "cali.k3.pdf", width = 12, height = 4)
plot_admix(Q3, col=c('#4169e1', '#cd327b', '#d3d2d2'),labs = loc$V2,
           labs_cex = 0.75, ylab_cex = 1, labs_sep=TRUE, labs_lwd = 1.5, labs_col="black", 
           labs_las = c(3), xlab_line=0.5, leg_omit = T, leg_width = 0, xlab = '')
dev.off()

png(file = "cali.k3.png", width = 12, height = 4, units = "in", res = 300)
plot_admix(Q3, col=c('#4169e1', '#cd327b', '#d3d2d2'),labs = loc$V2,
           labs_cex = 0.75, ylab_cex = 1, labs_sep=TRUE, labs_lwd = 1.5, labs_col="black", 
           labs_las = c(3), xlab_line=0.5, leg_omit = T, leg_width = 0, xlab = '')
dev.off()

