#################
# PCA of all California YST Samples
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated 7/25/2025
#################

library(ggplot2)
library(dplyr)

# reading in files from PLINK

pca <- read.delim("california_pca_structure.eigenvec",sep="", header=FALSE)
eigenval <- scan(paste("california_pca_structure.eigenval",sep=""))

pca <- pca[,-1]
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))


loc <- read.delim("alphabetized.calisamples.txt", sep="", header=FALSE)
cali <- read.delim("california.pop", sep="\t", header=FALSE)
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
ggplot(pca, mapping=aes(x=PC1, y=PC2, color=pop)) + geom_point() +
  theme_bw() + 
  scale_color_manual(values=c('DES'='hotpink', 'WAT'='blue')) +
  ggtitle("DES and WAT Clustering in YST California Populations")



