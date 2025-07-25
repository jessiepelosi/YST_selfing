###################
# ADMIXTURE plots for all California YST samples
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated 7/25/2025
###################

library(ggplot2)
library(devtools)
library(popkin)


# determining optimal K values based on CV error

cverror_cali <- read.delim("cverror_cali.out", sep="\t", header=FALSE)

ggplot(cverror_cali, mapping=aes(x=V1, y=V2)) + geom_line() + theme_bw()


# CV error curve suggested optimal K values are K=2 and K=3

Q2 <- read.delim("california_pca_structure.2.Q", sep="", header=FALSE)
Q3 <- read.delim("california_pca_structure.3.Q", sep="", header=FALSE)

loc <- read.delim("alphabetized.calisamples.txt", sep="", header=FALSE)
cali <- read.delim("california.pop", sep="\t", header=FALSE)
pop <- loc[loc$V1 %in% cali$V1, ]
pop <- loc$pop

#plotting K=2
plot_admix(Q2, col=c('lightblue', 'purple'), labs = loc2,
           labs_cex = 0.75, labs_sep=TRUE,labs_lwd = 3, labs_col="black", 
           labs_las = c(3), xlab_line=3) 

#plotting K=3
plot_admix(Q3, col=c('forestgreen', 'purple', 'lightblue'), labs = loc2,
           labs_cex = 0.75, labs_sep=TRUE,labs_lwd = 3, labs_col="black", 
           labs_las = c(3), xlab_line=3)
