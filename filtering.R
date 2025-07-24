####################
# Filtering samples with less than 2500 loci
# Abby Pearse, Taylor Curry, Jessie Pelosi
# Last updated July 25th, 2025
####################

library(ggplot2)
library(dplyr)

gstacks <- read.csv("stacks_coverage_4refmap.csv", header=TRUE)

mean_nloci <- mean(gstacks$n_loci)
stdv_nloci <- sd(gstacks$n_loci)

gstacks.filtered2500 <- filter(gstacks, n_loci > 2500)
gstacks.filtered2500.mean <- mean(gstacks.filtered2500$n_loci)

ggplot(gstacks.filtered2500, aes(n_loci)) + geom_histogram() + 
  geom_vline(xintercept = mean_nloci, color='blue') +
  geom_vline(xintercept=gstacks.filtered2500.mean, color='red') +
  geom_vline(xintercept=2500, color='forestgreen') +
  xlim(c(0, 30000))
