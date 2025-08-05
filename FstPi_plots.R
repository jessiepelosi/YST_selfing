#############################
#Plotting FST PI windows for California Invasion (1mb)
#Taylor Curry, Abby Pearse, Jessie Pelosi
#Dlugosch lab
#July 31, 2025
############################

library(dplyr)
library(ggplot2)
library(gridExtra)

#processing data and plotting Fst & Pi for DES
#pi ratio 
pi_des_out_plots <- pi_data_outcrossing_des
pi_des_out_plots$CHROM <- gsub("Scaffold_","",pi_des_out_plots$CHROM)
pi_des_out_plots$chr <- factor(pi_des_out_plots$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
pi_des_out_plots_data <- pi_des_out_plots %>%
  arrange(chr, BIN_MID) %>%
  group_by(chr) %>%
  mutate(pos_max_this_chr = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(bp_add = lag(cumsum(pos_max_this_chr), default = 0)) %>%
  group_by(chr) %>%
  mutate(bp_cum = BIN_MID + bp_add - min(BIN_MID)) %>% 
  ungroup()

pi_axis_set <- pi_des_out_plots_data %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

des_plot_pi <- ggplot(pi_des_out_plots_data, aes(x = bp_cum, y = -log10(pval_pi)))+
  geom_point(aes(color = as.factor(chr)), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = rep(c("black", "grey"), length.out = length(unique(pi_des_out_plots_data$chr)))) +
  geom_hline(yintercept = quantile(-log10(pi_des_out_plots$pval_pi), 0.95), color = "red", linetype = "dashed") +
  scale_x_continuous(label = pi_axis_set$chr, breaks = pi_axis_set$center) +
  labs(x = "Chromosome Number", y = "-log10(p-value)", title = "DES Pi ratio (SC/SI)") +
  theme_minimal() +
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+ 
  annotate("text", x = 921973854457, y = 13, label = "MIK2/LRRK", color="blue",)+
  geom_segment(aes(x=921973854457,yend = 12), color="blue")

#fst
FST_des_out_plots <- fst_des_outcrossing
FST_des_out_plots$CHROM <- gsub("Scaffold_","",FST_des_out_plots$CHROM)
FST_des_out_plots$chr <- factor(FST_des_out_plots$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
FST_des_out_plots_data <- FST_des_out_plots %>%
  arrange(chr, BIN_MID) %>%
  group_by(chr) %>%
  mutate(pos_max_this_chr = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(bp_add = lag(cumsum(pos_max_this_chr), default = 0)) %>%
  group_by(chr) %>%
  mutate(bp_cum = BIN_MID + bp_add - min(BIN_MID)) %>% 
  ungroup()

fst_axis_set <- FST_des_out_plots_data %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

des_plot_fst <- ggplot(FST_des_out_plots_data, aes(x = bp_cum, y = -log10(pval_fst))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = rep(c("black", "grey"), length.out = length(unique(FST_des_out_plots_data$chr)))) +
  geom_hline(yintercept = quantile(-log10(FST_des_out_plots$pval_fst), 0.95), color = "red", linetype = "dashed") +
  scale_x_continuous(label = fst_axis_set$chr, breaks = fst_axis_set$center) +
  labs(x = "Chromosome Number",y = "-log10(p-value)",title = "FST") +
  theme_minimal() +
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+ 
  annotate("text", x = 973649404706, y = 18, label = "MIK2/LRRK", color="blue",)+
  geom_segment(aes(x=973649404706, yend = 17), color="blue")

#arrange DES plots
des_genome <- grid.arrange(des_plot_pi, des_plot_fst, ncol = 1)

#processing data and plotting Fst & Pi for WAT
#pi ratio
pi_wat_out_plots <- pi_data_outcrossing_wat
pi_wat_out_plots$CHROM <- gsub("Scaffold_","",pi_wat_out_plots$CHROM)
pi_wat_out_plots$chr <- factor(pi_wat_out_plots$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
pi_wat_out_plots_data <- pi_wat_out_plots %>%
  arrange(chr, BIN_MID) %>%
  group_by(chr) %>%
  mutate(pos_max_this_chr = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(bp_add = lag(cumsum(pos_max_this_chr), default = 0)) %>%
  group_by(chr) %>%
  mutate(bp_cum = BIN_MID + bp_add - min(BIN_MID)) %>% 
  ungroup()

pi_axis_set <- pi_wat_out_plots_data %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

wat_plot_pi <- ggplot(pi_wat_out_plots_data, aes(x = bp_cum, y = -log10(pval_pi))) +
  geom_point(aes(color = as.factor(chr)), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = rep(c("black", "grey"), length.out = length(unique(pi_wat_out_plots_data$chr)))) +
  geom_hline(yintercept = quantile(-log10(pi_wat_out_plots$pval_pi), 0.95), color = "red", linetype = "dashed") +
  scale_x_continuous(label = pi_axis_set$chr, breaks = pi_axis_set$center) +
  labs(x = "Chromosome Number", y = "-log10(p-value)", title = "WAT Pi ratio (SC/SI)") +
  theme_minimal() +
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank())+ 
  annotate("text", x = 227136501026, y = 11, label = "MIK2/LRRK", color="blue",)+
  geom_segment(aes(x=227136501026, yend=10.5), color="blue")

#fst
FST_wat_out_plots <- fst_wat_outcrossing
FST_wat_out_plots$CHROM <- gsub("Scaffold_","",FST_wat_out_plots$CHROM)
FST_wat_out_plots$chr <- factor(FST_wat_out_plots$CHROM, levels = c("1", "2", "3", "4", "5", "6", "7", "8"))
FST_wat_out_plots_data <- FST_wat_out_plots %>%
  arrange(chr, BIN_MID) %>%
  group_by(chr) %>%
  mutate(pos_max_this_chr = max(BIN_MID)) %>%
  ungroup() %>%
  mutate(bp_add = lag(cumsum(pos_max_this_chr), default = 0)) %>%
  group_by(chr) %>%
  mutate(bp_cum = BIN_MID + bp_add - min(BIN_MID)) %>% 
  ungroup()

fst_axis_set <- FST_wat_out_plots_data %>%
  group_by(chr) %>%
  summarize(center = mean(bp_cum))

wat_plot_fst <- ggplot(FST_wat_out_plots_data, aes(x = bp_cum, y = -log10(pval_fst))) + 
  geom_point(aes(color = as.factor(chr)), alpha = 0.7, size = 1.5) +
  scale_color_manual(values = rep(c("black", "grey"), length.out = length(unique(FST_wat_out_plots_data$chr)))) +
  geom_hline(yintercept = quantile(-log10(FST_wat_out_plots$pval_fst), 0.95), color = "red", linetype = "dashed") +
  scale_x_continuous(label = fst_axis_set$chr, breaks = fst_axis_set$center) +
  labs(x = "Chromosome Number",y = "-log10(p-value)",title = "FST") +
  theme_minimal() +
  theme(legend.position = "none",panel.grid.major.x = element_blank(),panel.grid.minor.x = element_blank()) + 
  annotate("text", x = 279313501262, y = 9, label = "MIK2/LRRK", color="blue",)+
  geom_segment(aes(x=279313501262, yend = 8.5), color="blue")

#arrange WAT plots
wat_genome <- grid.arrange(wat_plot_pi, wat_plot_fst, ncol = 1)

#final plot arrange
all <- grid.arrange(des_genome,wat_genome, ncol=1)
