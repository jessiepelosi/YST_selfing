############
#Linkage Disequilibrium Decay
#Centaurea sol. California Invasion 
#Taylor Curry, Abby Peasrse, Jessie Pelosi
#June 6, 2025
##########

library(ggplot2)
library(ggpmisc)
library(dplyr)
library(readr)
library(pheatmap)
library(viridis)

###DES
#set path and read in data
des_ld_bins <- read_tsv("des_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(des_ld_bins, aes(distance, avg_R2)) +
  geom_line() +
  xlab("Distance (bp)") +
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000)+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("DES")

###WAT
# set path and read in data
wat_ld_bins <- read_tsv("wat_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(wat_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") +
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth()+
  theme_bw()+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("WAT")

####COL
#set path and read in data
col_ld_bins <- read_tsv("col_Scaffold_1.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(col_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("COL")

###DIA
#set path and read in data
dia_ld_bins <- read_tsv("dia_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(dia_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("DIA")

###GIL
#set path and read in data
gil_ld_bins <- read_tsv("gil_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(gil_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("GIL")

###GOL
#set path and read in data
gol_ld_bins <- read_tsv("gol_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(gol_ld_bins, aes(distance, avg_R2)) +
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr,)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("GOL")

###LEB
#set path and read in data
leb_ld_bins <- read_tsv("leb_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(leb_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("LEB")

###NEE
#set path and read in data
Nee_ld_bins <- read_tsv("nee_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(Nee_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") +
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000)+
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("NEE")

###ORC
#set path and read in data
orc_ld_bins <- read_tsv("orc_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(orc_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("ORC")

###ORO
#set path and read in data
oro_ld_bins <- read_tsv("oro_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(oro_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("ORO")

###RB
#set path and read in data
rb_ld_bins <- read_tsv("rb_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(rb_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("RB")

###RES
#set path and read in data
res_ld_bins <- read_tsv("res_1naa.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(res_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("RES")

###SIE
#set path and read in data
sie_ld_bins <- read_tsv("sie_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(sie_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") +
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("SIE")

###TRI
#set path and read in data
tri_ld_bins <- read_tsv("tri_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(tri_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("TRI")

###UKI
#set path and read in data
uki_ld_bins <- read_tsv("uki_1na.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(uki_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("UKI")

###VET
#set path and read in data
vet_ld_bins <- read_tsv("vet_1nad.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(vet_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("VET")

###YRE
#set path and read in data
yre_ld_bins <- read_tsv("yre_1nad.ld_decay_bins", show_col_types = FALSE)
# plot LD decay
ggplot(yre_ld_bins, aes(distance, avg_R2)) + 
  geom_line() +
  xlab("Distance (bp)") + 
  ylab(expression(italic(r)^2)) +
  facet_wrap(~chr)+
  xlim(0,50000) + 
  geom_smooth(method = "lm")+
  theme_bw()+
  stat_poly_eq(use_label(c("eq")))+
  ggtitle("YRE")

#calculate mean and std dev
mean_ld_col <- col_ld_bins %>%
 group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_des <- des_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_wat <- wat_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_dia <- dia_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_gil <- gil_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_gol <- gol_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_leb <- leb_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_nee <- Nee_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_orc <- orc_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_oro <- oro_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_rb <- rb_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_res <- res_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_sie <- sie_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_tri <- tri_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_uki <- uki_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_vet <- vet_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

mean_ld_yre <- yre_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))

#create data frame for pheatmap to plot
LD_means <- as.data.frame(cbind(mean_ld_des$mean,
               mean_ld_wat$mean,
               mean_ld_col$mean,
               mean_ld_dia$mean,
               mean_ld_gil$mean,
               mean_ld_gol$mean,
               mean_ld_leb$mean,
               mean_ld_nee$mean,
               mean_ld_orc$mean,
               mean_ld_oro$mean,
               mean_ld_rb$mean,
               mean_ld_res$mean,
               mean_ld_sie$mean,
               mean_ld_tri$mean,
               mean_ld_uki$mean,
               mean_ld_vet$mean,
               mean_ld_yre$mean))

colnames(LD_means) <- c("DES",
                    "WAT",
                    "COL",
                    "DIA",
                    "GIL",
                    "GOL",
                    "LEB",
                    "NEE",
                    "ORC",
                    "ORO",
                    "RB",
                    "RES",
                    "SIE",
                    "TRI",
                    "UKI",
                    "VET",
                    "YRE")

rownames(LD_means) <- c("Chromosome 1",
                     "Chromosome 2",
                     "Chromosome 3",
                     "Chromosome 4",
                     "Chromosome 5",
                     "Chromosome 6",
                     "Chromosome 7",
                     "Chromosome 8")
#Plot heatmap
pheatmap(LD_means,
         cellwidth = 25,
         cellheight =25,
         border_color = NA,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F, 
         color = viridis(20),
         display_numbers = F,
         number_color = "WHITE",
         fontsize_number = 9,
         angle_col = 0
)
