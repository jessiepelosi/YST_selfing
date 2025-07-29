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

###DES
# set path
des_my_bins <- "des_1na.ld_decay_bins"
# read in data
des_ld_bins <- read_tsv(des_my_bins, show_col_types = FALSE)
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
# set path
wat_my_bins <- "wat_1na.ld_decay_bins"
# read in data
wat_ld_bins <- read_tsv(wat_my_bins, show_col_types = FALSE)
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
#set path
col_my_bins <- "col_Scaffold_1.ld_decay_bins"
# read in data
col_ld_bins <- read_tsv(col_my_bins, show_col_types = FALSE)
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
#set path
dia_my_bins <- "dia_1na.ld_decay_bins"
# read in data
dia_ld_bins <- read_tsv(dia_my_bins, show_col_types = FALSE)
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
#set path
gil_my_bins <- "gil_1na.ld_decay_bins"
# read in data
gil_ld_bins <- read_tsv(gil_my_bins, show_col_types = FALSE)
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
#set path
gol_my_bins <- "gol_1na.ld_decay_bins"
# read in data
gol_ld_bins <- read_tsv(gol_my_bins, show_col_types = FALSE)
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
#set path
leb_my_bins <- "leb_1na.ld_decay_bins"
# read in data
leb_ld_bins <- read_tsv(leb_my_bins, show_col_types = FALSE)
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
#set path
Nee_my_bins <- "nee_1na.ld_decay_bins"
# read in data
Nee_ld_bins <- read_tsv(Nee_my_bins, show_col_types = FALSE)
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
#set path
orc_my_bins <- "orc_1na.ld_decay_bins"
# read in data
orc_ld_bins <- read_tsv(orc_my_bins, show_col_types = FALSE)
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
#set path
oro_my_bins <- "oro_1na.ld_decay_bins"
# read in data
oro_ld_bins <- read_tsv(oro_my_bins, show_col_types = FALSE)
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
#set path
rb_my_bins <- "rb_1na.ld_decay_bins"
# read in data
rb_ld_bins <- read_tsv(rb_my_bins, show_col_types = FALSE)
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
res_my_bins <- "res_1naa.ld_decay_bins"
# read in data
res_ld_bins <- read_tsv(res_my_bins, show_col_types = FALSE)
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
#set path
sie_my_bins <- "sie_1na.ld_decay_bins"
# read in data
sie_ld_bins <- read_tsv(sie_my_bins, show_col_types = FALSE)
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
#set path
tri_my_bins <- "tri_1na.ld_decay_bins"
# read in data
tri_ld_bins <- read_tsv(tri_my_bins, show_col_types = FALSE)
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
uki_my_bins <- "uki_1na.ld_decay_bins"
# read in data
uki_ld_bins <- read_tsv(uki_my_bins, show_col_types = FALSE)
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
#set path
vet_my_bins <- "vet_1nad.ld_decay_bins"
# read in data
vet_ld_bins <- read_tsv(vet_my_bins, show_col_types = FALSE)
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
#set path
yre_my_bins <- "yre_1nad.ld_decay_bins"
# read in data
yre_ld_bins <- read_tsv(yre_my_bins, show_col_types = FALSE)
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

mean_ld_dia <- dia_ld_bins %>%
  group_by(chr) %>% 
  summarize(mean = mean(avg_R2), sd = sd(avg_R2))
