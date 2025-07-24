#############################
#GoGetter Heatmap for California Invasion
#Taylor Curry, Abby Pearse, Jessie Pelosi
#Dlugosch lab
#July 24, 2025
############################

library(pheatmap)
library(viridis)
library(grid)
library(dplyr)

#Data
des_gg <- read.delim("des_gene_id_window_seqs.fa.blast.besthits.tsv.freqcounts-gene.tsv")
wat_gg <- read.delim("wat_gene_id_window_seqs.fa.blast.besthits.tsv.freqcounts-gene.tsv")
cds_gg <- read.delim("clean_cds_from_genomic.fna.blast.besthits.tsv.freqcounts-gene.tsv")

des_wat_gg_matrix <- full_join(des_gg,wat_gg, by = "GOSlimTerm")
des_wat_cds_gg_matrix <- full_join(des_wat_gg_matrix,cds_gg, by = "GOSlimTerm")
colnames(des_wat_cds_gg_matrix) <- c("GOSlimTerm","DES","WAT","GENOME")

rownames(des_wat_cds_gg_matrix) <- des_wat_cds_gg_matrix$GOSlimTerm
des_wat_cds_gg_matrix <- des_wat_cds_gg_matrix %>% replace(is.na(.), 0)
des_wat_cds_gg_matrix <- des_wat_cds_gg_matrix %>% select(-GOSlimTerm)

#pheatmap
heatmap_data <- des_wat_cds_gg_matrix
ordered <- heatmap_data[order(des_wat_cds_gg_matrix$GENOME, decreasing = T),]
heat_matrix <- as.matrix(ordered)

draw_colnames_75 <- function (coln, gaps, ...) {
  coord <- pheatmap:::find_coordinates(length(coln), gaps)
  x     <- coord$coord - 0.5 * coord$size
  res   <- grid::textGrob(
    coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"),
    vjust = 0.75, hjust = 1, rot = 75, gp = grid::gpar(...)
  )
  return(res)
}
assignInNamespace(
  x = "draw_colnames",
  value = "draw_colnames_75",
  ns = asNamespace("pheatmap")
)

quantile_breaks <- function(xs, n = 50) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}
mat_breaks <- quantile_breaks(heat_matrix, n = 20)

#plot with pheatmap
pheatmap(heat_matrix,
         breaks = mat_breaks,
         cellwidth = 25,
         cellheight = 9,
         fontsize = 8,
         cluster_rows = F,
         cluster_cols = F, 
         border_color = "NA",
         col = viridis(20),
         angle_col = "45",
         display_numbers = F,
         main ="GOSlim Heatmap (Frequency)"
         )

