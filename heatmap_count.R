#############################
#GoGetter Gene count heatmap 
#Centaurea solstitialis, California Invasion
#Taylor Curry, Abby Pearse, Jessie Pelosi
#Dlugosch lab
#July 24, 2025
############################

library(pheatmap)
library(dplyr)
library(viridis)
library(grid)

#Read in data and construct matrix
des_genes_gg <- read.delim("des_gene_id_window_seqs.fa.blast.besthits.tsv.rawcounts-gene.tsv")
wat_genes_gg <- read.delim("wat_gene_id_window_seqs.fa.blast.besthits.tsv.rawcounts-gene.tsv")
full_genome_genes_gg <- read.delim("clean_cds_from_genomic.fna.blast.besthits.tsv.rawcounts-gene.tsv")

des_wat_genes_gg <- full_join(des_genes_gg,wat_genes_gg, by = "GOSlimTerm")
all_genes_gg <- full_join(des_wat_genes_gg,full_genome_genes_gg, by = "GOSlimTerm")

colnames(all_genes_gg) <- c("GOSlimTerm","DES","WAT","Genome")
rownames(all_genes_gg) <- all_genes_gg$GOSlimTerm
all_genes_gg <- all_genes_gg %>% replace(is.na(.), 0)
all_genes_gg <- all_genes_gg %>% select(-GOSlimTerm)

#Chi Sq Test
all_genes_gg_x2<- chisq.test(all_genes_gg)
all_genes_gg_residuals<- as.data.frame(all_genes_gg_x2$residuals)

#heatmap
ordered <- all_genes_gg_residuals[order(all_genes_gg_residuals$Genome, decreasing = T),]
all_genes_gg_residuals <- as.matrix(ordered)

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

quantile_breaks <- function(xs, n = 96) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}

mat_breaks <- quantile_breaks(all_genes_gg_residuals, n = 20)

#Plot heatmap
pheatmap(all_genes_gg_residuals,
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
         main ="GOSlim Heatmap"
         )

