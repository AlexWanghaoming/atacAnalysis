# analysis relationship between histone modification and gene expression 

library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")

library(ChIPseeker)
pmt <- getPromoters(txdb, downstream = 1000, upstream = 1000)
ge <- genes(txdb)

###  relationship with genes exp
exp <- read.table("/home/wanghm/whm/1FG/FG_expr.v2018")
sex_8d.exp <- exp[,c(33,34)]  # get sex-8d gene expression
sex_8d.exp$sex_8d <- rowMeans(sex_8d.exp)  # add new columns----average of R1 and R2

gene_mapping1 <- read.table("/home/wanghm/whm/chip_seq/new_old_gene_mapping.txt", col.names = c("new", "old1"))
gene_mapping2 <- read.table("/home/wanghm/whm/1FG/rres2ramph1_gene_mapping.txt", col.names = c("rres", "ramph1"))
gene.id <- merge(gene_mapping1,gene_mapping2, by.x="old1", by.y="ramph1")
sex_8d_allgeneID.exp <- merge(sex_8d.exp, gene.id, by.x="row.names", by.y="rres")

histone_peakBindin_comparison <- function(type) {      ## type is one of: h3k4_me, h3k9_ac, h3k9_me, h3k27_ac
  ## check chip-binding genes
  if (type=="h3k4_me") {
    peak <- peak2
  }else if(type == "h3k9_ac"){
    peak <- peak3
  }else if(type == "h3k9_me"){
    peak <- peak4
  }else if(type == "h3k27_ac"){
    peak <- peak5
  }
  library(tidyverse)
  chip_binding.genes <- subsetByOverlaps(ge, peak) %>%
    mcols() %>%
    rownames()
  chip_binding.exp <- sex_8d_allgeneID.exp[sex_8d_allgeneID.exp$new %in% chip_binding.genes,][,4]
  other.exp <- sex_8d_allgeneID.exp[!sex_8d_allgeneID.exp$new %in% chip_binding.genes,][,4]
  level <- c(rep("peaks",length(chip_binding.exp)),rep("no_peaks",length(other.exp)))
  cpm <- c(chip_binding.exp, other.exp)
  gg.df <- data.frame(level, value)
  random_palette <- c("Set1", "Set2", "Accent", "Dark2", "Set3", "Pastel2", "Paired", "Pastel1")
  col <- sample(random_palette,1)
  library(ggstatsplot)
  p <- ggbetweenstats(data = gg.df, x = level, y = value,notch = T, xlab = "peak binding or not", ylab = "cpm",
                      title = type, palette = col) + ylim(0,200)
  return(p)
}

p1 <- histone_peakBindin_comparison("h3k4_me")
p2 <- histone_peakBindin_comparison("h3k9_ac")
p3 <- histone_peakBindin_comparison("h3k9_me")
p4 <- histone_peakBindin_comparison("h3k27_ac")

library(cowplot)
cowplot::plot_grid(p1, p2, p3, p4, labels = c("a","b","c","d"))











