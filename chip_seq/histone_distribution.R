# anaylysis and plot profiles of histone modification distribution in genome

library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")

library(ChIPseeker)
pmt <- getPromoters(txdb, downstream = 1000, upstream = 1000)
ge <- genes(txdb)

peak1 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K14_Ac_TBI.peaks");colnames(mcols(peak1))[1:2] <- c("p", "score")
peak2 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K4_Me_sex.peaks");colnames(mcols(peak2))[1:2] <- c("p", "score")
peak3 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K9_Ac_sex.peaks");colnames(mcols(peak3))[1:2] <- c("p", "score")
peak4 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K9_Me_sex.peaks");colnames(mcols(peak4))[1:2] <- c("p", "score")
peak5 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K27_Ac_sex.peaks");colnames(mcols(peak5))[1:2] <- c("p", "score")
peak6 <- readPeakFile("/home/wanghm/whm/chip_seq/hyp24h-H3K9ME3.peaks");colnames(mcols(peak6))[1:2] <- c("p", "score")
peak7 <- readPeakFile("/home/wanghm/whm/chip_seq/h3.peaks");colnames(mcols(peak7))[1:2] <- c("p", "score")
peak8 <- readPeakFile("/home/wanghm/whm/chip_seq/h3K27_Me_sex.peaks");colnames(mcols(peak8))[1:2] <- c("p", "score")
## Tss distribution
tag.mat1 <- getTagMatrix(peak1, windows = pmt)
tag.mat2 <- getTagMatrix(peak2, windows = pmt)
tag.mat3 <- getTagMatrix(peak3, windows = pmt)
tag.mat4 <- getTagMatrix(peak4, windows = pmt)
tag.mat5 <- getTagMatrix(peak5, windows = pmt)
tag.mat6 <- getTagMatrix(peak6, windows = pmt)
tag.mat7 <- getTagMatrix(peak7, windows = pmt)
tag.mat8 <- getTagMatrix(peak8, windows = pmt)

# make tag list obj
tag.list1 <- list(tag.mat2, tag.mat3, tag.mat4, tag.mat5, tag.mat8,tag.mat6, tag.mat7)
names(tag.list1) <- c("sex-h3k4_ME", "sex-h3k9_AC", "sex-h3k9_ME","sex-h3k27_AC", "sex-h3k27_ME", "hypha24h-h3k9_ME", "h3")
plotAvgProf(tag.list1, xlim = c(-1000, 1000))

# plot ChIP-seq signals genome distribution
library(ggplot2)
library(tidyverse)
library(cowplot)

## assign environment of covplot - > environment of mod_cov.plot 
environment(mod_cov.plot) <- environment(covplot) 

## add ATAC-seq peak info
atac_sex_peak <- readPeakFile("/home/wanghm/whm/ATAC/rmPCR_duplicate/sex-7d_all/sex-7d_all_peaks.narrowPeak")
colnames(mcols(atac_sex_peak))[2] <- c("score")
atac_hypha_peak <- readPeakFile("/home/wanghm/whm/ATAC/rmPCR_duplicate/hypha_all/hypha24_all_peaks.narrowPeak")
colnames(mcols(atac_hypha_peak))[2] <- c("score")

plot_histone_profile <- function(chr, ac_col, me_col, h3_col, atac_col){
  # histone ac
  p1 <- mod_cov.plot(peak1, weightCol = "score", chrs = chr, co = ac_col, title = "h3k14_ac_TBI")
  p2 <- mod_cov.plot(peak3, weightCol = "score", chrs = chr, co = ac_col, title = "h3k9_ac_sex")
  p3 <- mod_cov.plot(peak5, weightCol = "score", chrs = chr, co = ac_col, title = "h3k27_ac_sex")
  # histone me
  p4 <- mod_cov.plot(peak2, weightCol = "score", chrs = chr, co = me_col, title = "h3k4_Me_sex")
  p5 <- mod_cov.plot(peak4, weightCol = "score", chrs = chr, co = me_col, title = "h3k9_Me_sex")
  p6 <- mod_cov.plot(peak6, weightCol = "score", chrs = chr, co = me_col, title = "h3k9_Me_hyp24h")
  p7 <- mod_cov.plot(peak8, weightCol = "score", chrs = chr, co = me_col, title = "h3k27_Me_sex")
  
  p8 <- mod_cov.plot(peak7, weightCol = "score", chrs = chr, co = h3_col, title = "h3")
  p9 <- mod_cov.plot(atac_sex_peak, weightCol = "score", chrs = chr, co = atac_col, title = "sex_ATAC")
  p10 <- mod_cov.plot(atac_hypha_peak, weightCol = "score", chrs = chr, co = atac_col, title = "hypha_ATAC") + ylim(0, 1000)
  
  
  plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10, ncol = 1, align = "v")
  # save_plot(filename = "" , p, ncol = 1)
}
plot_histone_profile(chr = "chr1", ac_col ="#106f8e", me_col = "#ea725c", h3_col = "#e7d163", atac_col = "#3a7d44")
