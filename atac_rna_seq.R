# get sex/asex specific genes
sex_specific <- read.table("/home/wanghm/whm/1FG/xulab_current/sex_specific_newID", stringsAsFactors = F)
asex_specific <- read.table("/home/wanghm/whm/1FG/xulab_current/asex_specific_newID", stringsAsFactors = F)

tpm <- function(counts,lengths) {
  rate <- counts/lengths
  rate/sum(rate) * 1e6
}
# read exp matrix
a <- read.table("/home/wanghm/whm/ATAC/heatmap.data/atac_promoter_count_sex", header = T)
a <- a[,c(1,6,7)]
b <- read.table("/home/wanghm/whm/ATAC/heatmap.data/atac_promoter_count_hypha", header = T)
b <- b[,c(1,6,7)]

hypha.lib.size = 30572943
sex.lib.size = 84460558
library(edgeR) ## calculate RPKM of each promoter regions and corresponding gene expression
pmt_rpkm1 <- rpkm(y = a[,3], gene.length = a$Length, lib.size = sex.lib.size)
atacExp.df1 <- data.frame(geneID=a$Geneid, ATAC_pmt_exp.sex=pmt_rpkm1)

pmt_rpkm2 <- rpkm(y = b[,3], gene.length = a$Length, lib.size = hypha.lib.size)
atacExp.df2 <- data.frame(geneID=b$Geneid, ATAC_pmt_exp.hypha=pmt_rpkm2)
rna_seq <- read.table("/home/wanghm/wanghm_R/1.res/expr_newID.txt", header = T)

library(magrittr)
# dd <- merge(atacExp.df1, rna_seq, by.x="geneID", by.y = "new") %>% merge(., atacExp.df2, by.x="geneID", by.y="geneID")
# new.exp <- read.table("/home/wanghm/wanghm_R/1.res/expr_newID.txt", header = T, stringsAsFactors = F)
# top_500 <- new.exp[order(new.exp$sex_8d,decreasing = T),][1:500,]
# as <- dd[dd$geneID %in%top_500$new,]
### squeeze out stage-specific genes from dd
dd <- dd[dd$geneID %in% sex_specific$V1,]   # nrow(dd)
dd[order(dd$sex_8d, decreasing = T),]
heat.mat <- as.data.frame(dd[,c(3,2)])
scale(dd[,c(4,5,3,2)])
heat.mat <- heat.mat[order(heat.mat$sex_8d,decreasing = T),]
range(heat.mat$sex_8d)
bk <- c(seq(-1.6,-0.1,by=0.01), seq(0, 1.6,by = 0.01))
pheatmap::pheatmap(heat.mat,color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                  colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                  breaks = bk, cluster_rows = T, cluster_cols = F,scale = "")

### plot heatmap of sex peak score and corresponding gene expression
fc_score.sex <- read.table("/home/wanghm/whm/novo_atac/peak_calling/novo_sex_peaks.gappedPeak.bed")$V4
sex.gr$score <- fc_score.sex
sex.geneID <- sex.gr[sex.gr$annotation == "Promoter"]$geneId
sex.peakscore <- sex.gr[sex.gr$annotation == "Promoter"]$score
score.df <- data.frame(sex.geneID, sex.peakscore)
score_exp.df <- merge(score.df, rna_seq, by.x="sex.geneID", by.y="new")

bk <- c(seq(-1.5,-0.1,by=0.01), seq(1.5,by = 0.01))
pheatmap::pheatmap(scale(score_exp.df[,c(2,3)]),color = c(colorRampPalette(colors = c("blue","white"))(length(bk)/2),
                                                          colorRampPalette(colors = c("white","red"))(length(bk)/2)),
                                                          breaks = bk, cluster_rows = F, cluster_cols = F)

