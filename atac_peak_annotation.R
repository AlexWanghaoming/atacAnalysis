# a = ls()
# rm(list = a[!(a%in%c("dbObj", "dbObj_count", "dbObj_contrast", "dbObj_analysis", "comp1.deseq2", "comp1.edgeR"))])

library(GenomicFeatures)
atac.txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/xulab_current/PH1-fusion-version2_represent_transcript_rm_Mt.gff")

## plot peak profile
library(ChIPseeker)
pmt <- getPromoters(atac.txdb, downstream = 500, upstream = 1000)
#### build a terminator GRanges
terminator <- resize(genes(atac.txdb),width = 1000, fix = "end")
end(terminator[strand(terminator)=="+"]) <- end(terminator[strand(terminator)=="+"]) + 501
start(terminator[strand(terminator)=="-"]) <- start(terminator[strand(terminator)=="-"]) - 501
terminator <- terminator[start(terminator)>0]

## peak annotation Using ChIPpeakAnno package
library(ChIPpeakAnno)
library(ChIPseeker)
sex_peak <- toGRanges("/home/wanghm/whm/novo_atac/peak_calling/novo_sex_peaks.gappedPeak.bed")
hypha_peak <- toGRanges("/home/wanghm/whm/ATAC/peak_calling_hmmr/open_bed/hypha_all_peaks.bed")
sex_sp_peak <- sex_peak[!(sex_peak %in% subsetByOverlaps(sex_peak, hypha_peak, minoverlap = 100))]
hypha_sp_peak <-hypha_peak[!(hypha_peak %in% subsetByOverlaps(hypha_peak, sex_peak, minoverlap = 100))]
# sex_peak_summit <- toGRanges("/home/wanghm/whm/novo_atac/peak_calling/novo_sex_summits.bed")
# hypha_peak_summit <- toGRanges("/home/wanghm/whm/ATAC/peak_calling_hmmr/hypha_all_summits.bed")
# sex_sp_summit <- sex_peak_summit[which(sex_peak %in% sex_sp_peak)]
# hypha_sp_summit <- hypha_peak_summit[which(hypha_peak %in% hypha_sp_peak)]

# get peak distribution matrix using score as weight
sex_peak.mat <- getTagMatrix(peak = sex_peak, windows = pmt)
hypha_peak.mat <- getTagMatrix(peak = hypha_peak, windows = pmt)
sex_sp.mat <- getTagMatrix(peak = sex_sp_peak, windows = pmt)
hypha_sp.mat <- getTagMatrix(peak=hypha_sp_peak, windows = pmt)
# for TES
sex_peak.tts.mat <- getTagMatrix(peak = sex_peak, windows = terminator)
hypha_peak.tts.mat <- getTagMatrix(peak = hypha_peak, windows = terminator)
sex_sp.tts.mat <- getTagMatrix(peak = sex_sp_peak, windows = terminator)
hypha_sp.tts.mat <- getTagMatrix(peak=hypha_sp_peak, windows = terminator)
# plot peak distribution profile around TSS
ll = list(hypha_peak.mat, sex_peak.mat)
names(ll) <- c("hypha-24", "sex-7d")
p1 <- plotAvgProf(ll, xlim = c(-500, 1000), lwd=2) + ylim(0.00050, 0.00080)
ll2 = list(hypha_peak.tts.mat, sex_peak.tts.mat)
names(ll2) <- c("hypha-24", "sex-7d")
p2 <- plotAvgProf(ll2, xlim = c(-1000, 500), lwd=2) +  ylim(0.00050, 0.00080)

cowplot::plot_grid(p1,p2)
ChIPseeker:::plotAvgProf.internal
# peak distribution on exon, intron, UTR, region
library(ChIPseeker)
library(magrittr)
# ?annotatePeak
sex.gr <- annotatePeak(sex_peak, tssRegion = c(-500,300), TxDb = atac.txdb, 
                            genomicAnnotationPriority=c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), 
                            addFlankGeneInfo = T) %>% as.GRanges()
sex_sp.gr <- annotatePeak(sex_sp_peak, tssRegion = c(-500,300), TxDb = atac.txdb, 
                       genomicAnnotationPriority=c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), 
                       addFlankGeneInfo = T) %>% as.GRanges()
hypha.gr <- annotatePeak(hypha_peak, tssRegion = c(-500,300), TxDb = atac.txdb, 
                         genomicAnnotationPriority=c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), 
                         addFlankGeneInfo = T) %>% as.GRanges()
hypha_sp.gr <- annotatePeak(hypha_sp_peak, tssRegion = c(-500,300), TxDb = atac.txdb, 
                         genomicAnnotationPriority=c("Promoter", "5UTR", "3UTR", "Exon", "Intron", "Intergenic"), 
                         addFlankGeneInfo = T) %>% as.GRanges()
# output gff
ranges(sex.gr) <- ranges(sex_peak)
ranges(hypha.gr) <- ranges(hypha_peak)
ranges(sex_sp.gr) <- ranges(sex_sp_peak)
ranges(hypha_sp.gr) <- ranges(hypha_sp_peak)
rtracklayer::export.gff3(sex.gr, con = "/home/wanghm/whm/ATAC/peak_calling_hmmr/cal_rpkm/sex_peakAnno.gff")
rtracklayer::export.gff3(hypha.gr, con = "/home/wanghm/whm/ATAC/peak_calling_hmmr/cal_rpkm/hypha_peakAnno.gff")
rtracklayer::export.gff3(sex_sp.gr, con = "/home/wanghm/whm/ATAC/peak_calling_hmmr/cal_rpkm/sex_sp_peakAnno.gff")
# rtracklayer::export.gff3(hypha_sp.gr, con = "/home/wanghm/whm/ATAC/peak_calling/cal_rpkm/hypha_sp_peakAnno.gff")


# extract annotation info
sex_Anno <- sex.gr$annotation
sex_sp_Anno <- sex_sp.gr$annotation
hypha_Anno <- hypha.gr$annotation
hypha_sp_Anno <- hypha_sp.gr$annotation
library(dplyr)

# vars selection
a1 <- length(starts_with("Exon", vars = sex_Anno))
a2 <- length(starts_with("Intron", vars = sex_Anno))
a3 <- length(starts_with("Promoter", vars = sex_Anno))
a4 <- length(starts_with("3' UTR", vars = sex_Anno))
a5 <- length(starts_with("5' UTR", vars = sex_Anno))
a6 <- length(starts_with("Downstream", vars = sex_Anno))
a7 <- length(starts_with("Distal Intergenic",vars = sex_Anno))
# vars selection
aa1 <- length(starts_with("Exon", vars = sex_sp_Anno))
aa2 <- length(starts_with("Intron", vars = sex_sp_Anno))
aa3 <- length(starts_with("Promoter", vars = sex_sp_Anno))
aa4 <- length(starts_with("3' UTR", vars = sex_sp_Anno))
aa5 <- length(starts_with("5' UTR", vars = sex_sp_Anno))
aa6 <- length(starts_with("Downstream", vars = sex_sp_Anno))
aa7 <- length(starts_with("Distal Intergenic",vars = sex_sp_Anno))
# vars selection
b1 <- length(starts_with("Exon", vars = hypha_Anno))
b2 <- length(starts_with("Intron", vars = hypha_Anno))
b3 <- length(starts_with("Promoter", vars = hypha_Anno))
b4 <- length(starts_with("3' UTR", vars = hypha_Anno))
b5 <- length(starts_with("5' UTR", vars = hypha_Anno))
b6 <- length(starts_with("Downstream", vars = hypha_Anno))
b7 <- length(starts_with("Distal Intergenic",vars = hypha_Anno))
# vars selection
bb1 <- length(starts_with("Exon", vars = hypha_sp_Anno))
bb2 <- length(starts_with("Intron", vars = hypha_sp_Anno))
bb3 <- length(starts_with("Promoter", vars = hypha_sp_Anno))
bb4 <- length(starts_with("3' UTR", vars = hypha_sp_Anno))
bb5 <- length(starts_with("5' UTR", vars = hypha_sp_Anno))
bb6 <- length(starts_with("Downstream", vars = hypha_sp_Anno))
bb7 <- length(starts_with("Distal Intergenic",vars = hypha_sp_Anno))

## data frame
sex.df <- data.frame(group = c("Promoter", "Exon", "Intron", "5'UTR", "3'UTR", "Intergenic"), 
           counts=c(a3,a1,a2,a5,a4,a6+a7), 
           percentage=round(sapply(c(a3,a1,a2,a5,a4,a6+a7), function(x) x/sum(a3,a1,a2,a5,a4,a6+a7)),3))
hypha.df <- data.frame(group = c("Promoter", "Exon", "Intron", "5'UTR", "3'UTR", "Intergenic"), 
           counts=c(b3,b1,b2,b5,b4,b6+b7),
           percentage=round(sapply(c(b3,b1,b2,b5,b4,b6+b7), function(x) x/sum(b3,b1,b2,b5,b4,b6+b7)),3))
sex.df.sp <- data.frame(group = c("Promoter", "Exon", "Intron", "5'UTR", "3'UTR", "Intergenic"), 
                     counts=c(aa3,aa1,aa2,aa5,aa4,aa6+aa7), 
                     percentage=round(sapply(c(aa3,aa1,aa2,b5,b4,b6+b7), function(x) x/sum(aa3,aa1,aa2,aa5,aa4,aa6+aa7)),3))
hypha.df.sp <- data.frame(group = c("Promoter", "Exon", "Intron", "5'UTR", "3'UTR", "Intergenic"), 
                       counts=c(bb3,bb1,bb2,bb5,bb4,bb6+bb7),
                       percentage=round(sapply(c(bb3,bb1,bb2,bb5,bb4,bb6+bb7), function(x) x/sum(bb3,bb1,bb2,bb5,bb4,bb6+bb7)),3))
library(tidyverse)
df <- rbind(sex.df,sex.df.sp,hypha.df,hypha.df.sp) %>% add_column(., type=c(rep("sex", 12), rep("hypha",12))) %>%add_column(., class=c(rep("sex",6),
                                                                                                                                        rep("sex_sp",6),
                                                                                                                                        rep("hypha",6),
                                                                                                                                        rep("hypha_sp",6)))
## plot stack barplot of peak annotation sites
# p <- ggbarplot(df,x = "type",y="counts", fill = "class", color="white",
#                palette = RColorBrewer::brewer.pal(6, "Set3"), xlab = "",ylab = "Peak counts", 
#                width = 0.5)
# ggpar(p, legend = "right", legend.title = "")
ggplot(data = df, aes(x=class, y=counts, fill=group)) + geom_bar(stat = "identity") +theme_classic() + scale_fill_brewer(palette = "Set2")

#  Annotation using ChIPpeakAnno
# aCR.sex <- assignChromosomeRegion(sex_peak, TxDb = atac.txdb)
# aCR.hypha <- assignChromosomeRegion(hypha_peak, TxDb = atac.txdb)
# library(tidyverse)
# aCR_hypha.df <- as.data.frame(aCR.hypha$percentage)
# aCR_sex.df <- as.data.frame(aCR.sex$percentage)


# make annoData for downstream analysis
annoData <- toGRanges(atac.txdb)

# if peaks located on promoter region of a genes, annotated it!
# sexPeak.anno_promoters <- annotatePeakInBatch(sex_peak, AnnotationData = annoData, 
#                                     output = "nearestBiDirectionalPromoters", bindingRegion = c(-1000, 500))

# pie1(table(sexPeak.anno_promoters$insideFeature))  # most annotated peaks are in upstream of genes

###  relationship with genes exp
exp <- read.table("/home/wanghm/whm/1FG/xulab_current/FG_expr.v2018")
sex_8d.exp <- exp[,c(33,34)]  # get sex-8d gene expression
sex_8d.exp$sex_8d <- rowMeans(sex_8d.exp)  # add new columns----average of R1 and R2
# merge diff-versions geneID together with gene.exp data
gene.map <- read.table("/home/wanghm/whm/1FG/xulab_current/gene_mapping.txt")
sex_8d_allgeneID.exp <- merge(sex_8d.exp,gene.map, by.x="row.names", by.y="Row.names")

gene.exp <- read.table("/home/wanghm/wanghm_R/1.res/expr_newID.txt", header = T)

sex.gr.filter <- sex.gr[sex.gr$distanceToTSS<300 & sex.gr$distanceToTSS>-500]
gene_score.df <- data.frame(newgeneID=sex.gr.filter$geneId, score=sex.gr.filter$score)
gene_score_uniq.df <- gene_score.df[!duplicated(gene_score.df$newgeneID),]

hypha.gr.filter <- hypha.gr[hypha.gr$distanceToTSS<300 & hypha.gr$distanceToTSS>-500]
gene_score.df2 <- data.frame(newgeneID=hypha.gr.filter$geneId, score=hypha.gr.filter$score)
gene_score_uniq.df2 <- gene_score.df2[!duplicated(gene_score.df2$newgeneID),]

atac_exp.df.sex <- merge(gene_score_uniq.df, gene.exp, by.x="newgeneID", by.y="new")
atac_exp.df.hypha <- merge(gene_score_uniq.df2, gene.exp, by.x="newgeneID", by.y="new")

atac_exp.df <- atac_exp.df.hypha[atac_exp.df.hypha$hypha >1,]
plot(scale(atac_exp.df$score), scale(atac_exp.df$hypha), xlim=c(-2,2), ylim=c(-2,2))

lm(sex_8d ~ score, atac_exp.df)
cor(atac_exp.df$hypha , atac_exp.df$score)

# get peak binding genes
sex_peakBindingGene_with_exp.merged_df <- merge(gene_score_uniq.df, sex_8d_allgeneID.exp, by.x="newgeneID", by.y="new")

# get other genes 
diff_geneIN_all <- setdiff(sex_8d_allgeneID.exp$new, sex_peakBindingGene_with_exp.merged_df$newgeneID)
nobinding_gene.df <- sex_8d_allgeneID.exp[sex_8d_allgeneID.exp$new %in% diff_geneIN_all,]

## plot gene.exp data
# check P-value
wilcox.test(sex_peakBindingGene_with_exp.merged_df$sex_8d, nobinding_gene.df$sex_8d, alternative = "greater")
boxplot(sex_peakBindingGene_with_exp.merged_df$sex_8d, nobinding_gene.df$sex_8d, ylim=c(0,100))

## diff-peak analysis: sex VS hypha24 
# load("/home/wanghm/wanghm_R/atac_chip/atac_diff.RData")
diff_edgeR.gr <- toGRanges("/home/wanghm/wanghm_R/atac_chip/diffPeak/sex_vs_hypha_edgeR.bed", format="MACS")
diffPeak.anno_promoter <- annotatePeakInBatch(a, AnnotationData = annoData, 
                    output = "nearestBiDirectionalPromoters", bindingRegion = c(-500, 300))

## go enrichment analysis and GSEA
go_enricher_analysis <- function(FGRRES_geneID){
  library(clusterProfiler)
  go.ids <- read.table("/home/wanghm/whm/go/FG.RR.27.GO.ids",skip = 1)
  func1 <- function(x){
    term <- strsplit(x[2], ",")[[1]]
    gene <- rep(x[1],length(term))
    subdataframe <- data.frame(term=term, gene=gene)
    return(subdataframe)
  }
  ss <- apply(go.ids, 1, func1)
  goID2gene <- do.call(rbind, ss)
  termid <- as.numeric(as.character(goID2gene[,1]))
  termid <- sprintf("GO:%07d", termid)        ##### complemented numeric to 7 wei with zero
  term2gene <- data.frame(term=termid, gene=goID2gene[,2])
  term2name <- read.delim("/home/wanghm/whm/go/term2name.txt", header = F, sep = "\t", stringsAsFactors = F)
  res <- enricher(gene = FGRRES_geneID, TERM2GENE = term2gene, TERM2NAME = term2name)
  # res <- GSEA(geneList = FGRRES_geneID, TERM2GENE = term2gene,pvalueCutoff = 0.8)
  return(res)
}

## extract target genes
target.gene <- unique()
sex_8d_allgeneID.exp2 <- sex_8d_allgeneID.exp[order(sex_8d_allgeneID.exp$sex_8d,decreasing = T),]
target_gene.fgrres <- sex_8d_allgeneID.exp2[sex_8d_allgeneID.exp2$new %in% target.gene,][,c(1,4)]
target_gene.fgrres <- target_gene.fgrres[!duplicated(target_gene.fgrres$Row.names),]
geneList <- target_gene.fgrres[,2]
names(geneList) <- as.character(target_gene.fgrres[,1])
# run GO enrichment
go.ids <- read.table("/home/wanghm/whm/go/FG.RR.27.GO.ids",skip = 1)

geneList2 <- geneList[names(geneList) %in% go.ids$V1]

ge <- read.delim("/home/wanghm/whm/ATAC/peak_calling/cal_rpkm/sp_gene.id")[,1]
ge_hypha <- read.delim("/home/wanghm/whm/ATAC/peak_calling/cal_rpkm/hysp_gene.id")[,1]
ge2 <- read.delim("/home/wanghm/whm/1FG/xulab_current/gene_mapping.txt")
gg <- ge2[ge2$new %in% ge,][,1]
hh <- ge2[ge2$new %in% ge_hypha,][,1]

target.gene <- sex_8d_allgeneID.exp[sex_8d_allgeneID.exp$new %in% a2[,1],]$Row.names
res <- go_enricher_analysis(target.gene)

gseaplot(res, geneSetID = geneList)
dotplot(res)
barplot(res)
enrichplot::cnetplot(res)
enrichplot::emapplot(res)
enrichplot::gseaplot(x = res, geneSetID =rownames(res[2,]))
res[2,]
# plot 
p1 <- enrichplot::dotplot(res)
p2 <- enrichplot::emapplot(res)
cowplot::plot_grid(p1, p2, align = "h")

##################################################################################

# Using data m6A_fra_df from script: strand_specific_m6A.R 

detags <- toGRanges("/home/wanghm/wanghm_R/atac_chip/diffPeak/hypha_sex_edgeR_detags.bed", format="MACS")

fast_slowGenome.gr <- read.table("/home/wanghm/whm/1FG/2speed-state.yellow-purple.txt", col.names = c("seqnames","start","end", "color"))%>% toGRanges()

## sp boxplot

aa <- read.delim("/home/wanghm/whm/ATAC/peak_calling_hmmr/cal_rpkm/sex_sp_gene.id")
a2 <- read.delim("/home/wanghm/whm/1FG/xulab_current/sex_specific_newID")
# length(intersect(a2$FG1G01250,aa$FG1G00090))
# length(intersect(aa$FG1G00090,a2$FG1G01250))

bb <- read.delim("/home/wanghm/whm/ATAC/peak_calling_hmmr/cal_rpkm/hysp_gene.id")
exp <- read.table("/home/wanghm/wanghm_R/1.res/expr_newID.txt", header = T)
bo1 <- exp[exp$new %in% bb[,1],][,3]
bo2 <- exp[!exp$new %in% bb[,1],][,3]
bo1.sex <- exp[exp$new %in% aa[,1],][,2]
bo2.sex <- exp[!exp$new %in% aa[,1],][,2]

box.df <- data.frame(tpm=c(bo1, bo2, bo1.sex, bo2.sex), 
                     class=c(rep("hypha", length(bo1)+length(bo2)), rep("sex", length(bo1.sex)+length(bo2.sex))),
                     type=c(rep("hypha_sp", length(bo1)), rep("hypha_no_sp", length(bo2)), rep("sex_sp", length(bo1.sex)), rep("sex_no_sp", length(bo2.sex))))
library(ggsignif)
options(scipen = 200)
ggplot(box.df, aes(x=type, y=tpm)) + geom_violin() + geom_boxplot(fill=c("#6a0f49", "#be95c4", "darkgreen", "lightgreen") ,width=0.25) + 
  scale_y_continuous(trans = "log10", breaks = c(1, 10, 100, 1000, 10000)) + ylab("TPM") + xlab("") + 
  geom_signif(comparisons = list(c("hypha_no_sp", "hypha_sp"), c("sex_no_sp", "sex_sp")),test = "wilcox.test", test.args = "less") + 
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(colour = "black", size = 1))
