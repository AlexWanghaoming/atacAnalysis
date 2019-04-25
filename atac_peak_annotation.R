# a = ls()
# rm(list = a[!(a%in%c("dbObj", "dbObj_count", "dbObj_contrast", "dbObj_analysis", "comp1.deseq2", "comp1.edgeR"))])

library(GenomicFeatures)
atac.txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")

## plot peak profile
library(ChIPseeker)
pmt <- getPromoters(atac.txdb, downstream = 1000, upstream = 1000)

# get peak distribution matrix using score as weight
sex_peak.mat <- getTagMatrix(peak = sex_peak, windows = pmt, weightCol = "score")
hypha_peak.mat <- getTagMatrix(peak = hypha24_peak, windows = pmt, weightCol = "score")

# plot peak distribution profile around TSS
ll = list(hypha_peak.mat, sex_peak.mat)
names(ll) <- c("hypha-24", "sex-7d")
plotAvgProf(ll, xlim = c(-1000, 1000))

## peak annotation Using ChIPpeakAnno package
library(ChIPpeakAnno)
sex_peak <- toGRanges("/home/wanghm/whm/ATAC/rmPCR_duplicate/sex-7d_all/sex-7d_all_peaks.narrowPeak")
hypha_peak <- toGRanges("/home/wanghm/whm/ATAC/rmPCR_duplicate/hypha_all/hypha24_all_peaks.narrowPeak")
# sex_rep1 <- toGRanges("/home/wanghm/whm/ATAC/sex/peak/sex-rep1_peaks/sex-rep1_peaks.narrowPeak")
# sex_rep2 <- toGRanges("/home/wanghm/whm/ATAC/sex/peak/sex-rep2_peaks/sex-rep2_peaks.narrowPeak")
# sex_rep3 <- toGRanges("/home/wanghm/whm/ATAC/sex/peak/sex-rep3_peaks/sex-rep3_peaks.narrowPeak")
# sex_rep4 <- toGRanges("/home/wanghm/whm/ATAC/sex/peak/sex-rep4_peaks/sex-rep4_peaks.narrowPeak")

# ol <- findOverlapsOfPeaks(sex_peak, hypha24_peak)
# makeVennDiagram(ol)

# peak distribution on exon, intron, UTR, region
aCR <- assignChromosomeRegion(hypha24_peak, TxDb = atac.txdb)
barplot(aCR$percentage)

# make annoData for downstream analysis
annoData <- toGRanges(atac.txdb)

# if peaks located on promoter region of a genes, annotated it!
sexPeak.anno_promoters <- annotatePeakInBatch(sex_peak, AnnotationData = annoData, 
                                    output = "nearestBiDirectionalPromoters", bindingRegion = c(-1000, 500))

pie1(table(sexPeak.anno_promoters$insideFeature))  # most annotated peaks are in upstream of genes

###  relationship with genes exp
exp <- read.table("/home/wanghm/whm/1FG/FG_expr.v2018")
sex_8d.exp <- exp[,c(33,34)]  # get sex-8d gene expression
sex_8d.exp$sex_8d <- rowMeans(sex_8d.exp)  # add new columns----average of R1 and R2
# merge diff-versions geneID together with gene.exp data
gene_mapping1 <- read.table("/home/wanghm/whm/chip_seq/new_old_gene_mapping.txt", col.names = c("new", "old1"))
gene_mapping2 <- read.table("/home/wanghm/whm/1FG/rres2ramph1_gene_mapping.txt", col.names = c("rres", "ramph1"))
gene.id <- merge(gene_mapping1,gene_mapping2, by.x="old1", by.y="ramph1")
sex_8d_allgeneID.exp <- merge(sex_8d.exp, gene.id, by.x="row.names", by.y="rres")

gene_score.df <- data.frame(newgeneID=sexPeak.anno_promoters$feature, score=sexPeak.anno_promoters$score)
gene_score_uniq.df <- gene_score.df[!duplicated(gene_score.df$newgeneID),]

# get peak binding genes
sex_peakBindingGene_with_exp.merged_df <- merge(gene_score_uniq.df, sex_8d_allgeneID.exp, by.x="newgeneID", by.y="new")

# get other genes 
diff_geneIN_all <- setdiff(sex_8d_allgeneID.exp$new, sex_peakBindingGene_with_exp.merged_df$newgeneID)
nobinding_gene.df <- sex_8d_allgeneID.exp[sex_8d_allgeneID.exp$new %in% diff_geneIN_all,]

## plot gene.exp data
# check P-value
wilcox.test(sex_peakBindingGene_with_exp.merged_df$sex_8d, nobinding_gene.df$sex_8d, alternative = "greater")
boxplot(sex_peakBindingGene_with_exp.merged_df$sex_8d, nobinding_gene.df$sex_8d, ylim=c(0,100))

# library(squash)
# plot(scale(sex_peakBindingGene_with_exp.merged_df$sex_8d),
#       scale(sex_peakBindingGene_with_exp.merged_df$score), xlim = c(0,5))
# 
# dd <- data.frame(exp=scale(sex_peakBindingGene_with_exp.merged_df$sex_8d), 
#                  score=scale(sex_peakBindingGene_with_exp.merged_df$score))
# m = lm(score ~ exp, data = dd)
# plot(score ~ exp, data = dd)
# abline(m)
# summary(m)

## diff-peak analysis: sex VS hypha24 
# load("/home/wanghm/wanghm_R/atac_chip/atac_diff.RData")
diff_deseq2.gr <- toGRanges("/home/wanghm/wanghm_R/atac_chip/diffPeak/sex_vs_hypha_deseq2.bed", format="MACS")
diff_edgeR.gr <- toGRanges("/home/wanghm/wanghm_R/atac_chip/diffPeak/sex_vs_hypha_edgeR.bed", format="MACS")

a <- subsetByOverlaps(diff_deseq2.gr, diff_edgeR.gr)   # get intersection of edgeR and deseq2

diffPeak.anno_promoter <- annotatePeakInBatch(a, AnnotationData = annoData, 
                    output = "nearestBiDirectionalPromoters", bindingRegion = c(-1000, 500))

## go enrichment analysis
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
  return(res)
}
## extract target genes
target.gene <- unique(diffPeak.anno_promoter$feature)
target_gene.fgrres <- sex_8d_allgeneID.exp[sex_8d_allgeneID.exp$new %in% target.gene,]$Row.names

# run GO enrichment
res <- go_enricher_analysis(target_gene.fgrres)

# plot 
p1 <- enrichplot::dotplot(res)
p2 <- enrichplot::emapplot(res)
cowplot::plot_grid(p1, p2, align = "h")

##################################################################################

# Using data m6A_fra_df from script: strand_specific_m6A.R 





