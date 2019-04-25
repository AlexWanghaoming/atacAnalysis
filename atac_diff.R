library(GenomicFeatures)
library(Rsamtools)
library(RColorBrewer)
library(rtracklayer)
library(ggplot2)
atac.txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")

library(ChIPseeker)
pmt <- getPromoters(atac.txdb, downstream = 1000, upstream = 1000)
# ge <- genes(atac.txdb)

peak1.atac <- readPeakFile("/home/wanghm/whm/ATAC/no_split/hypha.peaks")
peak2.atac <- readPeakFile("/home/wanghm/whm/ATAC/no_split/sex-7d.peaks")

ChIPseeker::covplot(peak1.atac, weightCol = "X1000")   # plot ATAC-seq signals genome distribution
## Tss distribution
tag_atac.mat1 <- getTagMatrix(peak1.atac, windows = pmt)
tag_atac.mat2 <- getTagMatrix(peak2.atac, windows = pmt)

# two plots in one panel
ll = list(tag_atac.mat1, tag_atac.mat2)
names(ll) <- c("hypha-24", "sex-7d")
plotAvgProf(ll, xlim = c(-1000, 1000))
# plotAvgProf(tag_atac.mat2, xlim = c(-1000, 1000))

### peaks of replicates
sex_rep1 <- readPeakFile("/home/wanghm/whm/ATAC/sex/sex-rep1_peaks.narrowPeak")
sex_rep2 <- readPeakFile("/home/wanghm/whm/ATAC/sex/sex-rep2_peaks.narrowPeak")
sex_rep3 <- readPeakFile("/home/wanghm/whm/ATAC/sex/sex-rep3_peaks.narrowPeak")
sex_rep4 <- readPeakFile("/home/wanghm/whm/ATAC/sex/sex-rep4_peaks.narrowPeak")

hypha_rep1 <- readPeakFile("/home/wanghm/whm/ATAC/hypha24/hypha24_peaks-rep1.narrowPeak")
hypha_rep2 <- readPeakFile("/home/wanghm/whm/ATAC/hypha24/hypha24_peaks-rep2.narrowPeak")

## plot width of peaks for each replicates
par(mfrow=c(2,1))
boxplot(width(sex_rep1), width(sex_rep2), width(sex_rep3), width(sex_rep4), 
        names=c("sex-rep1", "sex-rep2", "sex-rep3", "sex-rep4"), notch=T, ylim=c(0,1500), col="red")
boxplot(width(hypha_rep1), width(hypha_rep2), names=c("hypha-rep1", "hypha-rep2"), notch=T, ylim=c(0, 1500), col="green")

## diff peak analysis
library(DiffBind)
dbObj <- dba(sampleSheet = "/home/wanghm/wanghm_R/atac_chip/diffPeak/diffbind_sample.csv")
dbObj_count <- dba.count(dbObj, bUseSummarizeOverlaps = TRUE)

# check sample replicate 
dba.plotPCA(dbObj_count)

# get overlap peaks of replicates 
# dba.overlap(dbObj, dbObj$masks$hypha, mode = DBA_OLAP_RATE)   # calculate overlap
dba.plotVenn(dbObj, dbObj$masks$hypha)
dba.plotVenn(dbObj, dbObj$masks$sex)

# Establishing a contrast 
dbObj_contrast <- dba.contrast(dbObj_count, categories = DBA_FACTOR, minMembers = 2)
# identify diff peaks using deseq2, edgeR model respectly
dbObj_analysis <- dba.analyze(dbObj_contrast, method = DBA_ALL_METHODS)

# plot
dba.plotMA(dbObj_analysis, fold= 0, method = DBA_ALL_METHODS)
dba.plotVenn(dbObj_analysis, contrast = 1, method = DBA_ALL_METHODS)
dba.plotVolcano(dbObj_analysis, fold = 1)
dba.plotHeatmap(dbObj_analysis, contrast = 1, correlations = FALSE, scale="row")


## output diff-peaks identification results 
# ?dba.report
# dbObj_analysis$config print default th for dba.report
comp1.edgeR <- dba.report(dbObj_analysis,contrast = 1, method = DBA_EDGER)    # identify diff-peak FDR<0.05
comp1.deseq2 <- dba.report(dbObj_analysis,contrast = 1, method = DBA_DESEQ2)
out.edgeR <- as.data.frame(comp1.edgeR)
write.table(out.edgeR, file = "atac_chip/diffPeak/sex_vs_hypha_edgeR.txt", sep = "\t", quote = F, row.names = F)
out.deseq2 <- as.data.frame(comp1.deseq2)
write.table(out.deseq2, file = "atac_chip/diffPeak/sex_vs_hypha_deseq2.txt", sep = "\t", quote = F, row.names = F)

## save bed out 
edge.bed <- out.edgeR[, c("seqnames", "start", "end", "strand")]
write.table(edge.bed, file = "atac_chip/diffPeak/sex_vs_hypha_edgeR.bed", sep = "\t", quote = F, row.names = F)
deseq2.bed <- out.deseq2[, c("seqnames", "start", "end", "strand")]
write.table(deseq2.bed, file = "atac_chip/diffPeak/sex_vs_hypha_deseq2.bed", sep = "\t", quote = F, row.names = F)


