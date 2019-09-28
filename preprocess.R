library(ChIPseeker)
### peaks of replicates
sex_rep1 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/sex-rep1_sorted_peaks.gappedPeak.bed")
sex_rep2 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/sex-rep2_sorted_peaks.gappedPeak.bed")
sex_rep3 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/sex-rep3_sorted_peaks.gappedPeak.bed")
sex_rep4 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/sex-rep4_sorted_peaks.gappedPeak.bed")

hypha_rep1 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/hypha24-rep1_sorted_peaks.gappedPeak.bed")
hypha_rep2 <- readPeakFile("/home/wanghm/whm/ATAC/peak_calling/open_bed/hypha24-rep2_sorted_peaks.gappedPeak.bed")

# overlap between replicates
library(ChIPpeakAnno)
makeVennDiagram(list(hypha_rep1, hypha_rep2), NameOfPeaks = c("rep1", "rep2"))
makeVennDiagram(list(sex_rep1, sex_rep2, sex_rep3, sex_rep4), NameOfPeaks = c("rep1", "rep2", "rep3", "rep4"))

# overlap between hypha and sexual
library(VennDiagram)
venn.diagram(x = list(hypha=1:5347, sex=1784:8655), filename = "sexVshypha_Vnn.pdf")

## peak length distribution
sex.peak <- read.table("/home/wanghm/whm/novo_atac/peak_calling/novo_sex_peaks.gappedPeak.bed", stringsAsFactors = F)
hypha.peak <- read.table("/home/wanghm/whm/ATAC/peak_calling_hmmr/open_bed/hypha_all_peaks.bed", stringsAsFactors = F)
peak.wid.sex <-  sex.peak$V3 - sex.peak$V2 + 1
mean(peak.wid.sex)
sex.mean <- 1000
peak.wid.hypha <- hypha.peak$V3 - hypha.peak$V2 + 1
mean(peak.wid.hypha)
hypha.mean <- 1314

mean <- data.frame(mean=c(sex.mean, hypha.mean), type=c("sex","hypha"))
library(ggpubr)
dd <- data.frame(width = c(peak.wid.sex, peak.wid.hypha), 
                 type = factor(c(rep("sex", length(peak.wid.sex)), rep("hypha", length(peak.wid.hypha)))))
# ggdensity(dd, x="width", fill = "type",color = "type", xlab = "Peak width", 
#           ylab = "Density", palette = "npg", xlim=c(0, 4000), add="mean")

ggplot(dd,aes(x=width, group=factor(type), fill=type, colour=type)) + geom_density(alpha=0.5)+
  geom_vline(data=mean, aes(xintercept=mean,color=type),linetype="dashed", size=1)+
  xlim(0,4000) + theme_classic() + xlab("Peak width") + ylab("Density")








