# Author: Wang Haoming 

# a <- ls()
# rm(list = a[a!="gAlignment.obj"])

library(ATACseqQC)
bamfile1 <- "/home/wanghm/whm/ATAC/hypha24h_atac_sorted.bam"
bamfile2 <-"/home/wanghm/whm/novo_atac/sex_atac_sorted_uniq_specific.bam"
bamfile.labels <- gsub(".bam", "", basename(bamfile2))

pdf(file = "libComplexity.pdf")
estimateLibComplexity(readsDupFreq(bamfile2))
dev.off()

# plot fragment size distribution 
environment(fragSizeDist.whm) <- environment(fragSizeDist)
fragSize <- fragSizeDist.whm(bamfile2, bamfile.labels)
fragSize <- fragSizeDist(bamfile2, bamfile.labels)
library(ATACseqQC)
## PT(promoter/transcripts) score
library(GenomicFeatures)
atac.txdb <-makeTxDbFromGFF(file = "/home/wanghm/whm/1FG/PH1-fusion-version2_represent_transcript_rm_Mt.gff")
txs <- transcripts(atac.txdb)

# pt <- PTscore(gAlignment.obj, txs)
# plot(pt$log2meanCoverage, pt$PT_score, xlab="log2 mean coverage", ylab="Promoter vs Transcript")
# 
# nfr <- NFRscore(gAlignment.obj, txs)
# plot(nfr$log2meanCoverage, nfr$NFR_score, 
#      xlab="log2 mean coverage",
#      ylab="Nucleosome Free Regions score",
#      main="NFRscore for 200bp flanking TSSs",
#      xlim=c(-10, 0), ylim=c(-5, 5))

## make BSgenome.fg.xulab and install it from tarball

# seed_files <- "/home/wanghm/whm/ATAC/bs_genome/BSgenome.fg-seed"
# BSgenome::forgeBSgenomeDataPkg("/home/wanghm/whm/ATAC/bs_genome/BSgenome.fg-seed")
## then quit R and turn to command line:
# cd to path of BSgenome.fg.xulab_v1 -> R CMD build <pkgdir> -> R CMD check <tarball> -> R CMD INSTALL <tarball>

## plot Footprints
library(MotifDb)
MCM1 <- query(MotifDb, c("HAP2"))
opi1 <- as.list(MCM1)
seqlev <- paste0("chr", 1:4)
library(BSgenome.fg.xulab.v1)

sigs <- factorFootprints(bamfile1, pfm = opi1[[1]], genome = BSgenome.fg.xulab.v1::BSgenome.fg.xulab.v1, 
                         min.score = "90%", seqlev = seqlev, upstream = 100, downstream = 100)

# load("/home/wanghm/wanghm_R/atac_chip/atacQC.RData")

