# bamCompare -b1 /home/wanghm/whm/pacbio_data/BAM-27/B2_Ip.bam -b2 /home/wanghm/whm/pacbio_data/BAM-27/B2_In.bam --extendReads -o B2_log2Ratio.bw -p 5
#computeMatrix scale-regions -S B2_log2Ratio.bw S2_log2Ratio.bw -R /home/wanghm/whm/1FG/FG_noCDS.gtf -b 500 -a 500 -o B2&S2Matrix --regionBodyLength 1830 --outFileSortedRegions S2.bed

## plot gene body atac-seq peaks distribution using deeptools

bamCoverage -b ./hypha24h_atac_sorted.bam -o hypha_coverage.bw -p 5
computeMatrix scale-regions -S hypha_coverage.bw -R ./only_trans.gff -b 1200 -a 1200 --regionBodyLength 1000 -o hypha_mat
plotProfile -m hypha_mat -o hypha.pdf --dpi 720 --yMin 100 --yMax 200

bamCoverage -b ./sex-7d_atac_sorted.bam -o sex-7d_coverage.bw -p 5
computeMatrix scale-regions -S sex-7d_coverage.bw -R ./only_trans.gff -b 1200 -a 1200 --regionBodyLength 1000 -o sex-7d_mat
plotProfile -m sex-7d_mat -o sex-7d.pdf --dpi 720
