## ChIP-seq data with INPUT and IP data
cat bw_input.conf | while read line;do
    array=($line)
    name=${array[2]}
    IP=${array[0]}
    IN=${array[1]}
    bamCompare -b1 $IP -b2 $IN --scaleFactorsMethod None -p 5 --binSize 1 --normalizeUsing RPKM --extendReads --ignoreDuplicates -o ${name}_log2Ratio.bw
done

## ATAC-seq data without INPUT bam

# make index
#ls ./rmPCR_duplicate/*_specific.bam |xargs -i samtools index {}

bamCoverage --normalizeUsing RPKM --extendReads -b ./rmPCR_duplicate/sex-7d_uniq_specific.bam -o ./sex_coverage.bw -p 5
bamCoverage --normalizeUsing RPKM --extendReads -b ./rmPCR_duplicate/hypha24h_uniq_specific.bam -o ./hypha_coverage.bw -p 5

## plot gene body atac-seq peaks distribution using deeptools

# TSS
#mkdir -p ./TSS
computeMatrix reference-point --referencePoint TSS -S ./sex_coverage.bw -R ./only_trans.gff -b 1200 -a 500 --skipZeros -o ./TSS/sex_mat --outFileSortedRegions ./TSS/sex-regions.bed
plotProfile -m ./TSS/sex_mat -o ./TSS/sex-profile_tss.pdf --dpi 720
plotHeatmap -m ./TSS/sex_mat -o ./TSS/sex-heatmap_tss.pdf --dpi 720 

computeMatrix reference-point --referencePoint TSS -S ./hypha_coverage.bw -R ./only_trans.gff -b 1200 -a 500 --skipZeros -o ./TSS/hypha_mat --outFileSortedRegions ./TSS/hypha-regions.bed
plotProfile -m ./TSS/hypha_mat -o ./TSS/hypha-profile_tss.pdf --dpi 720
plotHeatmap -m ./TSS/hypha_mat -o ./TSS/hypha-heatmap_tss.pdf --dpi 720 

# gene body 
#mkdir -p ./genebody
computeMatrix scale-regions -S ./sex_coverage.bw -R ./only_trans.gff -b 500 -a 500 --regionBodyLength 1000 --skipZeros -o ./genebody/sex_mat --outFileSortedRegions ./genebody/sex-region.bed
plotProfile -m ./genebody/sex_mat -o genebody/sex-profile_tss.pdf --dpi 720
plotHeatmap -m ./genebody/sex_mat -o genebody/sex-heatmap_tss.pdf --dpi 720

computeMatrix scale-regions -S ./hypha_coverage.bw -R ./only_trans.gff -b 500 -a 500 --regionBodyLength 1000 --skipZeros -o ./genebody/hypha_mat --outFileSortedRegions ./genebody/hypha-region.bed
plotProfile -m ./genebody/hypha_mat -o genebody/hypha-profile_tss.pdf --dpi 720
plotHeatmap -m ./genebody/hypha_mat -o genebody/hypha-heatmap_tss.pdf --dpi 720

# TES
#mkdir -p ./TES
computeMatrix reference-point --referencePoint TES -S ./sex_coverage.bw -R ./only_trans.gff -b 1200 -a 500 --skipZeros -o ./TES/sex_mat --outFileSortedRegions ./TES/sex-regions.bed
plotProfile -m ./TES/sex_mat -o ./TES/sex-profile_tss.pdf --dpi 720
plotHeatmap -m ./TES/sex_mat -o ./TES/sex-heatmap_tss.pdf --dpi 720 

computeMatrix reference-point --referencePoint TES -S ./hypha_coverage.bw -R ./only_trans.gff -b 1200 -a 500 --skipZeros -o ./TES/hypha_mat --outFileSortedRegions ./TES/hypha-regions.bed
plotProfile -m ./TES/hypha_mat -o ./TES/hypha-profile_tss.pdf --dpi 720
plotHeatmap -m ./TES/hypha_mat -o ./TES/hypha-heatmap_tss.pdf --dpi 720 


