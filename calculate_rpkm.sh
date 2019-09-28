#!/bin/bash

# used for calculating libsize of ATAC-seq and rpkm of peaks
#HYPHA_LIB_SIZE=samtools idxstats ../../rmPCR_duplicate/hypha24h_uniq_specific.bam | cut -f 3 | awk '{a=a+$0}END{print a}'
#SEX_LIB_SIZE=samtools idxstats ../../rmPCR_duplicate/sex-7d_uniq_specific.bam | cut -f 3 | awk '{a=a+$0}END{print a}'

## count histone chip-seq reads on ATAC-seq peaks
featureCounts -t sequence_feature -g geneId -f -O -T 5 --minOverlap 20 -a ./sex_peakAnno.gff -o sex_readCount /home/wanghm/whm/novo_atac/sex_atac_sorted_uniq_specific.bam \
    /home/wanghm/whm/chip_seq/sex-7d/h3_IP_sorted.bam /home/wanghm/whm/chip_seq/sex-7d/h3K27_Ac_sex_sorted.bam /home/wanghm/whm/chip_seq/sex-7d/h3K27_Me_IP_sorted.bam /home/wanghm/whm/chip_seq/sex-7d/h3K9_Ac_sex_sorted.bam /home/wanghm/whm/chip_seq/sex-7d/h3K9_Me_sex_sorted.bam /home/wanghm/whm/chip_seq/sex-7d/h3K4_Me_sex_sorted.bam
