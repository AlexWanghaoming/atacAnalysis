#!/bin/bash
# Author: Wang haoming
#GENOME=~/whm/ATAC/rmPCR_duplicate/hmmrATAC/genome.info
#BAM=(`ls ~/whm/ATAC/*/*uniq_specific.bam`)
#BAI=(`ls ~/whm/ATAC/*/*uniq_specific.bam.bai`)
#
#for ((i=0; i<${#BAM[@]} && i<${#BAI[@]}; i++)); do
#    kk=`basename ${BAM[$i]}`
#    prefix=${kk%_uniq*}
#    java -jar ~/whm/ATAC/rmPCR_duplicate/hmmrATAC/HMMRATAC_V1.2.5_exe.jar -u 20 -l 2 --score fc -o ${prefix} -b ${BAM[$i]} -i ${BAI[$i]} -g ${GENOME} 
#done

for i in ./*.gappedPeak;do
    cut -f 1,7,8,13 ${i} > open_bed/${i}.bed
done
#!/bin/bash
#cat sex-rep* | sort -k1,1 -k2,2n | mergeBed -i stdin > location.bed

#for i in sex-rep*;do
#    PREFIX=${i%_sorted*}
#    awk -v p="$PREFIX" '{print $0"\t", p"-"NR}' ${i} > ${i}_id
#done

# merge peaks and label from which they merge
cat sex*_id | sort -k1,1 -k2,2n | mergeBed -i stdin -o collapse -c 4,5 > sex_location.bed
