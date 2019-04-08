# make bam index
samtools index h3K14_Ac_TBI-3d_IN_sorted.bam
samtools index h3K14_Ac_TBI-3d_IP_sorted.bam

# bam -> bed
bamToBed -i sex-IN_sorted.bam > sex-IN_sorted.bed

# download khmer & calculate effective genome fraction(egf) for epic2 parament: -egf
# sudo apt-get install khmer
unique-kmers.py -k 100 xulab_ph1_rm_Mt.fasta > res  ## ph1 effective genome fraction: 366696207/37947603 = 0.96702306

# call broad peaks using epic2
for bam in h3K27_Ac/h3K27_Ac_sex_sorted.bam h3K4_Me/h3K4_Me_sex_sorted.bam h3K9_Ac/h3K9_Ac_sex_sorted.bam h3K9_Me/h3K9_Me_sex_sorted.bam;do
        # bamToBed -i $bam > ${bam%%.*}.bed
        epic2 -t ${bam%%.*}.bed -c sex-IN_sorted.bed -egf 0.967 --chromsizes chromsize.txt -o ${bam%%_sorted.bam}.peaks
done
