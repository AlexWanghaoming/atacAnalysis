## used for plot corration plot between multi bam files
multiBamSummary bins --bamfiles ./sex-rep1_sorted_uniq_specific.bam ./sex-rep2_sorted_uniq_specific.bam ./sex-rep3_sorted_uniq_specific.bam ./sex-rep4_sorted_uniq_specific.bam -o sex.npz -p 5
plotCorrelation --corData sex.npz -c spearman -p scatterplot -o sex_cor.pdf
