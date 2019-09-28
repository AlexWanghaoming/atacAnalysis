### used for plot with R circlize package
## calculate bin's reads coverage in windows
BIN_SIZE=20000
for i in ../sex/h3_IP_sorted.bam ../sex/h3K27_Ac_sex_sorted.bam ../sex/h3K27_Me_IP_sorted.bam ../sex/h3K4_Me_sex_sorted.bam ../sex/h3K9_Ac_sex_sorted.bam ../sex/h3K9_Me_sex_sorted.bam; do
    mosdepth -t 4 -b 20000 ${i%_sorted*} ${i}
    gzip -d ${i%_sorted*}.regions.bed.gz
done
