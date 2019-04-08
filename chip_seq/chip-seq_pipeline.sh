#!/bin/bash

# chip-seq Analysis Pipeline
# Author: Wang Haoming

## remove Mt sequence
awk '/>/{p=($0!~/Mt/)} p' xulab_ph1.fasta > xulab_ph1_rm_Mt.fasta
samtools faidx xulab_ph1_rm_Mt.fasta

## build index
bowtie2-build xulab_ph1_rm_Mt.fasta xulab_ph1_rm_Mt.fasta

cat mapping_Input.conf | while read line;do
array=($line)
name=${array[2]}
L=${array[0]}
R=${array[1]}
L_paired=${L%.fastq.gz}-paired.fastq
L_unpaired=${L%.fastq.gz}-unpaired.fastq
R_paired=${R%.fastq.gz}-paired.fastq
R_unpaired=${R%.fastq.gz}-unpaired.fastq
## filter raw data 
java -jar /share/nas3/xujr/local_software/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 15 -phred33 $L $R $L_paired $L_unpaired $R_paired $R_unpaired LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20

## bowtie2 alignment
conf_filename=`ls *.conf`
if [ -f $conf_filename ];then
	echo "reading fastq from conf"
else
	echo "Error: conf file not exit!"	
	exit 1
fi

bowtie2 -x xulab_ph1_rm_Mt.fasta -p 15 -1 $L_paired -2 $R_paired -U $L_unpaired -U $R_unpaired -S ${name}.sam

## sam -> bam & sort bam file
samtools view -bS -@ 12 ${name}.sam | samtools sort -@ 12 -o ${name}_sorted.bam
done 
