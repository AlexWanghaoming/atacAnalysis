#!/bin/bash

# ATAC-seq Analysis Pipeline
# Author: Wang Haoming

## many ATAC-seq reads may mapping to mitochondria, so remove Mt sequence firstly
awk '/>/{p=($0!~/Mt/)} p' xulab_ph1.fasta > xulab_ph1_rm_Mt.fasta
samtools faidx xulab_ph1_rm_Mt.fasta

## build index
bowtie2-build xulab_ph1_rm_Mt.fasta xulab_ph1_rm_Mt.fasta

## bowtie2 alignment
conf_filename=`ls *.conf`
if [ -f $conf_filename ];then
        echo "reading fastq from conf"
else
        echo "Error: conf file not exit!"
        exit 1
fi

cat mapping_input.conf | while read line;do
array=($line)
NAME=${array[2]}
LEFT=${array[0]}
RIGHT=${array[1]}

# remove adapters if do not have adapter.fasta
NGmerge -a -u 41 -n 15 -1 $LEFT -2 $RIGHT -o $NAME

bowtie2 --very-sensitive -x xulab_ph1_rm_Mt.fasta -p 15 -1 ${NAME}_1.fastq.gz -2 ${NAME}_2.fastq.gz | samtools view -bS | samtools sort -o ${NAME}_sorted.bam

done