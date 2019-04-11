#!/bin/bash
for i in `ls output*.bam`;do
    samtools sort ${i} -n -@ 4 -o ${i%%_sorted.bam}_sort.bam
    Genrich -t ${i%%_sorted.bam}_sort.bam -o ${i%%_sort.bam}.peaks -r -j
done
