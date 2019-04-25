#!/bin/bash

# split paired-end bam randomly via split_bam.py
python3 split_bam -d ./ -i sex-7d_atac_sorted.bam -f 0.5
