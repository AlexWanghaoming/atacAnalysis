#!/bin/bash
# rm h3 noise in histone chip seq of h3k27ac, h3k9ac, h3k9me, h3k4me ...
cat bw_input.conf | while read line;do
    array=($line)
    name=${array[2]}
    IP=${array[0]}
    IN=${array[1]}
    bamCompare -b1 $IP -b2 $IN --scaleFactorsMethod None --operation ratio -p 5 -bs 1 --normalizeUsing RPKM --extendReads --ignoreDuplicates --minMappingQuality 20 -o ${name}_1bp_ratio.bw
done
# remove noise and calculate ratio in bins
bigwigCompare -b1 sex-h3k27AC_1bp_ratio.bw -b2 sex-h3_1bp_ratio.bw --operation ratio -bs 20000 -p 4 -o h3k27Ac_rmH3.bedgraph -of bedgraph
