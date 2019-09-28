#!/bin/bash
## used for deeptools plot profile: expr_newID.txt is a dataframe with three columns, we sort it by gene exp and split it equally and output three levels of gff

#sort -k3n,3 expr_newID.txt | cut -f 1 | split -l 4019

level=(low media high)
k=0
for i in `ls xa*`;do
    awk -F "\"" 'ARGIND==1{a[$1]=$1}ARGIND==2{if($2 in a){print $0}}' ${i} ./only_trans.gff > ${level[k]}_exp_hypha.gff
    k=$k+1
done
