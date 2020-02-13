#!/bin/bash

touch Genes1000.gene
for LINE in $(cat SNP_ids_1000.txt)
do
echo $LINE > temp.temp
CHROM=$(cut -d_ -f1 temp.temp)
POS=$(cut -d_ -f2 temp.temp)
awk -v pos=$POS -v chrom=$CHROM '$1 == chrom  && $4 > pos-1000 && $5 < pos+1000 {print$10}' Arabidopsis_thaliana.TAIR10.43.gene >> Genes1000.gene
echo $LINE
done
