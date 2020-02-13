#!/bin/bash

r2=$1
if [ $r2 == K ]; then
	echo NO LD TRIMMING... USING K MATRIX INSTEAD

	#something about changing directories and running python and loading the K_Matrix and loading that as the X

else	
	cd ~/Clim_GWAS/PLINK
	sudo ./plink --bfile k2029 --geno 0.05 --maf 0.01 --indep-pairwise 500 100 $r2 --out SNPs_$r2
	sudo mv SNPs_$r2.prune.in ~/Clim_GWAS/Clim_GWAS_2/Temp_Files
	cd ~/Clim_GWAS/Clim_GWAS_2/Scripts
	./SNP_PROCESSING.py $1 

	#Then run the processing script to conform this to the right shape

fi

