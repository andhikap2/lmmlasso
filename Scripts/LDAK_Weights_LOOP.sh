''' Get weightings'''

cd ~/Clim_GWAS/PLINK
for i in {3001..3100}
do
	./ldak5.linux --calc-weights sections --bfile k2029 --section $i
done


#using parallel
parallel --gnu -j 32 './ldak5.linux --calc-weights sections --bfile k2029 --section {}' ::: `echo {1..1000}`

#Joining weights across sections
./ldak5.linux --bfile k2029 --join-weights sections/

#!/bin/bash
#$ -S /bin/bash
#$ -t 1-4533

number=$SGE_TASK_ID
../ldak5.linux --calc-weights sections --bfile k2029 --section $number

''' Get kinships '''

./ldak5.linux --bfile k2029 
./ldak5.linux --calc-kins-direct k2029_kinship_matrix_neg0.25 --bfile k2029 --weights weights.all --power -0.25 --kinship-raw YES #try out other power values as well


./ldak5.linux --calc-kins-direct k2029_kinship_matrix_neg0.25_chr1 --bfile k2029 --chr 1 --weights weights.all --power -0.25 --kinship-raw YES #calculate a chr1 kinship [do this after you figure out which power to use]
#APNO

