#!/bin/bash

#Parameters: $1= r2 threshold; $2= Th (induction temperature threshold), $3= N (monthly period) 



for r2 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
do
	./GENOTYPES.sh $r2 &
	echo Running GENOTYPES.sh
	BACK_PID=$!
	wait $BACK_PID

	./SNP_PROCESSING.py $r2&
	SNP_PROCESSING_ID=$!
done


for period in 1 24 168 720
do
	for temp in 0 3 5 7
	do
		./MICROCLIMATE.R $temp $period &
		echo Running MICROCLIMATE.R
	done
done
q


wait


for period in 1 24 168 720
do
	for temp in 0 3 5 7
	do
		for r2 in 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9
		do
			./LMMLASSO_PIPELINE.py $r2 $temp $period
		done
	done
done
