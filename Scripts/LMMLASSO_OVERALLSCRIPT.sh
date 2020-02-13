#!/bin/bash

#Parameters: $1= r2 threshold; $2= Th (induction temperature threshold), $3= N (monthly period) 


for period in 1 24 168 720
do
	./MICROCLIMATE.R $2 $period
	MICROCLIM_ID=$!

	./GENOTYPES.sh $1 &
	BACK_PID=$!

	wait $BACK_PID

	./SNP_PROCESSING.py $1
	SNP_PROCESSING_ID=$!

	wait

	./LMMLASSO_PIPELINE.py $1 $2 $period

done






'''./MICROCLIMATE.R $2 $3 &
MICROCLIM_ID=$!


./GENOTYPES.sh $1 &
BACK_PID=$!

wait $BACK_PID
./SNP_PROCESSING.py $1
SNP_PROCESSING_ID=$!

wait

./LMMLASSO_PIPELINE.py $1 $2 $3
'''
