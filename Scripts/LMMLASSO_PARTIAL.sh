#!/bin/bash

for r2 in 0.1 0.2 0.3;  
do
	for temp in 7 5 3 0;
	do
		for period in 720 168 24 1;
		do
			./LMMLASSO_PIPELINE.py $r2 $temp $period
		done
	done
done
