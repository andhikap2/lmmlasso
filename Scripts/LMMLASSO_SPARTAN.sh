#!/bin/sh


for period in 1 24 168 720
do
	for temp in 0 3 5 7
	do
		for par in 0 250 500 1000
		do
			for r2 in K 0.1
			do
				python3 LMMLASSO_PIPELINE_SPARTAN.py $r2 $temp $period $par
			done
		done
	done
done
