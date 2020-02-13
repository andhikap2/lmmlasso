#!/bin/bash
#Making microclimate files

for Th in 0 3 5 7
do 
	for N in 1 24 168 720
	do
		for PAR_Th in 0 250 500 1000
		do
			./MICROCLIMATE.R $Th $N $PAR_Th
		done
	done
done
