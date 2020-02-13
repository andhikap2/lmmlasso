#!/bin/bash
#SBATCH --job-name=serial_job_test    # Job name
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --mem=100gb                   # Job memory request
#SBATCH --time=03:30:00               # Time limit hrs:min:sec
pwd; hostname; date
 
module load R

Rscript Your_model_script.R
 
date
