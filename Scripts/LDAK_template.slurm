#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Fri Sep 06 2019 16:37:28 GMT+1000 (Australian Eastern Standard Time)

# Partition for the job:
#SBATCH --partition=physical

# The name of the job:
#SBATCH --job-name="LDAK Sections"

# The project ID which this job should run under:
#SBATCH --account="punim0543"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1



# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL

# Use this email address:
#SBATCH --mail-user=andhikap@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-5:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The job command(s):
#need to do 3101 to 4533
./ldak5.linux --calc-weights sections --bfile k2029 --section FRAGMENT


k=0
for FRAG in 2583; do k=$(($k+1))
sed "s/FRAGMENT/$FRAG/" LDAK_template.slurm > $k.ladksection.slurm
done 


for TEMPERATURE in 0 3 5 7; do for Pth in 0 250 500 1000; do for N in 24; do for SNP in 0.1 K; do k=$(($k+1))
sed "s/RSQUARE/$SNP/" LMMLASSO_SLURM_template_minmax.slurm > $k.minmax.slurm
sed -i "s/TEMP/$TEMPERATURE/" $k.minmax.slurm
sed -i "s/PERIOD/$N/" $k.minmax.slurm
sed -i "s/PHOTO/$Pth/" $k.minmax.slurm
done
done
done
done

for i in *ladksection.slurm; do sbatch $i
done



