#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Fri Aug 23 2019 17:36:12 GMT+1000 (Australian Eastern Standard Time)

# Partition for the job:
#SBATCH --partition=physical

# The name of the job:
#SBATCH --job-name="LMMLASSO"

# The project ID which this job should run under:
#SBATCH --account="punim0543"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --ntasks=1


# The amount of memory in megabytes per process in the job:
#SBATCH --mem-per-cpu=81920



# Send yourself an email when the job:
# aborts abnormally (fails)
#SBATCH --mail-type=FAIL


# Use this email address:
#SBATCH --mail-user=andhikap@student.unimelb.edu.au

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-12:0:00

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

# The modules to load:
module load Python/3.7.3-spartan_gcc-8.1.0

# The job command(s):
python3 keras_dl_SPARTAN.py HL1 HL2 HL3 HL4 EP 