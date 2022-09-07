#!/bin/bash

#SBATCH --job-name=basic_job_test        # Job name
#SBATCH --mail-type=BEGIN,END,FAIL             # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=bob.davies@york.ac.uk     # Where to send mail  
#SBATCH --ntasks=24                       # Run on a single CPU
#SBATCH --time=12:00:00                  # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --output=array_job_%A_%a.log     # Standard output and error log
#SBATCH --account=phys-gspt-2019         # Project account

sim_name="input.in"
mpirun -n 24 /users/rd1042/scratch/stella/stella $sim_name

