#!/bin/bash
#SBATCH -N 2
#SBATCH -A FUA35_STELTURB
#SBATCH -p skl_fua_prod
#SBATCH --time 20:00:00
#SBATCH --job-name=stella
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bob.davies@york.ac.uk

NUMPROC=48

#cd ${SLURM_SUBMIT_DIR}
#echo "working on dir "$PATH

module purge
module load intel/pe-xe-2018--binary
module load intelmpi/2018--binary
module load zlib/1.2.8--gnu--6.1.0
module load szip/2.1--gnu--6.1.0
module load hdf5/1.10.4--intel--pe-xe-2018--binary
module load netcdf/4.6.1--intel--pe-xe-2018--binary
module load netcdff/4.4.4--intel--pe-xe-2018--binary
module load fftw/3.3.7--intelmpi--2018--binary
module load lapack
module load blas
module load mkl/2018--binary
# Bob: This seems to cause error
#module unload hdf5/1.10.4--intel--pe-xe-2018--binary


EXEFILE=~/stella_leapfrog/stella
#
OUTFILE=$SLURM_JOBID"_"$SLURM_JOB_NAME"_stella.out"
ERRFILE=$SLURM_JOBID"_"$SLURM_JOB_NAME"_stella.err"

mpirun -errfile-pattern $ERRFILE -outfile-pattern $OUTFILE -envall -genv -n $NUMPROC $EXEFILE input.in
