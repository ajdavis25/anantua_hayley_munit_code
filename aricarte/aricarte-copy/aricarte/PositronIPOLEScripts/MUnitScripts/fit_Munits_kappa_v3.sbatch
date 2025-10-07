#!/bin/bash
#SBATCH -n 36 # Number of cores requested
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH -t 3-0:00:00 # Runtime in minutes
#SBATCH -p narayan # Partition to submit to
#SBATCH --mem=8000 # Memory per cpu in MB (see also --mem-per-cpu)
#SBATCH -o /n/home11/aricarte/slurm_output/%j.out          # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /n/home11/aricarte/slurm_output/%j.err          # File to which STDERR will be written, %j inserts jobid
#SBATCH --mail-type=ALL                                             # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=angelo.ricarte@cfa.harvard.edu                  # Email to which notifications will be sent

export OMP_NUM_THREADS=48

module load Anaconda/5.0.1-fasrc02
module load hdf5
module load GCC/8.2.0-2.31.1 GSL/2.5

python fit_Munits_kappa_v3.py ${1} ${2} ${3}
