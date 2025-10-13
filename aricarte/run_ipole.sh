#!/usr/bin/env bash
# make sure GSL libraries are available
export LD_LIBRARY_PATH=/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:$LD_LIBRARY_PATH

# set OpenMP threads to what Slurm gave us
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "[$(date)] starting ipole with $OMP_NUM_THREADS threads"

# run ipole
/work/vmo703/aricarte/aricarte-copy/aricarte/ipole+e-/ipole "$@"

echo "[$(date)] finished ipole"
