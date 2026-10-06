#!/bin/bash
#SBATCH --job-name=p4_confirm_sane
#SBATCH --partition=anantuabhg
#SBATCH --account=anantuabhg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=2-00:00:00
#SBATCH --output=/work/vmo703/igrmonty_logs/%x_%j.out
#SBATCH --error=/work/vmo703/igrmonty_logs/%x_%j.err
set -euo pipefail
module purge
module load gcc/11.3.0 hdf5/1.12.0
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
/work/vmo703/scratch/p4_staging/build/grmonty_p4b -par /work/vmo703/scratch/p4_staging/confirm/confirm_bias0.5.par
