#!/bin/bash

#SBATCH --job-name=pbfix_spec
#SBATCH --partition=compute1,compute2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=/work/vmo703/scratch/pb_fix_ab/slurm_logs/%x_%j.out
#SBATCH --error=/work/vmo703/scratch/pb_fix_ab/slurm_logs/%x_%j.err

set -euo pipefail

# Spectrum-level A/B for constant_beta_paper_literal: full Ns=1e6 grmonty run
# of one par file (PAR_FILE env), same seed/M_unit/config as its pair partner.
# Adapted from the project's single-parfile template (p3_test.slurm).

BIN=/work/vmo703/scratch/pb_fix_ab/build/grmonty_pbfix
PAR_FILE="${PAR_FILE:?PAR_FILE env var must be set to an absolute parfile path}"

if [[ ! -f "$PAR_FILE" ]]; then
  echo "[ERROR] PAR_FILE '$PAR_FILE' not found" >&2
  exit 1
fi

module purge
module load gcc/11.3.0
module load hdf5/1.12.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"

echo "[INFO] $(date) starting pbfix spectrum run: $PAR_FILE"
"$BIN" -par "$PAR_FILE"
echo "[INFO] $(date) finished pbfix spectrum run: $PAR_FILE"
