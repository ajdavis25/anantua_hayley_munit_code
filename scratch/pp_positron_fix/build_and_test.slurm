#!/bin/bash

#SBATCH --job-name=pp1_positron_tests
#SBATCH --partition=compute1,compute2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/work/vmo703/scratch/pp_positron_fix/slurm_logs/%x_%j.out
#SBATCH --error=/work/vmo703/scratch/pp_positron_fix/slurm_logs/%x_%j.err

set -euo pipefail

# Builds the staged igrmonty tree (PP-1 e+e- brems prefactor fix + tightened
# pair tests) and runs the internal regression suite. The suite signals
# all-tests-passed by exit(42). Staged copy keeps the repo's live ./grmonty
# binary untouched.

STAGE=/work/vmo703/scratch/pp_positron_fix/igrmonty

module purge
module load gcc/11.3.0
module load hdf5/1.12.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"

cd "$STAGE"
echo "[INFO] $(date) building on $(hostname), tree $(git describe --dirty --always 2>/dev/null || echo unknown)"
make clean >/dev/null 2>&1 || true
make -j"$SLURM_CPUS_PER_TASK"

echo "[INFO] $(date) running --run-tests"
set +e
./grmonty --run-tests
rc=$?
set -e

if [[ $rc -ne 42 ]]; then
  echo "[ERROR] test harness exited rc=$rc (expected 42 = all tests passed)" >&2
  exit 1
fi
echo "[INFO] $(date) all tests passed (rc=42)"
