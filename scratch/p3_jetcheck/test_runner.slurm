#!/bin/bash

#SBATCH --job-name=grmonty_d6tests
#SBATCH --partition=compute1,compute2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:20:00
#SBATCH --output=/work/vmo703/scratch/p3_jetcheck/slurm_logs/%x_%j.out
#SBATCH --error=/work/vmo703/scratch/p3_jetcheck/slurm_logs/%x_%j.err

set -euo pipefail

# Runs grmonty's internal regression suite (--run-tests), which now includes the
# D6 deep-KN Compton sampler ground-truth tests (src/tests.c,
# test_deep_kn_compton_sampling). The harness signals success by exit(42).

GRMONTY_BIN="${GRMONTY_BIN:-/work/vmo703/scratch/p3_jetcheck/grmonty_d6}"

if [[ ! -x "$GRMONTY_BIN" ]]; then
  echo "[ERROR] GRMONTY NOT FOUND: $GRMONTY_BIN" >&2
  exit 1
fi

module purge
module load gcc/11.3.0
module load hdf5/1.12.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"

echo "[INFO] $(date) running grmonty --run-tests on $(hostname) with $GRMONTY_BIN"

set +e
"$GRMONTY_BIN" --run-tests
rc=$?
set -e

if [[ $rc -ne 42 ]]; then
  echo "[ERROR] test harness exited rc=$rc (expected 42 = all tests passed)" >&2
  exit 1
fi
echo "[INFO] $(date) all tests passed (rc=42)"
