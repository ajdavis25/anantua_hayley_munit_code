#!/bin/bash

#SBATCH --job-name=grmonty_p3test
#SBATCH --partition=compute1,compute2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --output=/work/vmo703/scratch/p3_jetcheck/slurm_logs/%x_%j.out
#SBATCH --error=/work/vmo703/scratch/p3_jetcheck/slurm_logs/%x_%j.err

set -euo pipefail

# Adapted from /work/vmo703/igrmonty/wjet_single_par.slurm (the project's established
# single-parfile template) for the session's D5/H2-fixed test binary. All paths are on
# /work (Lustre) deliberately -- an earlier version of this script staged the binary,
# parfiles, and log dir under the login node's /tmp, which is node-local (rootfs) and
# invisible to compute nodes; those jobs (778072/778073) were canceled before start.
# Right-sized for a Phase 3 validation run (Ns=1000, minutes), not production scale.
# The binary is a copy of build_archive/grmonty (md5 d36cbc97...) -- the live
# /work/vmo703/igrmonty/grmonty is untouched and still serves the tuning campaign.

GRMONTY_BIN="${GRMONTY_BIN:-/work/vmo703/scratch/p3_jetcheck/grmonty_h2d5}"

PAR_FILE="${PAR_FILE:?PAR_FILE env var must be set to an absolute parfile path}"
JOB_LABEL="${JOB_LABEL:-$(basename "${PAR_FILE%.*}")}"

if [[ ! -f "$PAR_FILE" ]]; then
  echo "[ERROR] PAR_FILE '$PAR_FILE' not found" >&2
  exit 1
fi

if [[ ! -x "$GRMONTY_BIN" ]]; then
  echo "[ERROR] GRMONTY NOT FOUND: $GRMONTY_BIN" >&2
  exit 1
fi

module purge
module load gcc/11.3.0
module load hdf5/1.12.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"

echo "[INFO] $(date) starting GRMONTY P3 test job '$JOB_LABEL' on $(hostname)"
echo "[INFO] PAR_FILE = $PAR_FILE"
echo "[INFO] GRMONTY_BIN = $GRMONTY_BIN"

"$GRMONTY_BIN" -par "$PAR_FILE"

echo "[INFO] $(date) finished GRMONTY P3 test job '$JOB_LABEL'"
