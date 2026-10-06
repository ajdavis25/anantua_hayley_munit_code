#!/bin/bash

#SBATCH --job-name=pbfix_dumps
#SBATCH --partition=compute1,compute2
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:30:00
#SBATCH --output=/work/vmo703/scratch/pb_fix_ab/slurm_logs/%x_%j.out
#SBATCH --error=/work/vmo703/scratch/pb_fix_ab/slurm_logs/%x_%j.err

set -euo pipefail

# Field-level A/B for the constant_beta_paper_literal fix: dump every zone's
# assigned Thetae (GRMONTY_DEBUG_DUMP_THETAE exits before any Monte Carlo)
# for the legacy (12*pi-hot) and paper-literal P_B forms on the MAD P3.1
# config. Compared afterwards against the independent numpy port.
# Live /work/vmo703/igrmonty/grmonty untouched; binary built in this scratch.

BIN=/work/vmo703/scratch/pb_fix_ab/build/grmonty_pbfix
DIR=/work/vmo703/scratch/pb_fix_ab

module purge
module load gcc/11.3.0
module load hdf5/1.12.0

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export LD_LIBRARY_PATH="/work/vmo703/aricarte/gsl-1.16/lib:/work/vmo703/aricarte/gsl-1.16/lib64:${LD_LIBRARY_PATH:-}"

for tag in legacy paperlit; do
  echo "[INFO] $(date) thetae dump: $tag"
  GRMONTY_DEBUG_DUMP_THETAE="$DIR/out/thetae_mad_wjet_${tag}.bin" \
    "$BIN" -par "$DIR/pars/mad_wjet_${tag}.par"
done

echo "[INFO] $(date) dumps done"
ls -la "$DIR/out"
