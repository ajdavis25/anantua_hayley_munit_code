#!/bin/bash
#SBATCH --job-name=critbeta_betacrit_0.1
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=10-00:00:00

#SBATCH --array=0-143%20
#SBATCH --output=/work/vmo703/logs/%x_%A_%a.out
#SBATCH --error=/work/vmo703/logs/%x_%A_%a.err

# --- safety ---
set -euo pipefail

# --- environment ---
module purge
module load gcc/11.3.0
source /work/vmo703/ipole_venv/bin/activate

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- workdir ---
cd /work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts

# --- run ---
python auto_munit_finder.py
