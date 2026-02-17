#!/bin/bash
#SBATCH --job-name=critbeta_betacrit
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=80
#SBATCH --cpus-per-task=1
#SBATCH --time=3-00:01:00
#SBATCH --output=/work/vmo703/logs/%x_%j.out
#SBATCH --error=/work/vmo703/logs/%x_%j.err

# --- safety ---
set -euo pipefail

# --- environment ---
module purge
module load gcc/11.3.0
source /work/vmo703/ipole_venv/bin/activate

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# --- workdir ---
cd /work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts

# --- run ---
mkdir -p /work/vmo703/logs

INPUT_CSV="/work/vmo703/data/paper_data_critbeta_both.csv"
export INPUT_CSV
TOTAL_ROWS="$(python - <<'PY'
import csv
import os

path = os.environ["INPUT_CSV"]
with open(path, "r", newline="") as fh:
    n = sum(1 for _ in csv.reader(fh)) - 1
print(max(n, 0))
PY
)"

if [[ "${TOTAL_ROWS}" -le 0 ]]; then
  echo "[error] no data rows found in ${INPUT_CSV}"
  exit 1
fi

MAX_CONCURRENT="${SLURM_NTASKS:-80}"
LAST_ROW="$((TOTAL_ROWS - 1))"

echo "[info] launching rows 0..${LAST_ROW} with up to ${MAX_CONCURRENT} concurrent workers"

for ROW in $(seq 0 "${LAST_ROW}"); do
  srun --exclusive -N1 -n1 -c1 \
    python auto_munit_finder.py --row "${ROW}" \
    > "/work/vmo703/logs/row_${SLURM_JOB_ID}_${ROW}.out" \
    2> "/work/vmo703/logs/row_${SLURM_JOB_ID}_${ROW}.err" &

  while (( $(jobs -pr | wc -l) >= MAX_CONCURRENT )); do
    sleep 1
  done
done

wait
echo "[info] all rows finished"
