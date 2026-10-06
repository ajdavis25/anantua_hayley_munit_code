#!/bin/bash
#SBATCH --job-name=p4_frames
#SBATCH --partition=compute1,compute2
#SBATCH --account=anantuabhg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --array=0-23%8
#SBATCH --output=/work/vmo703/igrmonty_logs/%x_%A_%a.out
#SBATCH --error=/work/vmo703/igrmonty_logs/%x_%A_%a.err
set -euo pipefail
module purge
module load gcc/11.3.0 hdf5/1.12.0
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
HERE=/work/vmo703/ipole_outputs/M87/frames_p4
BRANCH=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$HERE/branches.txt")
echo "[task $SLURM_ARRAY_TASK_ID] branch $BRANCH"
while read -r par; do
  img=$(grep '^outfile' "$par" | awk '{print $2}')
  if [ -s "$img" ]; then echo "[skip] $img exists"; continue; fi
  echo "[render] $par"
  /work/vmo703/ipole/ipole/ipole -par "$par" 2>&1 | tail -1
done < "$HERE/pars/$BRANCH.list"
echo "[done] $BRANCH"
