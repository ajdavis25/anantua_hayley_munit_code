#!/bin/bash
#SBATCH --job-name=ipole_array
#SBATCH --output=/work/vmo703/slurm_logs/ipole_array_%A_%a.out
#SBATCH --error=/work/vmo703/slurm_logs/ipole_array_%A_%a.err
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00  # adjust to your liking
#SBATCH --array=0-95%10  # adjust based on lines of data in csv

cd /work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts

# Activate your env if needed
source /work/vmo703/ipole_venv/bin/activate

echo "[$(date)] SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID"

# Run Python script with array index
python3 ipole_positron_test_unlocked.py $SLURM_ARRAY_TASK_ID
