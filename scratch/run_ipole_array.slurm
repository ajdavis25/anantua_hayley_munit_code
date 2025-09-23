#!/bin/bash
#SBATCH --job-name=ipole_array
#SBATCH --output=/work/vmo703/slurm_logs/ipole_array_%A_%a.out
#SBATCH --error=/work/vmo703/slurm_logs/ipole_array_%A_%a.err
#SBATCH --partition=compute1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-00:00:00  # adjust to your liking
#SBATCH --array=0-<N-1>%10  # re-insert this dynamically

cd /work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts

# Activate your env if needed
source /work/vmo703/ipole_venv/bin/activate

# Run Python script with array index
python3 render_ipole_outputs.py $SLURM_ARRAY_TASK_ID
