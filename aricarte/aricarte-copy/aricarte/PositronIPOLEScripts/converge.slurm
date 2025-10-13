#!/bin/bash
#SBATCH --job-name=ipole_converge
#SBATCH --output=/work/vmo703/logs/converge_%A_%a.out
#SBATCH --error=/work/vmo703/logs/converge_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=7-00:00:00
#SBATCH --partition=compute2
#SBATCH --array=0-11%10

# Load your Python venv
source /work/vmo703/ipole_venv/bin/activate

# Path to your CSV
CSV=/work/vmo703/data/example.csv

# Figure out row count (skip 1 header rows)
NROWS=$(($(wc -l < $CSV) - 1))

# Set threads to match allocated cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "Running row $SLURM_ARRAY_TASK_ID out of $NROWS"

# Run Python with row index
python /work/vmo703/aricarte/aricarte-copy/aricarte/PositronIPOLEScripts/shaum_munit_finder.py --start_row $SLURM_ARRAY_TASK_ID
