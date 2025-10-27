#!/bin/bash
#SBATCH --job-name=ipole_converge
#SBATCH --output=/work/vmo703/logs/converge_%A_%a.out
#SBATCH --error=/work/vmo703/logs/converge_%A_%a.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=80
#SBATCH --time=7-00:00:00
#SBATCH --partition=compute2
#SBATCH --array=0-184%10

# load your python venv
source /work/vmo703/ipole_venv/bin/activate

# path to your csv
CSV=/work/vmo703/data/sgra_batch.csv

# figure out row count (skip 1 header rows)
NROWS=$(($(wc -l < $CSV) - 1))

# set threads to match allocated cores
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "running row $SLURM_ARRAY_TASK_ID out of $NROWS"

# run python with row index
python /work/vmo703/SgrA/auto_munit_finder_sgrA.py --start_row $SLURM_ARRAY_TASK_ID
