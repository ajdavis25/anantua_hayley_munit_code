#!/usr/bin/env python3
import sys
import pandas as pd
import numpy as np
from pathlib import Path
from ipole_positron_test_unlocked import makePositronImages
from ipole_many_models import runIPOLE  # already imported inside

# --- your model lookup table ---
model_settings = {
    "RBETA":        {"electronModel": 2, "sigma_transition": 2.0},
    "RBETAWJET":    {"electronModel": 2, "sigma_transition": 1.0},
    "CRITBETA":     {"electronModel": 4, "sigma_transition": 2.0},
    "CRITBETAWJET": {"electronModel": 4, "sigma_transition": 1.0},
}

# --- load CSV with metadata ---
csv_path = Path("/work/vmo703/data/munits_table.csv")
df = pd.read_csv(csv_path, skiprows=2)
df.columns = df.columns.str.strip()  # clean up spaces
print("CSV columns:", df.columns.tolist())

# --- where all your dump files live ---
dump_dir = Path("/work/vmo703/grmhd_dump_samples/")

# --- check if we are running under SLURM array ---
task_id = None
if len(sys.argv) > 1:
    try:
        task_id = int(sys.argv[1])
    except ValueError:
        print(f"⚠️ Invalid task ID argument: {sys.argv[1]}")
        sys.exit(1)

files = sorted(dump_dir.glob("*.h5"))

if task_id is not None:
    # Only pick the one file corresponding to this task_id
    files = [files[task_id]]
    print(f"SLURM task {task_id}: running on {files[0]}")

# --- loop over all dump files ---
for simFile in files:
    # pull timestep number from filename (assumes dump_00005000.h5 etc.)
    try:
        timestep = int(simFile.stem.split("_")[-1])
    except ValueError:
        print(f"⚠️ Skipping {simFile} (couldn’t parse timestep)")
        continue

    # find corresponding CSV row
    row = df[df["Timestep"] == timestep]
    if row.empty:
        print(f"⚠️ No CSV row found for timestep {timestep}, skipping.")
        continue
    row = row.iloc[0]

    # get model parameters
    model = row["Model"]  # should be one of RBETA, CRITBETA, etc.
    if model not in model_settings:
        print(f"⚠️ Model {model} not in model_settings, skipping {simFile}")
        continue
    params = model_settings[model]

    # grab Munit (or MunitUsed if you prefer)
    Munit = row["Munit"]
    # Munit = row["MunitUsed"]

    MunitOffset = row["MunitOffset"]
    MunitSlope = row["MunitSlope"]

    # build output base name
    nameBase = f"/work/vmo703/ipole_outputs/{model}_{timestep}.h5"

    # run the batch
    positronRatios = np.linspace(0, 1, 2)

    makePositronImages(
        simFile,
        Munit,
        MunitOffset,
        MunitSlope,
        positronRatios=positronRatios,
        nameBase=nameBase,
        Rhigh=20,
        **params  # inject electronModel and sigma_transition
    )
    print(f"✅ Finished {model} timestep {timestep}")
