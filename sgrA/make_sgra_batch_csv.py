#!/usr/bin/env python3
"""
make_sgra_batch_csv.py
----------------------
Scans /work/vmo703/grmhd_dump_samples/ and generates a CSV listing
all available timesteps for MAD/SANE and both spin signs, expanded
over the four electron-heating models.
"""

import re, pandas as pd
from pathlib import Path

samples_dir = Path("/work/vmo703/grmhd_dump_samples")
out_csv = Path("/work/vmo703/data/sgra_batch.csv")

# regex pattern: e.g., Ma+0.94_5000.h5
pattern = re.compile(r"(Ma|Sa)([+-]\d+\.\d+)_(\d+)\.h5")

rows = []
for f in samples_dir.glob("*.h5"):
    m = pattern.match(f.name)
    if not m:
        continue
    evolution = "MAD" if m.group(1) == "Ma" else "SANE"
    spin = m.group(2)
    timestep = int(m.group(3))
    for model in ["RBETA", "CRITBETA", "RBETAwJET", "CRITBETAwJET"]:
        rows.append({
            "Timestep": timestep,
            "MAD/SANE": evolution,
            "Model": model,
            "Spin": spin,
            "Rhigh": 20
        })

df = pd.DataFrame(rows).sort_values(["MAD/SANE", "Spin", "Timestep", "Model"])
out_csv.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_csv, index=False)

print(f"[done] wrote {len(df)} rows → {out_csv}")
print(df.head(10))
