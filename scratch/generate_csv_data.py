#!/usr/bin/env python3

import pandas as pd
from itertools import product

# configuration space
dump_indices = [4000, 5000, 6000]
timesteps = {4000: 20000, 5000: 25000, 6000: 30000}

states = ["SANE", "MAD"]
models = ["RBETA", "RBETAwJET", "CRITBETA", "CRITBETAwJET"]
spins = [-0.5, 0.94]
rhighs = [40, 80, 160]
positron_states = [0, 1]

# build rows
rows = []
job_id = 0

for dump_index, state, model, spin, Rhigh, pos in product(
    dump_indices, states, models, spins, rhighs, positron_states
):
    rows.append({
        "job_id": job_id,
        "dump_index": dump_index,
        "timestep": timesteps[dump_index],
        "state": state,
        "model": model,
        "spin": spin,
        "Rhigh": Rhigh,
        "pos": pos,
        "notes": ""
    })
    job_id += 1

# create DataFrame and write CSV
df = pd.DataFrame(rows, columns=[
    "job_id",
    "dump_index",
    "timestep",
    "state",
    "model",
    "spin",
    "Rhigh",
    "pos",
    "notes"
])

output_file = "paper_data.csv"
df.to_csv(output_file, index=False)

print(f"wrote {len(df)} rows to {output_file}")
