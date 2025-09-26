#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
munit_finder.py
---------------
scans ~/grmhd_dump_samples for .h5 dumps,
parses metadata from filenames, expands over models + positron fractions,
computes MunitUsed, and write a CSV matching munits_table.csv format
"""

import re, pandas as pd
from pathlib import Path


DUMPS_DIR = Path("/work/vmo703/grmhd_dump_samples")
OUTPUT_CSV = Path("/work/vmo703/data/munits.csv")
MODELS = ["RBETA", "RBETAWJET", "CRITBETA", "CRITBETAWJET"]
POSITRON_FRACS = [0, 1]


def parse_filename(filename: str):
    """
    parse metadata from filename
    parse filename like 'Ma-0.5_4000.h5' -> (MAD/SANE, spin, timestep)
    """
    m = re.match(r"([MS])a([+-]?\d+(?:\.\d+)?)_(\d+)\.h5", filename)
    if not m:
        return None
    mad_sane = "MAD" if m.group(1) == "M" else "SANE"
    spin = float(m.group(2))
    timestep = int(m.group(3))
    return mad_sane, spin, timestep


def get_base_munit(mad_sane: str) -> float:
    """
    return base Munit for MAD/SANE
    """
    if mad_sane == "MAD":
        return 7.49e24
    elif mad_sane == "SANE":
        return 1.83e27
    else:
        raise ValueError(f"Unknown MAD/SANE: {mad_sane}")


def get_offsets_and_slopes(model: str):
    """
    return (MunitOffset, MunitSlope) for given model
    """
    if model in ["RBETA", "RBETAWJET"]:
        return 1.3e25 if model == "RBETA" else 8.0e24, 2.8 if model == "RBETA" else 2.9
    elif model in ["CRITBETA", "CRITBETAWJET"]:
        return 1.0e25 if model == "CRITBETA" else 6.0e24, 2.5 if model == "CRITBETA" else 2.6
    else:
        raise ValueError(f"Unknown model: {model}")
    

def compute_munit_used(munit_offset: float, munit_slope: float, base_munit: float, positron_frac: int):
    """
    compute MunitUsed = (offset + slope * munit) / (1 + 2 * positron_fraction)
    """
    return (munit_offset + munit_slope * base_munit) / (1 + 2 * positron_frac)


def main():
    rows = []

    for f in DUMPS_DIR.glob("*.h5"):
        parsed = parse_filename(f.name)
        if not parsed:
            print(f"skipping {f.name}: couldn't parse")
            continue

        mad_sane, spin, timestep = parsed

        for model in MODELS:
            offset, slope = get_offsets_and_slopes(model, mad_sane, spin)

            for positron_frac in POSITRON_FRACS:
                munit = get_base_munit(mad_sane)
                munit_used = compute_munit_used(offset, slope, munit, positron_frac)
                rows.append({
                    "Timestep": timestep,
                    "MAD/SANE": mad_sane,
                    "Model": model,
                    "Spin": spin,
                    "positron frac": positron_frac,
                    "MunitOffset": offset,
                    "MunitSlope": slope,
                    "Munit": munit,
                    "MunitUsed": munit_used,
                    "M87 GRMONTY spectrum?": "idk"
                }) 

    df = pd.DataFrame(rows)

    note = "finding singular munits for grmonty from 'ipole_positron_unlocked.py' equation using munits from fall 23 document"
    with open(OUTPUT_CSV, "w") as f:
        f.write(note + "\n")
    df.to_csv(OUTPUT_CSV, mode="a", index=False)

    print(f"wrote {len(df)} rows to {OUTPUT_CSV}")

if __name__ == "__main__":
    main()
