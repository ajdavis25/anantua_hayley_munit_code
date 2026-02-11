#!/usr/bin/env python3
import pandas as pd

INPUT_CSV  = "paper_data.csv"
OUTPUT_CSV = "paper_output.csv"


def main():
    df = pd.read_csv(INPUT_CSV)

    df = df.rename(columns={"MAD/SANE": "state"})

    group_keys = [
        "dump_index",
        "timestep",
        "state",
        "model",
        "spin",
        "Rhigh",
    ]

    rows = []
    row_id = 0

    for keys, g in df.groupby(group_keys):
        g0 = g[g["pos"] == 0]
        g1 = g[g["pos"] == 1]

        if len(g0) != 1 or len(g1) != 1:
            print(f"WARNING: skipping incomplete group {keys}")
            continue

        g0 = g0.iloc[0]
        g1 = g1.iloc[0]

        rows.append({
            "row_id": row_id,
            "dump_index": g0.dump_index,
            "timestep": g0.timestep,
            "state": g0.state,
            "model": g0.model,
            "spin": g0.spin,
            "Rhigh": g0.Rhigh,
            "Munit": g0.Munit,
            "MunitOffset": g0.MunitOffset,
            "MunitSlope": g0.MunitSlope,
            "MunitUsed_pos0": g0.MunitUsed,
            "MunitUsed_pos1": g1.MunitUsed,
            "converged": True,
            "notes": g0.notes if "notes" in g0 else "",
        })

        row_id += 1

    out = pd.DataFrame(rows)

    out = out[
        [
            "row_id",
            "dump_index",
            "timestep",
            "state",
            "model",
            "spin",
            "Rhigh",
            "Munit",
            "MunitOffset",
            "MunitSlope",
            "MunitUsed_pos0",
            "MunitUsed_pos1",
            "converged",
            "notes",
        ]
    ]

    out.to_csv(OUTPUT_CSV, index=False)
    print(f"wrote {len(out)} rows -> {OUTPUT_CSV}")


if __name__ == "__main__":
    main()
