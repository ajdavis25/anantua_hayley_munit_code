#!/usr/bin/env python3
import argparse
import csv
import re
from pathlib import Path

from runner import IPoleRunner

# electron models we care about
E_MODELS = ["RBETA", "RBETAwJET", "CRITBETA", "CRITBETAwJET"]
POSITIONS = [0, 1]

offset_slope_table = {
    # (electronModel, pos) : (offset, slope)
    ("RBETA", 0): (0.0, 0.0),
    ("RBETA", 1): (0.0, 0.0),
    ("RBETAWJET", 0): (0.0, 0.0),
    ("RBETAWJET", 1): (0.0, 0.0),
    ("CRITBETA", 0): (0.0, 0.0),
    ("CRITBETA", 1): (0.0, 0.0),
    # add more as needed
}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", "--folder", required=True,
        help="Folder with GRMHD .h5 dumps (flat, like Sa-0.5_4000.h5)"
    )
    parser.add_argument(
        "-o", "--outfile", required=True,
        help="Path to output CSV (munits_table style)"
    )
    parser.add_argument(
        "--ipole_path", required=True,
        help="Path to ipole binary or wrapper script (e.g. run_ipole.sh)"
    )
    parser.add_argument("--dump_range", nargs=2, type=int, metavar=("START","END"))

    # NEW: parameters for the optimization
    parser.add_argument(
        "--thetacam", type=float, default=163.0,
        help="Camera inclination angle in degrees (default: 163.0)"
    )
    parser.add_argument(
        "--flux_goal", type=float, default=0.5,
        help="Target flux in Jy (default: 0.5)"
    )

    args = parser.parse_args()

    folder = Path(args.folder)
    runner = IPoleRunner(ipole_path=args.ipole_path)

    thetacam = args.thetacam
    flux_goal = args.flux_goal

    # open CSV for writing
    with open(args.outfile, "w", newline="") as fout:
        writer = csv.writer(fout)
        # header matches your munits_table.csv
        writer.writerow([
            "timestep","flow","model","spin","pos",
            "MunitOffset","MunitSlope","Munit","MunitMax","MunitUsed","notes"
        ])

        for dump_file in sorted(folder.glob("*.h5")):
            # extract timestep and flow from filename
            stem = dump_file.stem  # e.g. Sa-0.5_4000 or Ma+0.94_6000
            flow, timestep = stem.split("_")
            timestep = int(timestep)

            # only keep in range if requested
            if args.dump_range:
                start, end = args.dump_range
                if not (start <= timestep <= end):
                    continue

            # parse spin robustly (handles both + and -)
            match = re.match(r"^(Sa|Ma)([+-]\d+\.?\d*)$", flow)
            if match:
                sim_type, spin_str = match.groups()
                spin = float(spin_str)
            else:
                sim_type, spin = flow, 0.0

            flow_label = {"Sa": "SANE", "Ma": "MAD"}.get(sim_type, sim_type)

            # loop over models and pos
            for electronModel in E_MODELS:
                for pos in POSITIONS:
                    Munit_val = runner.find_munit(
                        dump_file, thetacam=thetacam, flux_goal=flux_goal
                    )

                    offset, slope = offset_slope_table.get((electronModel, pos), (0.0, 0.0))
                    positronRatio = float(pos)  # adjust if pos maps differently
                    MunitUsed = offset + (slope * Munit_val) / (1 + 2 * positronRatio)

                    writer.writerow([
                        timestep,
                        flow_label,
                        electronModel,
                        spin,
                        pos,
                        offset,
                        slope,
                        Munit_val,
                        0.0,  # placeholder for MunitMax
                        MunitUsed,
                        ""
                    ])

if __name__ == "__main__":
    main()
