#!/usr/bin/env python3
import argparse, csv, re, h5py
from pathlib import Path

from munit_oskey import offset_slope_key

# electron models (order must match 4-tuple lists)
E_MODELS = ["RBETA", "RBETAwJET", "CRITBETA", "CRITBETAwJET"]
MODEL_IDX = {m: i for i, m in enumerate(E_MODELS)}
POSITIONS = [0, 1]

# baseline Munit constants
SANE_BASE_MUNIT = 1.83e27
MAD_BASE_MUNIT = 7.49e24


def get_sim_key(flow_label: str, spin: float) -> str:
    """exact keys expected by offset_slope_key: 'SANE-0.5', 'SANE+0.94', etc"""
    if abs(spin - 0.94) < 1e-3:
        return f"{flow_label}+0.94"
    if abs(spin + 0.94) < 1e-3:
        return f"{flow_label}-0.94"
    return f"{flow_label}{spin:+.1f}" # e.g. SANE-0.5 / MAD-0.5


def read_timestep_from_h5(dump_file: Path) -> int:
    """read true simulation time [M] from HDF5 header"""
    candidate_keys = ["/fluid_header/t", "/header/t", "t"]
    with h5py.File(dump_file, "r") as H:
        for key in candidate_keys:
            if key in H:
                return int(round(float(H[key][()])))
    raise KeyError(f"no timestep field found in {dump_file}")


def base_munit_for_flow(flow_label: str) -> float:
    return SANE_BASE_MUNIT if flow_label == "SANE" else MAD_BASE_MUNIT


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", required=True,
                        help="Folder with GRMHD .h5 dumps (flat, like Sa-0.5_4000.h5)")
    parser.add_argument("-o", "--outfile", required=True,
                        help="Path to output CSV (munits_table style)")
    parser.add_argument("--dump_range", nargs=2, type=int, metavar=("START", "END"))
    args = parser.parse_args()

    folder = Path(args.folder)

    rows = []
    seen = set() # dedupe (dump_index, flow, spin, model, pos)

    for dump_file in sorted(folder.glob("*.h5")):
        stem = dump_file.stem # e.g. Sa-0.5_4000 or Sa+0.94_5000
        try:
            flow_token, dump_index_str = stem.split("_")
            dump_index = int(dump_index_str)
        except Exception:
            print(f"[warn] unexpected filename format: {stem}")
            continue

        # parse flow & spin from filename
        m = re.match(r"^(Sa|Ma)([+-]\d+\.?\d*)$", flow_token)
        if not m:
            print(f"[warn] could not parse flow/spin from {stem}")
            continue
        sim_type, spin_str = m.groups()
        spin = float(spin_str)
        flow_label = "SANE" if sim_type == "Sa" else "MAD"

        # true timestep
        try:
            timestep = read_timestep_from_h5(dump_file)
        except KeyError as e:
            print(f"[warn] {e}")
            continue

        if args.dump_range:
            start, end = args.dump_range
            if not (start <= timestep <= end):
                print(f"[skip] outside dump_range: dump_index={dump_index}, timestep={timestep}")
                continue

        tkey = str(timestep)
        if tkey not in offset_slope_key:
            print(f"[skip] no offset/slope data for this timestep={timestep} (dump_index={dump_index})")
            continue

        sim_key = get_sim_key(flow_label, spin)
        params = offset_slope_key[tkey]

        if sim_key not in params["MunitOffset"]:
            print(f"[skip] sim_key {sim_key} missing at timestep={timestep}")
            continue

        offsets_list = params["MunitOffset"][sim_key]
        slopes_list = params["MunitSlope"][sim_key]
        base_M = base_munit_for_flow(flow_label)

        for model in E_MODELS:
            mi = MODEL_IDX[model]
            if mi >= len(offsets_list) or mi >= len(slopes_list):
                print(f"[warn] incomplete data for {sim_key} at timestep={timestep}, model={model}")
                continue

            offset = float(offsets_list[mi])
            slope = float(slopes_list[mi])

            for pos in POSITIONS:
                key = (dump_index, flow_label, round(spin, 2), model, pos)
                if key in seen:
                    continue
                seen.add(key)

                MunitUsed = offset + (slope * base_M) / (1 + 2 * pos)
                rows.append((dump_index, timestep, flow_label, model, spin, pos,
                             offset, slope, base_M, MunitUsed, ""))
                print(f"[info] wrote row: dump={dump_index}, t={timestep}, sim={sim_key}, model={model}, pos={pos}")

    # stable, human-friendly ordering
    flow_order = {"SANE": 0, "MAD": 1}

    def spin_rank(x: float) -> int:
        if abs(x - (-0.5)) < 1e-6: return 0
        if abs(x - (0.94)) < 1e-3: return 1
        return 2

    rows.sort(key=lambda r: (
        flow_order.get(r[2], 99),   # flow
        spin_rank(r[4]),            # spin order
        r[0],                       # dump_index
        MODEL_IDX[r[3]],            # model order
        r[5]                        # pos
    ))

    # write
    with open(args.outfile, "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(["dump_index", "timestep", "MAD/SANE", "model", "spin", "pos",
                         "MunitOffset", "MunitSlope", "Munit", "MunitUsed", "notes"])
        writer.writerows(rows)

    print(f"[done] wrote {len(rows)} rows to {args.outfile}")


if __name__ == "__main__":
    main()
