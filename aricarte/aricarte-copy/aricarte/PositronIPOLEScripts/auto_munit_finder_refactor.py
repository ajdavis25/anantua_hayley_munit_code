#!/usr/bin/env python3
import sys, io, os, re, argparse, subprocess, math
import numpy as np
import pandas as pd

# slurm-safe stdout
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, line_buffering=True)

# constants
TARGET_FLUX = 0.5
REL_TOL = 0.05

DATA_CSV = "/work/vmo703/data/paper_data.csv"
OUTPUT_CSV = "/work/vmo703/data/paper_output.csv"

REQUIRED_COLUMNS = [
    "dump_index", "timestep", "state", "model",
    "spin", "Rhigh", "Munit",
]

ELECTRON_MAP = {
    "RBETA": (2, 2.0),
    "RBETAwJET": (2, 0.5),
    "CRITBETA": (4, 2.0),
    "CRITBETAwJET": (4, 0.5),
}


# validation
def validate_input_csv(df):
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        sys.exit(f"[error] missing required columns: {missing}")

    if not set(df["state"]).issubset({"MAD", "SANE"}):
        sys.exit("[error] invalid values in state column")

    if not set(df["model"]).issubset(ELECTRON_MAP.keys()):
        sys.exit("[error] invalid model name detected")

    if df[REQUIRED_COLUMNS].isnull().any().any():
        sys.exit("[error] NaNs detected in required columns")


# helpers
def run_ipole_and_parse_flux(cmd):
    result = subprocess.run(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
        universal_newlines=True
    )
    text = result.stdout + "\n" + result.stderr

    m = re.findall(r'Ftot:\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*Jy', text)
    if m:
        return float(m[-1])

    m = re.findall(r'I,Q,U,V\s*\[Jy\]\s*:\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)', text)
    if m:
        return float(m[-1])

    return None


def bracket_and_bisect(simFile, nameBase, fpos, init_M, Rhigh,
                       ipole, inclination, freq_Hz,
                       electronModel, sigma_transition, sigma_cut,
                       max_iter=20):

    def run(M):
        cmd = [
            "python", "-c",
            (
                "from ipole_many_models import runIPOLE; "
                f"runIPOLE('{simFile}', '{nameBase}_{fpos:.2f}.h5', {M}, "
                f"ipoleExecutable='{ipole}', thetacam={inclination}, "
                f"Rhigh={Rhigh}, freq_Hz={freq_Hz}, "
                f"positronRatio={fpos}, electronModel={electronModel}, "
                f"sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"
            )
        ]
        return run_ipole_and_parse_flux(cmd)

    M, F = init_M, run(init_M)
    if F is None:
        return None, None

    if abs(F - TARGET_FLUX) / TARGET_FLUX < REL_TOL:
        return M, F

    lo_M = hi_M = M
    lo_F = hi_F = F

    for _ in range(10):
        if F < TARGET_FLUX:
            hi_M *= 3
            hi_F = run(hi_M)
        else:
            lo_M /= 3
            lo_F = run(lo_M)

        if lo_F is not None and hi_F is not None:
            if (lo_F - TARGET_FLUX) * (hi_F - TARGET_FLUX) <= 0:
                break

    if lo_F is None or hi_F is None:
        return None, None

    for _ in range(max_iter):
        mid_M = math.sqrt(lo_M * hi_M)
        mid_F = run(mid_M)
        if mid_F is None:
            return None, None

        if abs(mid_F - TARGET_FLUX) / TARGET_FLUX < REL_TOL:
            return mid_M, mid_F

        if mid_F < TARGET_FLUX:
            lo_M, lo_F = mid_M, mid_F
        else:
            hi_M, hi_F = mid_M, mid_F

    return mid_M, mid_F


def already_completed(row, out_df):
    if out_df is None or out_df.empty:
        return False

    mask = (
        (out_df["dump_index"] == row["dump_index"]) &
        (out_df["timestep"] == row["timestep"]) &
        (out_df["state"] == row["state"]) &
        (out_df["model"] == row["model"]) &
        (out_df["spin"] == row["spin"]) &
        (out_df["Rhigh"] == row["Rhigh"]) &
        (out_df["converged"] == True)
    )
    return mask.any()


def print_summary(input_df, output_df):
    total = len(input_df)
    completed = 0 if output_df is None else output_df["converged"].sum()
    print(f"\nsummary: {completed}/{total} completed\n")

    if output_df is None:
        print("all rows pending.")
        return

    merged = input_df.merge(
        output_df,
        on=["dump_index", "timestep", "state", "model", "spin", "Rhigh"],
        how="left"
    )

    missing = merged[merged["converged"] != True]
    print("missing configurations:")
    print(missing.groupby(["state", "model", "spin", "Rhigh"]).size())


# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated Munit finder for IPOLE runs")
    parser.add_argument("--row", type=int, default=None,
                        help="Row index to process (overridden by SLURM_ARRAY_TASK_ID if set)")
    parser.add_argument("--dry-run", action="store_true",
                        help="Print resolved configuration and exit without running IPOLE")
    parser.add_argument("--summary", action="store_true",
                        help="Print progress / missing configurations and exit")
    args = parser.parse_args()

    # load input CSV
    df = pd.read_csv(DATA_CSV).rename(columns=str.strip)
    validate_input_csv(df)

    # load output CSV if it exists
    out_df = pd.read_csv(OUTPUT_CSV) if os.path.exists(OUTPUT_CSV) else None

    # summary-only mode
    if args.summary:
        print_summary(df, out_df)
        sys.exit(0)

    # determine row index: SLURM overrides CLI
    if "SLURM_ARRAY_TASK_ID" in os.environ:
        row_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    elif args.row is not None:
        row_id = args.row
    else:
        sys.exit("[error] no row specified (use --row or SLURM_ARRAY_TASK_ID)")

    if row_id < 0 or row_id >= len(df):
        sys.exit(f"[error] row index {row_id} out of bounds (0–{len(df)-1})")

    row = df.iloc[row_id]

    # skip completed rows
    if already_completed(row, out_df):
        print(f"[skip] row {row_id} already converged", flush=True)
        sys.exit(0)

    # dry-run mode
    if args.dry_run:
        print(f"[dry-run] row {row_id}")
        print(row)
        sys.exit(0)

    # begin actual computation
    state = row["state"]
    model = row["model"]
    spin = row["spin"]
    timestep = row["timestep"]
    dump_index = row["dump_index"]
    Rhigh = row["Rhigh"]
    Munit = row["Munit"]
    A0 = 0.0
    B0 = 1.0

    sim_prefix = "Ma" if state == "MAD" else "Sa"
    simFile = f"/work/vmo703/grmhd_dump_samples/{sim_prefix}{spin:+}_{timestep}.h5"

    electronModel, sigma_cut = ELECTRON_MAP[model]

    nameBase = (
        f"/work/vmo703/ipole_outputs/"
        f"{sim_prefix}{spin:+}_{timestep}_{model}_Rhigh{Rhigh}"
    )

    denom0, denom1 = 1.0, 3.0
    init0 = A0 + B0 * Munit / denom0
    init1 = A0 + B0 * Munit / denom1

    ipole_exec = "/work/vmo703/aricarte/run_ipole.sh"

    M0, F0 = bracket_and_bisect(
        simFile, nameBase, 0.0, init0, Rhigh,
        ipole_exec, 163.0, 228e9,
        electronModel, 0.5, sigma_cut
    )

    M1, F1 = bracket_and_bisect(
        simFile, nameBase, 1.0, init1, Rhigh,
        ipole_exec, 163.0, 228e9,
        electronModel, 0.5, sigma_cut
    )

    if M0 is not None and M1 is not None:
        B = (M0 - M1) / (Munit * (1 - 1/3))
        A = M0 - B * Munit
        converged = True
    else:
        A, B = np.nan, np.nan
        converged = False

    out = {
        "row_id": row_id,
        "dump_index": dump_index,
        "timestep": timestep,
        "state": state,
        "model": model,
        "spin": spin,
        "Rhigh": Rhigh,
        "Munit": Munit,
        "MunitOffset": A,
        "MunitSlope": B,
        "MunitUsed_pos0": M0,
        "MunitUsed_pos1": M1,
        "converged": converged,
        "notes": row.get("notes", "")
    }

    pd.DataFrame([out]).to_csv(
        OUTPUT_CSV,
        mode="a",
        header=not os.path.exists(OUTPUT_CSV),
        index=False
    )
