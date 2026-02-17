#!/usr/bin/env python3
import sys, io, os, re, argparse, subprocess, math
import numpy as np
import pandas as pd
from ipole_many_models import runIPOLE

# ---------------------------------------------------------------------
# SLURM-safe stdout
# ---------------------------------------------------------------------
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, line_buffering=True)

# ---------------------------------------------------------------------
# constants
# ---------------------------------------------------------------------
TARGET_FLUX = 0.5
REL_TOL = 0.05

DATA_CSV = "/work/vmo703/data/paper_data_critbeta_both.csv"
OUTPUT_CSV = "/work/vmo703/data/paper_output_critbeta.csv"
PARFILE_PATH = "./runIPOLE.par"
DEFAULT_UNPOL = False

REQUIRED_COLUMNS = [
    "dump_index", "timestep", "state", "model",
    "spin", "Munit"
]

MODEL_MAP = {
    "RBETA": (2, 2.0),
    "RBETAWJET": (2, 0.5),
    "CRITBETA": (4, 2.0),
    "CRITBETAWJET": (4, 0.5),
}
RHIGH_MODELS = {"RBETA", "RBETAWJET"}
CRIT_MODELS = {"CRITBETA", "CRITBETAWJET"}


def _canonical_column_name(col):
    raw = col.strip()
    norm = re.sub(r"[^a-z0-9]+", "", raw.lower())
    mapping = {
        "jobid": "job_id",
        "dumpindex": "dump_index",
        "timestep": "timestep",
        "state": "state",
        "model": "model",
        "spin": "spin",
        "rhigh": "Rhigh",
        "rlow": "Rlow",
        "f": "f",
        "betacrit": "beta_crit",
        "pos": "pos",
        "munit": "Munit",
        "notes": "notes",
    }
    return mapping.get(norm, raw)


def normalize_input_columns(df):
    renamed = {_c: _canonical_column_name(_c) for _c in df.columns}
    out = df.rename(columns=renamed)
    dups = [c for c in out.columns if list(out.columns).count(c) > 1]
    if dups:
        raise ValueError(f"ambiguous columns after normalization: {sorted(set(dups))}")
    return out


def _value_missing(v):
    if pd.isna(v):
        return True
    if isinstance(v, str) and not v.strip():
        return True
    return False


def _series_missing(df, col):
    if col not in df.columns:
        return pd.Series(True, index=df.index)
    s = df[col]
    return s.isna() | s.astype(str).str.strip().eq("")


def _required_float(v, field, row_id, model):
    if _value_missing(v):
        raise ValueError(f"row {row_id} model={model} missing required '{field}'")
    try:
        return float(v)
    except (TypeError, ValueError):
        raise ValueError(f"row {row_id} model={model} invalid '{field}'={v!r}")


def sanity_check_model_knobs(df):
    model_upper = df["model"].astype(str).str.upper()
    miss_rhigh_mask = model_upper.isin(RHIGH_MODELS) & _series_missing(df, "Rhigh")
    miss_crit_mask = model_upper.isin(CRIT_MODELS) & (
        _series_missing(df, "f") | _series_missing(df, "beta_crit")
    )
    missing_rhigh = int(miss_rhigh_mask.sum())
    missing_crit = int(miss_crit_mask.sum())

    print(f"[sanity] RBETA/RBETAwJET rows missing Rhigh: {missing_rhigh}", flush=True)
    print(f"[sanity] CRITBETA/CRITBETAwJET rows missing f or beta_crit: {missing_crit}", flush=True)

    if missing_rhigh or missing_crit:
        if missing_rhigh:
            bad_idx = int(df.index[miss_rhigh_mask][0])
        else:
            bad_idx = int(df.index[miss_crit_mask][0])
        bad_model = str(df.iloc[bad_idx]["model"]).upper()
        raise ValueError(
            f"row {bad_idx} model={bad_model} missing required model-specific parameter(s) "
            f"(missing_rhigh={missing_rhigh}, missing_crit={missing_crit})"
        )

# ---------------------------------------------------------------------
# validation + idempotency helpers
# ---------------------------------------------------------------------
def validate_input_csv(df):
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(f"missing required columns: {missing}")

    if not set(df["state"]).issubset({"MAD", "SANE"}):
        raise ValueError("invalid values in state column")

    model_upper = df["model"].astype(str).str.upper()
    if not set(model_upper).issubset(set(MODEL_MAP)):
        raise ValueError("invalid model name detected")

    if df[REQUIRED_COLUMNS].isnull().any().any():
        raise ValueError("NaNs detected in required columns")


def already_completed(row_id, row, out_df):
    if out_df is None or out_df.empty:
        return False

    if "row_id" in out_df.columns:
        return ((out_df["row_id"] == row_id) & (out_df["converged"] == True)).any()

    rhigh = row.get("Rhigh", np.nan)
    mask = (
        (out_df["dump_index"] == row["dump_index"]) &
        (out_df["timestep"] == row["timestep"]) &
        (out_df["state"] == row["state"]) &
        (out_df["model"] == row["model"]) &
        (out_df["spin"] == row["spin"]) &
        (out_df["Rhigh"] == rhigh) &
        (out_df["converged"] == True)
    )
    return mask.any()


def resolve_model_params(row, row_id):
    model = str(row["model"]).upper()
    if model in RHIGH_MODELS:
        rhigh = _required_float(row.get("Rhigh"), "Rhigh", row_id, model)
        rlow = None
        if "Rlow" in row and not _value_missing(row.get("Rlow")):
            rlow = _required_float(row.get("Rlow"), "Rlow", row_id, model)
        return {"Rhigh": rhigh, "Rlow": rlow, "f": None, "beta_crit": None}

    if model in CRIT_MODELS:
        f_coeff = _required_float(row.get("f"), "f", row_id, model)
        beta_crit = _required_float(row.get("beta_crit"), "beta_crit", row_id, model)
        return {"Rhigh": None, "Rlow": None, "f": f_coeff, "beta_crit": beta_crit}

    raise ValueError(f"row {row_id} model={model} is unsupported")


def _fmt_token(v):
    if v is None:
        return "None"
    return f"{float(v):g}"


def _py_literal(v):
    if v is None:
        return "None"
    return repr(float(v))


def build_run_tag(sim_prefix, spin, dump_index, model, model_params):
    if model in CRIT_MODELS:
        return (
            f"{sim_prefix}{spin:+}_{dump_index}_{model}_"
            f"f{_fmt_token(model_params['f'])}_bc{_fmt_token(model_params['beta_crit'])}"
        )

    tag = f"{sim_prefix}{spin:+}_{dump_index}_{model}_Rhigh{_fmt_token(model_params['Rhigh'])}"
    if model_params["Rlow"] is not None:
        tag += f"_Rlow{_fmt_token(model_params['Rlow'])}"
    return tag

# ---------------------------------------------------------------------
# ipole helpers (unchanged physics)
# ---------------------------------------------------------------------
def _log_effective_params(row_id, model, run_tag, simFile, nameBase, MunitUsed, fpos, inclination,
                          model_params, freq_Hz, fov, npixel, counterjet, rmax_geo,
                          electronModel, sigma_transition, sigma_cut):
    pol_state = "off" if DEFAULT_UNPOL else "on"
    print(
        "[ipole] "
        f"row={row_id} run_tag={run_tag} model={model} sim={simFile} out={nameBase}_{fpos:.2f}.h5 "
        f"Munit={MunitUsed:.6e} incl={inclination} mode={electronModel} "
        f"Rhigh={_fmt_token(model_params['Rhigh'])} Rlow={_fmt_token(model_params['Rlow'])} "
        f"f={_fmt_token(model_params['f'])} beta_crit={_fmt_token(model_params['beta_crit'])} "
        f"sigma_transition={sigma_transition} "
        f"sigma_cut={sigma_cut} positronRatio={fpos} counterjet={counterjet} pol={pol_state} "
        f"parfile={PARFILE_PATH} freq_Hz={freq_Hz} fov={fov} npixel={npixel} rmax_geo={rmax_geo}",
        flush=True
    )


def _run_ipole_once(simFile, nameBase, MunitUsed, fpos, ipole, inclination,
                    freq_Hz, fov, npixel, counterjet, rmax_geo,
                    electronModel, sigma_transition, sigma_cut, row_id, model, run_tag,
                    model_params):

    _log_effective_params(
        row_id, model, run_tag, simFile, nameBase, MunitUsed, fpos, inclination, model_params,
        freq_Hz, fov, npixel, counterjet, rmax_geo, electronModel, sigma_transition, sigma_cut
    )

    result = subprocess.run(
        [
            "python", "-c",
            (
                "from ipole_many_models import runIPOLE; "
                f"runIPOLE('{simFile}', '{nameBase}_{fpos:.2f}.h5', {MunitUsed}, "
                f"ipoleExecutable='{ipole}', thetacam={inclination}, "
                f"Rhigh={_py_literal(model_params['Rhigh'])}, "
                f"Rlow={_py_literal(model_params['Rlow'])}, "
                f"freq_Hz={freq_Hz}, "
                f"fov={fov}, npixel={npixel}, counterjet={counterjet}, rmax_geo={rmax_geo}, "
                f"positronRatio={fpos}, electronModel={electronModel}, "
                f"beta_crit_coefficient={_py_literal(model_params['f'])}, "
                f"beta_crit={_py_literal(model_params['beta_crit'])}, "
                f"sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"
            )
        ],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        universal_newlines=True,
    )

    text = result.stdout + "\n" + result.stderr

    m = re.findall(r"Ftot:\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)\s*Jy", text)
    if m:
        return float(m[-1])

    m = re.findall(r"I,Q,U,V\s*\[Jy\]\s*:\s*([+\-]?\d+(?:\.\d+)?(?:[eE][+\-]?\d+)?)", text)
    if m:
        return float(m[-1])

    return None


def bracket_and_bisect(simFile, nameBase, fpos, init_M,
                       ipole, inclination, freq_Hz,
                       electronModel, sigma_transition, sigma_cut,
                       row_id, model, run_tag, model_params, max_iter=20):

    def run(M):
        return _run_ipole_once(
            simFile, nameBase, M, fpos, ipole,
            inclination, freq_Hz, 160.0, 320,
            0, 50, electronModel, sigma_transition, sigma_cut, row_id, model,
            run_tag, model_params
        )

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

# ---------------------------------------------------------------------
# main
# ---------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Automated Munit finder")
    parser.add_argument("--row", type=int, default=None)
    args = parser.parse_args()

    # load input
    df = normalize_input_columns(pd.read_csv(DATA_CSV).rename(columns=str.strip))
    validate_input_csv(df)
    df["model"] = df["model"].astype(str).str.upper()
    sanity_check_model_knobs(df)

    # load output if exists
    out_df = pd.read_csv(OUTPUT_CSV) if os.path.exists(OUTPUT_CSV) else None

    # resolve row index
    if "SLURM_ARRAY_TASK_ID" in os.environ:
        row_id = int(os.environ["SLURM_ARRAY_TASK_ID"])
    elif args.row is not None:
        row_id = args.row
    else:
        sys.exit("[error] no row specified")

    row = df.iloc[row_id]
    model_params = resolve_model_params(row, row_id)

    # skip completed
    if already_completed(row_id, row, out_df):
        print(f"[skip] row {row_id} already converged", flush=True)
        sys.exit(0)

    state = row["state"]
    model = row["model"]
    spin = row["spin"]
    timestep = row["timestep"]
    dump_index = row["dump_index"]
    Munit = row["Munit"]

    electronModel, sigma_cut = MODEL_MAP[model]

    sim_prefix = "Ma" if state == "MAD" else "Sa"
    simFile = f"/work/vmo703/grmhd_dump_samples/{sim_prefix}{spin:+}_{dump_index}.h5"
    run_tag = build_run_tag(sim_prefix, spin, dump_index, model, model_params)

    print(
        f"processing row {row_id}\n"
        f"  dump_index={dump_index}, "
        f"timestep={timestep}, "
        f"state={state}, "
        f"model={model}, "
        f"spin={spin}, "
        f"run_tag={run_tag}, "
        f"Rhigh={_fmt_token(model_params['Rhigh'])}, "
        f"Rlow={_fmt_token(model_params['Rlow'])}, "
        f"f={_fmt_token(model_params['f'])}, "
        f"beta_crit={_fmt_token(model_params['beta_crit'])}",
        flush=True
    )

    nameBase = (
        f"/work/vmo703/ipole_outputs/"
        f"{run_tag}"
    )

    A0, B0 = 0.0, 1.0
    init0 = A0 + B0 * Munit
    init1 = A0 + B0 * Munit / 3

    ipole_exec = "/work/vmo703/aricarte/run_ipole.sh"

    M0, F0 = bracket_and_bisect(
        simFile, nameBase, 0.0, init0,
        ipole_exec, 163.0, 228e9,
        electronModel, 0.5, sigma_cut, row_id, model, run_tag, model_params
    )

    M1, F1 = bracket_and_bisect(
        simFile, nameBase, 1.0, init1,
        ipole_exec, 163.0, 228e9,
        electronModel, 0.5, sigma_cut, row_id, model, run_tag, model_params
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
        "run_tag": run_tag,
        "spin": spin,
        "Rhigh": row.get("Rhigh", np.nan),
        "Rhigh_used": model_params["Rhigh"],
        "Rlow_used": model_params["Rlow"],
        "f_used": model_params["f"],
        "beta_crit_used": model_params["beta_crit"],
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

    print(f"[done] row {row_id} converged={converged}", flush=True)
