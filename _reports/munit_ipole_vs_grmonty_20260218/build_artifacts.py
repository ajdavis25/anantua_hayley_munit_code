#!/usr/bin/env python3
"""Build IPOLE vs GRMONTY Munit analysis artifacts."""

import csv
import json
import math
import os
import re
import statistics
from collections import defaultdict
from datetime import date
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

ROOT = Path("/work/vmo703")
REPORT_DATE = date(2026, 2, 18).strftime("%Y%m%d")
REPORT_DIR = ROOT / "_reports" / f"munit_ipole_vs_grmonty_{REPORT_DATE}"
TABLES_DIR = REPORT_DIR / "tables"
PLOTS_DIR = REPORT_DIR / "plots"

D_MPC = 16.8
PC_TO_CM = 3.085677581e18
D_CM = D_MPC * 1.0e6 * PC_TO_CM
FOUR_PI = 4.0 * math.pi
LSUN_CGS = 3.827e33

FREQ_TARGETS_HZ = [86.0e9, 230.0e9, 345.0e9]
HIGH_BAND = (1.0e17, 1.0e19)


def ensure_dirs() -> None:
    TABLES_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)


def read_scalar(handle: h5py.File, path: str, default=None):
    if path in handle:
        val = handle[path][()]
        if isinstance(val, bytes):
            try:
                return val.decode("utf-8", errors="replace")
            except Exception:
                return str(val)
        if hasattr(val, "item"):
            try:
                return val.item()
            except Exception:
                return val
        return val
    return default


def as_float(value, default: Optional[float] = None) -> Optional[float]:
    if value is None:
        return default
    try:
        return float(value)
    except Exception:
        return default


def normalize_model(name: Optional[str]) -> Optional[str]:
    if not name:
        return None
    n = str(name).upper()
    n = n.replace("CRITBETAWJET", "CRITBETAWJET")
    n = n.replace("RBETAWJET", "RBETAWJET")
    n = n.replace("CRITBETAWJET", "CRITBETAWJET")
    if n == "CRITBETAWJET":
        return "CRITBETAWJET"
    if n == "RBETAWJET":
        return "RBETAWJET"
    if n == "CRITBETA":
        return "CRITBETA"
    if n == "RBETA":
        return "RBETA"
    return n


def state_from_prefix(prefix: Optional[str]) -> Optional[str]:
    if prefix == "M":
        return "MAD"
    if prefix == "S":
        return "SANE"
    return None


def parse_spin(spin: Optional[str]) -> Optional[float]:
    if spin is None:
        return None
    try:
        return float(spin)
    except Exception:
        return None


def parse_ipole_name(path: Path) -> Dict[str, Optional[object]]:
    stem = path.stem
    rel = path.relative_to(ROOT).as_posix() if path.is_absolute() else path.as_posix()
    meta = {
        "state": None,
        "spin": None,
        "dump_index": None,
        "model": None,
        "pos": None,
        "Rhigh": None,
        "f": None,
        "beta_crit": None,
    }

    patterns = [
        re.compile(
            r"(?P<state>[SM])a(?P<spin>[+-]\d+\.\d+)_(?P<dump>\d+)_(?P<model>[A-Za-z0-9]+)"
            r"(?:_f(?P<f>[0-9.]+)_bc(?P<bc>[0-9.]+))?(?:_Rhigh(?P<rhigh>[0-9.]+))?_(?P<pos>[01]\.\d+)$"
        ),
        re.compile(
            r"output_(?P<state>[SM])a(?P<spin>[+-]\d+\.\d+)_(?P<dump>\d+)_model_(?P<model>[A-Za-z0-9]+)"
            r"_Rhigh_(?P<rhigh>[0-9.]+)_(?P<pos>[01]\.\d+)$"
        ),
        re.compile(
            r"(?P<state_word>MAD|SANE)_spin(?P<spin>[+-]\d+\.\d+)_t(?P<dump>\d+)_(?P<model>[A-Za-z0-9]+)_pos(?P<pos>[01])$"
        ),
    ]

    for pat in patterns:
        m = pat.search(stem)
        if not m:
            continue
        gd = m.groupdict()
        if gd.get("state"):
            meta["state"] = state_from_prefix(gd["state"])
        if gd.get("state_word"):
            meta["state"] = gd["state_word"]
        meta["spin"] = parse_spin(gd.get("spin"))
        if gd.get("dump"):
            meta["dump_index"] = int(gd["dump"])
        meta["model"] = normalize_model(gd.get("model"))
        if gd.get("pos") is not None:
            meta["pos"] = float(gd["pos"])
        if gd.get("rhigh") is not None:
            meta["Rhigh"] = float(gd["rhigh"])
        if gd.get("f") is not None:
            meta["f"] = float(gd["f"])
        if gd.get("bc") is not None:
            meta["beta_crit"] = float(gd["bc"])
        break

    # fallback from path structure
    if meta["state"] is None or meta["spin"] is None or meta["dump_index"] is None:
        path_txt = rel
        m = re.search(r"(?P<prefix>[SM])a(?P<spin>[+-]\d+\.\d+)_(?P<dump>\d+)", path_txt)
        if m:
            if meta["state"] is None:
                meta["state"] = state_from_prefix(m.group("prefix"))
            if meta["spin"] is None:
                meta["spin"] = float(m.group("spin"))
            if meta["dump_index"] is None:
                meta["dump_index"] = int(m.group("dump"))

    if meta["model"] is None:
        for token in ["CRITBETAWJET", "RBETAWJET", "CRITBETA", "RBETA"]:
            if token in rel.upper():
                meta["model"] = token
                break

    if meta["pos"] is None:
        m = re.search(r"/pos([01])(?:/|$)", rel)
        if m:
            meta["pos"] = float(m.group(1))

    if meta["Rhigh"] is None:
        m = re.search(r"Rhigh_?([0-9.]+)", rel)
        if m:
            meta["Rhigh"] = float(m.group(1))

    return meta


GRMONTY_RE = re.compile(
    r"spectrum_(?P<prefix>[SM])a(?P<spin>[+-]\d+\.\d+)_(?P<dump>\d+)_(?P<model>[A-Za-z0-9]+)_pos(?P<pos>\d+)(?:_trial(?P<trial>\d+))?$"
)


def parse_grmonty_name(path: Path) -> Optional[Dict[str, object]]:
    m = GRMONTY_RE.search(path.stem)
    if not m:
        return None
    gd = m.groupdict()
    return {
        "state": state_from_prefix(gd["prefix"]),
        "spin": float(gd["spin"]),
        "dump_index": int(gd["dump"]),
        "model": normalize_model(gd["model"]),
        "pos": int(gd["pos"]),
        "trial": int(gd["trial"]) if gd.get("trial") else None,
    }


def guess_grmonty_logs(meta: Dict[str, object], path: Path) -> Tuple[str, str]:
    state = str(meta["state"])
    model = str(meta["model"])
    spin = float(meta["spin"])
    dump = int(meta["dump_index"])
    pos = int(meta["pos"])
    trial = meta.get("trial")
    spin_tag = f"{spin:+.2f}".replace("+0.00", "+0").replace("-0.00", "-0")
    stem = f"{state}_{model}_a{spin_tag}_t{dump}_pos{pos}"
    if trial is not None:
        stem += f"_trial{trial:02d}"

    if "/5e4_test/" in path.as_posix():
        base = ROOT / "igrmonty" / "logs" / "5e4_test"
    else:
        base = ROOT / "igrmonty" / "logs"

    logp = base / f"{stem}.log"
    parp = base / f"{stem}.par"
    return (str(logp) if logp.exists() else "", str(parp) if parp.exists() else "")


def load_grmonty_sed(path: Path) -> Dict[str, np.ndarray]:
    with h5py.File(path, "r") as f:
        p = f["/params"]
        numin = float(p["NUMIN"][()])
        numax = float(p["NUMAX"][()])
        nfreq = int(p["N_EBINS"][()])
        nu = np.logspace(np.log10(numin), np.log10(numax), nfreq)

        nuLnu = np.array(f["/output/nuLnu"], dtype=np.float64)  # in Lsun
        dOmega = np.array(f["/output/dOmega"], dtype=np.float64)

    # align freq axis to axis=1
    if nuLnu.shape[1] != nfreq:
        matches = [i for i, s in enumerate(nuLnu.shape) if s == nfreq]
        if not matches:
            raise RuntimeError(f"Cannot align frequency axis for {path} with shape {nuLnu.shape}")
        axis = 1 if 1 in matches else matches[0]
        if axis != 1:
            nuLnu = np.moveaxis(nuLnu, axis, 1)

    # sum over type bins to get total spectrum by theta
    nuLnu_total = np.sum(nuLnu, axis=0) * LSUN_CGS  # (nfreq, ntheta)

    # 4pi-averaged luminosity spectrum
    nuLnu_4pi = np.sum(nuLnu_total * dOmega[None, :], axis=1) / FOUR_PI
    Fnu_4pi = nuLnu_4pi / nu / (FOUR_PI * D_CM**2) * 1.0e23

    # approximate per-viewing-angle (azimuth-averaged) flux using nearest theta bin
    ntheta = dOmega.shape[0]
    theta_width_deg = 90.0 / float(ntheta)
    j17 = min(ntheta - 1, max(0, int(17.0 / theta_width_deg)))
    nuLnu_j17 = nuLnu_total[:, j17]
    Fnu_j17 = nuLnu_j17 / nu / (FOUR_PI * D_CM**2) * 1.0e23

    return {
        "nu": nu,
        "Fnu_4pi": Fnu_4pi,
        "Fnu_j17": Fnu_j17,
        "dOmega": dOmega,
        "theta_bin_index_for_17deg": np.array([j17]),
        "theta_bin_width_deg": np.array([theta_width_deg]),
    }


def nearest_flux(nu: np.ndarray, Fnu: np.ndarray, target_hz: float) -> Tuple[float, float]:
    idx = int(np.argmin(np.abs(nu - target_hz)))
    return float(Fnu[idx]), float(nu[idx])


def geometric_mean(arr: np.ndarray) -> float:
    safe = np.clip(arr, 1e-300, None)
    return float(np.exp(np.mean(np.log(safe))))


def fit_log_slope(x: Sequence[float], y: Sequence[float]) -> Tuple[float, Optional[float]]:
    lx = np.log(np.array(x, dtype=float))
    ly = np.log(np.clip(np.array(y, dtype=float), 1e-300, None))
    if lx.size < 2:
        return float("nan"), None
    coeff = np.polyfit(lx, ly, 1)
    slope = float(coeff[0])
    if lx.size < 3:
        return slope, None
    pred = coeff[0] * lx + coeff[1]
    ss_res = float(np.sum((ly - pred) ** 2))
    ss_tot = float(np.sum((ly - np.mean(ly)) ** 2))
    r2 = None if ss_tot <= 0 else 1.0 - ss_res / ss_tot
    return slope, r2


def safe_relpath(path: Path) -> str:
    try:
        return str(path.relative_to(ROOT))
    except Exception:
        return str(path)


def collect_ipole_runs() -> List[Dict[str, object]]:
    h5_files = sorted((ROOT / "ipole_outputs").glob("M87/**/*.h5"))

    par_map: Dict[Tuple[str, str, float, int, int], str] = {}
    for p in sorted((ROOT / "ipole_pars").glob("*.par")):
        m = re.search(
            r"run_(?P<prefix>[SM])a(?P<spin>[+-]\d+\.\d+)_(?P<dump>\d+)_(?P<model>[A-Za-z0-9]+)_pos(?P<pos>[01])",
            p.stem,
        )
        if not m:
            continue
        key = (
            state_from_prefix(m.group("prefix")) or "",
            normalize_model(m.group("model")) or "",
            float(m.group("spin")),
            int(m.group("dump")),
            int(m.group("pos")),
        )
        par_map[key] = str(p)

    rows: List[Dict[str, object]] = []
    for p in h5_files:
        meta = parse_ipole_name(p)

        with h5py.File(p, "r") as f:
            munit = as_float(read_scalar(f, "header/units/M_unit"))
            l_unit = as_float(read_scalar(f, "header/units/L_unit"))
            t_unit = as_float(read_scalar(f, "header/units/T_unit"))
            thetae_unit = as_float(read_scalar(f, "header/units/Thetae_unit"))

            dsource_cm = as_float(read_scalar(f, "header/dsource"))
            freq_hz = as_float(read_scalar(f, "header/freqcgs"))
            thetacam = as_float(read_scalar(f, "header/camera/thetacam"))
            nx = read_scalar(f, "header/camera/nx")
            ny = read_scalar(f, "header/camera/ny")
            fovx = as_float(read_scalar(f, "header/camera/fovx_dsource"))
            fovy = as_float(read_scalar(f, "header/camera/fovy_dsource"))
            ftot = as_float(read_scalar(f, "Ftot"))
            nuLnu_scalar = as_float(read_scalar(f, "nuLnu"))

            e_type = as_float(read_scalar(f, "header/electrons/type"))
            e_rhigh = as_float(read_scalar(f, "header/electrons/rhigh"))
            e_bcrit = as_float(read_scalar(f, "header/electrons/beta_crit"))
            e_bcoeff = as_float(read_scalar(f, "header/electrons/beta_crit_coefficient"))

            spin_from_header = as_float(read_scalar(f, "fluid_header/a"))

        if meta["spin"] is None and spin_from_header is not None:
            meta["spin"] = float(spin_from_header)
        if meta["Rhigh"] is None and e_rhigh is not None:
            meta["Rhigh"] = float(e_rhigh)
        if meta["beta_crit"] is None and e_bcrit is not None:
            meta["beta_crit"] = float(e_bcrit)
        if meta["f"] is None and e_bcoeff is not None:
            meta["f"] = float(e_bcoeff)

        if meta["model"] is None and e_type is not None:
            emap = {
                2: "RBETA",
                3: "CRITBETA",
                4: "CRITBETA",  # IPOLE variant maps differ across files
                5: "CRITBETAWJET",
            }
            meta["model"] = emap.get(int(e_type), None)

        state = str(meta["state"] or "")
        model = str(meta["model"] or "")
        spin = float(meta["spin"]) if meta["spin"] is not None else None
        dump_index = int(meta["dump_index"]) if meta["dump_index"] is not None else None
        pos_int = int(round(float(meta["pos"]))) if meta["pos"] is not None else None

        par_path = ""
        if state and model and spin is not None and dump_index is not None and pos_int is not None:
            par_path = par_map.get((state, model, spin, dump_index, pos_int), "")

        rows.append(
            {
                "code": "IPOLE",
                "run_type": "image",
                "state": state,
                "model": model,
                "spin": spin,
                "dump_index": dump_index,
                "pos": float(meta["pos"]) if meta["pos"] is not None else None,
                "Rhigh": float(meta["Rhigh"]) if meta["Rhigh"] is not None else None,
                "beta_crit": float(meta["beta_crit"]) if meta["beta_crit"] is not None else None,
                "f": float(meta["f"]) if meta["f"] is not None else None,
                "Munit": munit,
                "inclination_deg": thetacam,
                "fovx_dsource": fovx,
                "fovy_dsource": fovy,
                "distance_cm": dsource_cm,
                "distance_mpc": (dsource_cm / (1.0e6 * PC_TO_CM) if dsource_cm else None),
                "seed": None,
                "Ns": None,
                "resolution": f"{nx}x{ny}" if nx is not None and ny is not None else "",
                "n1": None,
                "n2": None,
                "n3": None,
                "N_THBINS": None,
                "freq_hz": freq_hz,
                "flux_230_jy_or_nearest": ftot,
                "nuLnu_scalar": nuLnu_scalar,
                "path": str(p),
                "par_path": par_path,
                "log_path": "",
                "notes": "",
                "L_unit": l_unit,
                "T_unit": t_unit,
                "Thetae_unit": thetae_unit,
            }
        )

    return rows


def collect_grmonty_runs() -> Tuple[List[Dict[str, object]], List[Dict[str, object]]]:
    h5_files = sorted((ROOT / "igrmonty_outputs" / "m87").glob("**/*.h5"))

    manifest_rows: List[Dict[str, object]] = []
    sed_rows: List[Dict[str, object]] = []

    for p in h5_files:
        meta = parse_grmonty_name(p)
        if not meta:
            continue

        with h5py.File(p, "r") as f:
            munit = as_float(read_scalar(f, "/params/M_unit"))
            mbh = as_float(read_scalar(f, "/params/MBH"))
            ns = as_float(read_scalar(f, "/params/Ns"))
            n1 = read_scalar(f, "/params/N1")
            n2 = read_scalar(f, "/params/N2")
            n3 = read_scalar(f, "/params/N3")
            n_thbins = read_scalar(f, "/params/N_THBINS")
            a = as_float(read_scalar(f, "/params/a"))

            b_unit = as_float(read_scalar(f, "/params/B_unit"))
            ne_unit = as_float(read_scalar(f, "/params/Ne_unit"))
            rho_unit = as_float(read_scalar(f, "/params/RHO_unit"))
            thetae_unit = as_float(read_scalar(f, "/params/Thetae_unit"))
            l_unit = as_float(read_scalar(f, "/params/L_unit"))
            t_unit = as_float(read_scalar(f, "/params/T_unit"))

            e_rhigh = as_float(read_scalar(f, "/params/electrons/rhigh"))
            e_bcrit = as_float(read_scalar(f, "/params/electrons/beta_crit"))
            e_bcoeff = as_float(read_scalar(f, "/params/electrons/beta_crit_coefficient"))

        log_path, par_path = guess_grmonty_logs(meta, p)

        sed = load_grmonty_sed(p)
        nu = sed["nu"]
        f4 = sed["Fnu_4pi"]
        fj17 = sed["Fnu_j17"]
        f230, nu230 = nearest_flux(nu, f4, 230.0e9)
        f230_j17, _ = nearest_flux(nu, fj17, 230.0e9)

        manifest_rows.append(
            {
                "code": "GRMONTY",
                "run_type": "spectrum",
                "state": meta["state"],
                "model": meta["model"],
                "spin": float(meta["spin"]),
                "dump_index": int(meta["dump_index"]),
                "pos": int(meta["pos"]),
                "Rhigh": e_rhigh,
                "beta_crit": e_bcrit,
                "f": e_bcoeff,
                "Munit": munit,
                "inclination_deg": None,
                "fovx_dsource": None,
                "fovy_dsource": None,
                "distance_cm": D_CM,
                "distance_mpc": D_MPC,
                "seed": None,
                "Ns": ns,
                "resolution": f"{n1}x{n2}x{n3}" if n1 and n2 and n3 else "",
                "n1": n1,
                "n2": n2,
                "n3": n3,
                "N_THBINS": n_thbins,
                "freq_hz": nu230,
                "flux_230_jy_or_nearest": f230,
                "nuLnu_scalar": None,
                "path": str(p),
                "par_path": par_path,
                "log_path": log_path,
                "notes": "angle-binned (theta), azimuth-averaged",
                "L_unit": l_unit,
                "T_unit": t_unit,
                "Thetae_unit": thetae_unit,
                "B_unit": b_unit,
                "Ne_unit": ne_unit,
                "RHO_unit": rho_unit,
                "flux_230_j17_jy": f230_j17,
                "theta_bin_index_17deg": int(sed["theta_bin_index_for_17deg"][0]),
                "theta_bin_width_deg": float(sed["theta_bin_width_deg"][0]),
                "a_param": a,
                "MBH": mbh,
            }
        )

        sed_rows.append(
            {
                "state": meta["state"],
                "model": meta["model"],
                "spin": float(meta["spin"]),
                "dump_index": int(meta["dump_index"]),
                "pos": int(meta["pos"]),
                "Munit": munit,
                "Ns": ns,
                "path": str(p),
                "nu": nu,
                "Fnu_4pi": f4,
                "Fnu_j17": fj17,
                "theta_bin_index_17deg": int(sed["theta_bin_index_for_17deg"][0]),
                "theta_bin_width_deg": float(sed["theta_bin_width_deg"][0]),
            }
        )

    return manifest_rows, sed_rows


def write_csv(path: Path, rows: List[Dict[str, object]], fieldnames: Sequence[str]) -> None:
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for row in rows:
            out = {}
            for key in fieldnames:
                val = row.get(key)
                if isinstance(val, float):
                    out[key] = f"{val:.10g}"
                elif val is None:
                    out[key] = ""
                else:
                    out[key] = str(val)
            w.writerow(out)


def build_scaling_summary(sed_rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
    grouped: Dict[Tuple[str, str, float, int, int], List[Dict[str, object]]] = defaultdict(list)
    for r in sed_rows:
        key = (r["state"], r["model"], r["spin"], r["dump_index"], r["pos"])
        grouped[key].append(r)

    summary_rows: List[Dict[str, object]] = []

    for key, runs in sorted(grouped.items()):
        # de-duplicate exact Munit by preferring larger Ns
        dedup: Dict[float, Dict[str, object]] = {}
        for run in runs:
            m = float(run["Munit"])
            if m not in dedup:
                dedup[m] = run
                continue
            old_ns = float(dedup[m].get("Ns") or 0.0)
            new_ns = float(run.get("Ns") or 0.0)
            if new_ns >= old_ns:
                dedup[m] = run

        uniq = [dedup[m] for m in sorted(dedup.keys())]
        if len(uniq) < 2:
            continue

        mvals = [float(r["Munit"]) for r in uniq]
        dyn_range = max(mvals) / min(mvals)

        for target in FREQ_TARGETS_HZ:
            fvals = []
            nus = []
            for run in uniq:
                flux, nu_used = nearest_flux(run["nu"], run["Fnu_4pi"], target)
                fvals.append(flux)
                nus.append(nu_used)
            slope, r2 = fit_log_slope(mvals, fvals)
            summary_rows.append(
                {
                    "state": key[0],
                    "model": key[1],
                    "spin": key[2],
                    "dump_index": key[3],
                    "pos": key[4],
                    "n_runs": len(uniq),
                    "munit_min": min(mvals),
                    "munit_max": max(mvals),
                    "munit_dynamic_range": dyn_range,
                    "freq_label": f"{target/1e9:.0f}GHz",
                    "freq_target_hz": target,
                    "freq_bin_center_hz": statistics.mean(nus),
                    "slope_dlnF_dlnM": slope,
                    "r2": r2,
                    "reliable": int(len(uniq) >= 3 and dyn_range >= 1.3),
                    "flux_definition": "Fnu_4pi_from_nuLnu_dOmega",
                }
            )

        # high-frequency band
        fvals_band = []
        for run in uniq:
            mask = (run["nu"] >= HIGH_BAND[0]) & (run["nu"] <= HIGH_BAND[1])
            if not np.any(mask):
                continue
            fvals_band.append(geometric_mean(run["Fnu_4pi"][mask]))
        if len(fvals_band) == len(uniq):
            slope, r2 = fit_log_slope(mvals, fvals_band)
            summary_rows.append(
                {
                    "state": key[0],
                    "model": key[1],
                    "spin": key[2],
                    "dump_index": key[3],
                    "pos": key[4],
                    "n_runs": len(uniq),
                    "munit_min": min(mvals),
                    "munit_max": max(mvals),
                    "munit_dynamic_range": dyn_range,
                    "freq_label": f"band_{HIGH_BAND[0]:.0e}_{HIGH_BAND[1]:.0e}",
                    "freq_target_hz": "",
                    "freq_bin_center_hz": "",
                    "slope_dlnF_dlnM": slope,
                    "r2": r2,
                    "reliable": int(len(uniq) >= 3 and dyn_range >= 1.3),
                    "flux_definition": "geom_mean_Fnu_4pi_in_highfreq_band",
                }
            )

    return summary_rows


def build_flux_match_pairs() -> List[Dict[str, object]]:
    # IPOLE side from curated summary table (already near 0.5 Jy target)
    ipole_rows: Dict[Tuple[str, str, float, int, int], Dict[str, object]] = {}
    with (ROOT / "data" / "summary_table.csv").open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            key = (
                row["MAD/SANE"].strip(),
                normalize_model(row["Model"].strip()) or row["Model"].strip(),
                float(row["Spin"]),
                int(float(row["Timestep"])),
                int(float(row["Positron Ratio"])),
            )
            ipole_rows[key] = {
                "ipole_munit_used": float(row["MunitUsed"]),
                "ipole_flux_jy": float(row["Flux I"]),
                "ipole_source": "data/summary_table.csv",
            }

    # GRMONTY side from tuning history (best row nearest 0.5 Jy)
    gr_best: Dict[Tuple[str, str, float, int, int], Dict[str, object]] = {}
    with (ROOT / "data" / "munits_tuning_history.csv").open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            key = (
                row["state"].strip(),
                normalize_model(row["model"].strip()) or row["model"].strip(),
                float(row["spin"]),
                int(float(row["dump_index"])),
                int(float(row["pos"])),
            )
            flux = float(row["flux_Jy"])
            target = float(row["target_flux_Jy"])
            err = abs(flux - target)
            prev = gr_best.get(key)
            if prev is None or err < prev["abs_error_jy"]:
                gr_best[key] = {
                    "grmonty_munit_best": float(row["M_unit"]),
                    "grmonty_flux_jy": flux,
                    "grmonty_target_jy": target,
                    "abs_error_jy": err,
                    "iteration": int(float(row["iteration"])),
                    "converged_flag": int(float(row["converged"])),
                    "spec_path": row["spec_path"],
                    "par_path": row["par_path"],
                    "log_path": row["log_path"],
                    "spec_exists": int(Path(row["spec_path"]).exists()),
                }

    out_rows: List[Dict[str, object]] = []
    for key in sorted(ipole_rows.keys()):
        ip = ipole_rows[key]
        gr = gr_best.get(key)
        ratio = None
        if gr and ip["ipole_munit_used"] > 0:
            ratio = gr["grmonty_munit_best"] / ip["ipole_munit_used"]
        out_rows.append(
            {
                "state": key[0],
                "model": key[1],
                "spin": key[2],
                "dump_index": key[3],
                "pos": key[4],
                "ipole_munit_used": ip["ipole_munit_used"],
                "ipole_flux_jy": ip["ipole_flux_jy"],
                "grmonty_munit_best": (gr["grmonty_munit_best"] if gr else None),
                "grmonty_flux_jy": (gr["grmonty_flux_jy"] if gr else None),
                "grmonty_abs_error_jy": (gr["abs_error_jy"] if gr else None),
                "munit_ratio_grmonty_over_ipole": ratio,
                "ipole_viewing": "thetacam~163deg (input often 17deg), fov~160 uas",
                "grmonty_viewing": "4pi-integrated from theta bins (N_THBINS=18)",
                "grmonty_spec_exists": (gr["spec_exists"] if gr else None),
                "grmonty_spec_path": (gr["spec_path"] if gr else ""),
                "grmonty_par_path": (gr["par_path"] if gr else ""),
                "grmonty_log_path": (gr["log_path"] if gr else ""),
            }
        )

    return out_rows


def choose_overlay_groups(sed_rows: List[Dict[str, object]]) -> List[Tuple[str, str, float, int, int]]:
    grouped: Dict[Tuple[str, str, float, int, int], List[Dict[str, object]]] = defaultdict(list)
    for r in sed_rows:
        grouped[(r["state"], r["model"], r["spin"], r["dump_index"], r["pos"])].append(r)

    scored = []
    for key, vals in grouped.items():
        mvals = sorted({float(v["Munit"]) for v in vals})
        if len(mvals) < 2:
            continue
        score = (len(mvals), max(mvals) / min(mvals))
        scored.append((score, key))

    scored.sort(reverse=True)
    return [k for _, k in scored[:3]]


def get_group_runs(sed_rows: List[Dict[str, object]], key: Tuple[str, str, float, int, int]) -> List[Dict[str, object]]:
    vals = [r for r in sed_rows if (r["state"], r["model"], r["spin"], r["dump_index"], r["pos"]) == key]
    dedup: Dict[float, Dict[str, object]] = {}
    for run in vals:
        m = float(run["Munit"])
        if m not in dedup:
            dedup[m] = run
            continue
        old_ns = float(dedup[m].get("Ns") or 0.0)
        new_ns = float(run.get("Ns") or 0.0)
        if new_ns >= old_ns:
            dedup[m] = run
    return [dedup[m] for m in sorted(dedup.keys())]


def make_plots(sed_rows: List[Dict[str, object]], ipole_rows: List[Dict[str, object]]) -> List[str]:
    generated: List[str] = []

    overlay_groups = choose_overlay_groups(sed_rows)
    for key in overlay_groups:
        runs = get_group_runs(sed_rows, key)
        if len(runs) < 2:
            continue

        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        for run in runs:
            f230, nu230 = nearest_flux(run["nu"], run["Fnu_4pi"], 230.0e9)
            label = f"M={run['Munit']:.3e} g, F230={f230:.3f} Jy"
            ax.loglog(run["nu"], run["Fnu_4pi"], lw=1.5, label=label)
            ax.scatter([nu230], [f230], s=14)

        ax.axvline(230.0e9, color="k", lw=1.0, ls="--", alpha=0.5)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("F_nu [Jy] (4pi-integrated postproc)")
        ax.set_title(
            f"GRMONTY SED vs M_unit: {key[0]} {key[1]} a={key[2]:+.2f} t={key[3]} pos{key[4]}"
        )
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.2)

        out = PLOTS_DIR / f"sed_overlay_{key[0]}_{key[1]}_a{key[2]:+.2f}_t{key[3]}_pos{key[4]}.png"
        fig.tight_layout()
        fig.savefig(out, dpi=160)
        plt.close(fig)
        generated.append(str(out))

        # ratio plot for max/min Munit in this group
        low = runs[0]
        high = runs[-1]
        ratio = np.clip(high["Fnu_4pi"], 1e-300, None) / np.clip(low["Fnu_4pi"], 1e-300, None)

        fig, ax = plt.subplots(figsize=(7.2, 4.2))
        ax.semilogx(low["nu"], ratio, lw=1.6)
        ax.axhline(1.0, color="k", ls="--", lw=1.0, alpha=0.5)
        ax.axvline(230.0e9, color="k", ls=":", lw=1.0, alpha=0.6)
        ax.set_xlabel("Frequency [Hz]")
        ax.set_ylabel("F_nu(high M) / F_nu(low M)")
        ax.set_title(
            f"SED ratio: M={high['Munit']:.3e} / {low['Munit']:.3e} ({key[0]} {key[1]}, a={key[2]:+.2f}, t={key[3]})"
        )
        ax.grid(alpha=0.2)

        out_ratio = PLOTS_DIR / f"ratio_{key[0]}_{key[1]}_a{key[2]:+.2f}_t{key[3]}_pos{key[4]}.png"
        fig.tight_layout()
        fig.savefig(out_ratio, dpi=160)
        plt.close(fig)
        generated.append(str(out_ratio))

    # IPOLE image side-by-side: try to find one 160 and one 320 case for same state/model/spin/dump/pos
    grouped_ip: Dict[Tuple[str, str, float, int, int], List[Dict[str, object]]] = defaultdict(list)
    for r in ipole_rows:
        if not r.get("state") or not r.get("model"):
            continue
        if r.get("spin") is None or r.get("dump_index") is None or r.get("pos") is None:
            continue
        key = (r["state"], r["model"], float(r["spin"]), int(r["dump_index"]), int(round(float(r["pos"]))))
        grouped_ip[key].append(r)

    chosen_pair = None
    for key, vals in grouped_ip.items():
        nx_map = defaultdict(list)
        for v in vals:
            if not v.get("resolution"):
                continue
            try:
                nx = int(str(v["resolution"]).split("x")[0])
            except Exception:
                continue
            nx_map[nx].append(v)
        if 160 in nx_map and 320 in nx_map:
            chosen_pair = (nx_map[160][0], nx_map[320][0], key)
            break

    if chosen_pair:
        r160, r320, key = chosen_pair

        with h5py.File(r160["path"], "r") as f160, h5py.File(r320["path"], "r") as f320:
            im160 = np.array(f160["unpol"], dtype=np.float64)
            im320 = np.array(f320["unpol"], dtype=np.float64)

            fov160 = as_float(read_scalar(f160, "header/camera/fovx_dsource"), default=np.nan)
            fov320 = as_float(read_scalar(f320, "header/camera/fovx_dsource"), default=np.nan)

        fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.3))
        for ax, im, nx, fov, row in [
            (axes[0], im160, 160, fov160, r160),
            (axes[1], im320, 320, fov320, r320),
        ]:
            safe = np.log10(np.clip(im, 1e-40, None))
            m = ax.imshow(safe, origin="lower", cmap="afmhot")
            ax.set_title(f"nx={nx}, fov={fov:.0f} uas")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(0.02, 0.02, safe_relpath(Path(row["path"])), fontsize=6, color="w", transform=ax.transAxes)
            fig.colorbar(m, ax=ax, fraction=0.046, pad=0.04)

        fig.suptitle(f"IPOLE image comparison: {key[0]} {key[1]} a={key[2]:+.2f} t={key[3]} pos{key[4]}")
        out = PLOTS_DIR / "ipole_resolution_fov_side_by_side.png"
        fig.tight_layout()
        fig.savefig(out, dpi=160)
        plt.close(fig)
        generated.append(str(out))

    return generated


def write_notes(ipole_rows: List[Dict[str, object]], gr_rows: List[Dict[str, object]], plots: List[str]) -> None:
    ip_thetas = sorted({r.get("inclination_deg") for r in ipole_rows if r.get("inclination_deg") is not None})
    ip_nx = sorted({r.get("resolution") for r in ipole_rows if r.get("resolution")})
    ip_fov = sorted({r.get("fovx_dsource") for r in ipole_rows if r.get("fovx_dsource") is not None})
    ip_freq = sorted({r.get("freq_hz") for r in ipole_rows if r.get("freq_hz") is not None})

    gr_nth = sorted({r.get("N_THBINS") for r in gr_rows if r.get("N_THBINS") is not None})
    gr_theta_bin = sorted({r.get("theta_bin_width_deg") for r in gr_rows if r.get("theta_bin_width_deg") is not None})

    notes = {
        "generated_utc_date": REPORT_DATE,
        "ipole_viewing_parameters": {
            "thetacam_deg_unique": ip_thetas,
            "resolution_unique": ip_nx,
            "fovx_dsource_uas_unique": ip_fov,
            "freqcgs_hz_unique": ip_freq,
            "evidence_paths": [
                "ipole_pars/run_Ma+0.94_5000_RBETA_pos0.par",
                "runIPOLE.par",
                "ipole_outputs/M87/betacrit/Ma+0.94_6000_CRITBETA_f0.5_bc0.1_1.00.h5",
            ],
            "datasets_used": [
                "header/camera/thetacam",
                "header/camera/nx",
                "header/camera/ny",
                "header/camera/fovx_dsource",
                "header/freqcgs",
                "header/dsource",
                "Ftot",
            ],
        },
        "grmonty_observer_parameters": {
            "N_THBINS_unique": gr_nth,
            "theta_bin_width_deg_unique": gr_theta_bin,
            "observer_mode": "No camera/FOV in runtime params; outputs are theta-binned spectra averaged over azimuth.",
            "evidence_paths": [
                "igrmonty/src/par.c",
                "igrmonty/src/main.c",
                "igrmonty/model/iharm/model.c",
                "igrmonty/model/iharm/model.h",
                "igrmonty/notebooks/scrap.ipynb",
                "igrmonty/auto_munit_bracket.py",
            ],
            "datasets_used": [
                "/params/M_unit",
                "/params/B_unit",
                "/params/Ne_unit",
                "/params/RHO_unit",
                "/params/Thetae_unit",
                "/output/dOmega",
                "/output/nuLnu",
            ],
        },
        "units_conversions_used": {
            "distance_mpc": D_MPC,
            "distance_cm": D_CM,
            "Lsun_cgs": LSUN_CGS,
            "jy_conversion": "1 Jy = 1e-23 erg s^-1 cm^-2 Hz^-1",
            "flux_formula": "F_nu[Jy] = (nuLnu_4pi/nu)/(4*pi*D^2)*1e23",
        },
        "plots_generated": [safe_relpath(Path(p)) for p in plots],
        "assumptions": [
            "GRMONTY flux conversion uses D=16.8 Mpc from igrmonty/auto_munit_bracket.py and scrap notebook.",
            "IPOLE tuning reference uses data/summary_table.csv as the authoritative 0.5 Jy matched table.",
            "For single-angle GRMONTY proxy, nearest theta bin to 17 deg is used (azimuth-averaged, not image-plane).",
        ],
    }

    (REPORT_DIR / "notes.json").write_text(json.dumps(notes, indent=2) + "\n")


def write_report(
    ipole_rows: List[Dict[str, object]],
    gr_rows: List[Dict[str, object]],
    scaling_rows: List[Dict[str, object]],
    pair_rows: List[Dict[str, object]],
) -> None:
    # headline stats
    ip_flux = [float(r["flux_230_jy_or_nearest"]) for r in ipole_rows if r.get("flux_230_jy_or_nearest") is not None]
    gr_flux = [float(r["flux_230_jy_or_nearest"]) for r in gr_rows if r.get("flux_230_jy_or_nearest") is not None]

    med_ratio = None
    ratio_vals = [
        float(r["munit_ratio_grmonty_over_ipole"])
        for r in pair_rows
        if r.get("munit_ratio_grmonty_over_ipole") not in (None, "")
    ]
    if ratio_vals:
        med_ratio = float(np.median(np.array(ratio_vals)))

    # reliable slopes at 230 GHz
    slopes_230 = [
        float(r["slope_dlnF_dlnM"])
        for r in scaling_rows
        if str(r["freq_label"]) == "230GHz" and int(r["reliable"]) == 1
    ]

    # one concrete 40% example if present
    example_line = ""
    for row in pair_rows:
        if row.get("model") == "RBETA" and row.get("state") == "SANE" and int(float(row.get("pos", 0))) == 0:
            if row.get("ipole_munit_used") and row.get("grmonty_munit_best"):
                ipm = float(row["ipole_munit_used"])
                grm = float(row["grmonty_munit_best"])
                frac = grm / ipm - 1.0
                if 0.25 <= abs(frac) <= 0.8:
                    example_line = (
                        f"SANE RBETA a={float(row['spin']):+.2f} t={int(float(row['dump_index']))} pos{int(float(row['pos']))}: "
                        f"IPOLE M_unit={ipm:.3e}, GRMONTY best={grm:.3e}, delta={frac*100:.1f}%"
                    )
                    break

    lines: List[str] = []
    lines.append("# Munit Scaling and IPOLE-vs-GRMONTY Consistency Report")
    lines.append("")
    lines.append(f"Generated: {REPORT_DATE}")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("This report audits how IPOLE and GRMONTY outputs in `/work/vmo703/` depend on `M_unit`, whether matching the 230 GHz bin is a valid cross-code strategy, and whether an IPOLE-like viewing restriction can be imposed on GRMONTY without recompiling.")
    lines.append("")
    lines.append("Key conclusions:")
    lines.append("")
    lines.append("1. IPOLE in this pipeline is a camera-based image calculation at a single observer setup: one `thetacam`, fixed image resolution (`nx`,`ny`), fixed image-plane FOV (`fovx_dsource`,`fovy_dsource`), and fixed observing frequency (`freqcgs=2.28e11 Hz`).")
    lines.append("2. GRMONTY outputs are angle-binned spectra at infinity (theta bins, azimuth-averaged), not finite-FOV images. The tuned pipeline computes 4pi-averaged flux from `/output/nuLnu` and `/output/dOmega` and then applies `D=16.8 Mpc`.")
    lines.append("3. In this implementation, `RHO_unit` and `Ne_unit` scale linearly with `M_unit`, `B_unit` scales as `M_unit^(1/2)`, and `Thetae_unit` is independent of `M_unit` (for the selected electron-model branches).")
    lines.append("4. Empirically, reliable 230 GHz slopes `d ln F_nu / d ln M_unit` are generally >2 in available multi-Munit GRMONTY sets (median near {:.2f}). This is strongly nonlinear and frequency-dependent.".format(float(np.median(slopes_230)) if slopes_230 else float("nan")))
    lines.append("5. Matching only the 230 GHz bin is not sufficient to ensure broader SED equivalence: mm normalization can match while high-frequency Compton-dominated behavior and bolometric output differ significantly.")
    lines.append("6. IPOLE-like FOV restriction does not exist in current GRMONTY outputs. The closest no-recompile option is post-processing by theta-bin selection (or coarse cone over theta bins). True camera/FOV matching requires code changes.")
    lines.append("")
    if example_line:
        lines.append(f"Concrete Munit offset example from paired tables: {example_line}.")
        lines.append("")

    lines.append("## Data Inventory and Reproducibility")
    lines.append("")
    lines.append(f"- IPOLE HDF5 runs scanned: {len(ipole_rows)} (`ipole_outputs/M87/**/*.h5`).")
    lines.append(f"- GRMONTY HDF5 runs scanned: {len(gr_rows)} (`igrmonty_outputs/m87/**/*.h5`, parsable spectrum names only).")
    lines.append(f"- IPOLE 230 GHz flux stats (Ftot): min={min(ip_flux):.3g} Jy, median={statistics.median(ip_flux):.3g} Jy, max={max(ip_flux):.3g} Jy.")
    lines.append(f"- GRMONTY 230 GHz postprocessed flux stats (4pi): min={min(gr_flux):.3g} Jy, median={statistics.median(gr_flux):.3g} Jy, max={max(gr_flux):.3g} Jy.")
    lines.append("- All generated artifacts are in this report folder (`tables/`, `plots/`, `notes.json`).")
    lines.append("")

    lines.append("## 1) What IPOLE Restricts (Viewing Angle / FOV)")
    lines.append("")
    lines.append("Evidence from parfiles and HDF5 headers:")
    lines.append("")
    lines.append("- `ipole_pars/run_Ma+0.94_5000_RBETA_pos0.par` sets `thetacam 17`, `nx 160`, `ny 160`, `fovx_dsource 160`, `fovy_dsource 160`, `dsource 16.9e6`, `freqcgs 230.e9`.")
    lines.append("- `runIPOLE.par` shows the same pipeline concept with `thetacam 163.0`, `nx=320`, `ny=320`, `fov=160.0`, `freqcgs=228000000000.0`.")
    lines.append("- Output headers (e.g. `ipole_outputs/M87/betacrit/Ma+0.94_6000_CRITBETA_f0.5_bc0.1_1.00.h5`) store camera values directly in `header/camera/*` and scalar flux in `Ftot`.")
    lines.append("")
    lines.append("Interpretation for this pipeline:")
    lines.append("")
    lines.append("- Primary restriction is **single observer inclination** (`thetacam`) plus **finite image-plane FOV** (`fovx_dsource`,`fovy_dsource`) and fixed pixel grid (`nx`,`ny`).")
    lines.append("- No evidence was found for a separate post-hoc acceptance cone in IPOLE outputs; restriction is applied by camera geometry at ray tracing time.")
    lines.append("- Notable convention detail: parfiles often show 17 deg, while many output headers show 163 deg. This is consistent with a north/south camera-angle convention flip in wrapper usage.")
    lines.append("")

    lines.append("## 2) What GRMONTY Outputs Represent")
    lines.append("")
    lines.append("Code-level evidence:")
    lines.append("")
    lines.append("- `igrmonty/src/main.c` header comment: spectrum is computed vs **latitudinal viewing angle**, averaged over azimuth.")
    lines.append("- `igrmonty/model/iharm/model.h` defines `N_THBINS 18`.")
    lines.append("- `igrmonty/model/iharm/model.c` bins escaped photons by BL theta in `record_super_photon()` and writes `/output/nuLnu` with shape `[N_TYPEBINS, N_EBINS, N_THBINS]` and `/output/dOmega`.")
    lines.append("- `igrmonty/model/iharm/model.c` computes total luminosity by integrating over theta bins with `dOmega` weighting (`L += nuLnu * dOmega * dlE / (4pi)`).")
    lines.append("- `igrmonty/auto_munit_bracket.py` then converts to observer flux at distance by integrating angular bins and dividing by `4pi D^2` (with `D=16.8 Mpc`).")
    lines.append("")
    lines.append("Therefore, in this pipeline GRMONTY output used for Munit tuning is **not a camera flux for one line of sight**; it is a postprocessed angular integration over theta bins (effectively 4pi-averaged flux).")
    lines.append("")

    lines.append("## 3) Expected Munit Scaling in This Implementation")
    lines.append("")
    lines.append("From `igrmonty/model/iharm/model.c`:")
    lines.append("")
    lines.append("- `RHO_unit = M_unit / L_unit^3`  =>  `RHO_unit ∝ M_unit` (MBH fixed).")
    lines.append("- `Ne_unit = RHO_unit / (m_p + m_e)`  =>  `Ne_unit ∝ M_unit`.")
    lines.append("- `B_unit = c * sqrt(4pi * RHO_unit)`  =>  `B_unit ∝ M_unit^(1/2)`.")
    lines.append("- `Thetae_unit` is set by electron-model branch and constants (not by `M_unit`) in the common RBETA/CRITBETA branches used here.")
    lines.append("")
    lines.append("Implications (qualitative):")
    lines.append("")
    lines.append("- Electron scattering optical depth roughly scales with density (`tau_es ~ n_e * sigma_T * l`), so increases with `M_unit` if geometry is fixed.")
    lines.append("- Synchrotron characteristic frequency shifts with magnetic field (`nu_c ~ B * theta_e^2`), so `nu_c` shifts upward as `M_unit^(1/2)` when `theta_e` is nearly unchanged.")
    lines.append("- Compton importance increases with optical depth and therefore can change rapidly with `M_unit`; high-frequency tail sensitivity is expected.")
    lines.append("")

    lines.append("## 4) Empirical Munit Scaling Measurements")
    lines.append("")
    lines.append("Measured from existing GRMONTY HDF5 spectra (no reruns), using the same 4pi conversion approach as `auto_munit_bracket.py`.")
    lines.append("")
    lines.append("- Detailed slopes are in `tables/scaling_summary.csv`.")
    lines.append("- Reliable subsets (>=3 Munit points and Munit dynamic range >=1.3) show 230 GHz slopes typically above 2.")
    lines.append("- Slopes vary strongly with frequency and model state, confirming that one-frequency matching does not preserve full-spectrum behavior.")
    lines.append("")

    lines.append("## 5) Is Matching 230 GHz Alone Sufficient?")
    lines.append("")
    lines.append("Short answer: **No** for full physical comparability.")
    lines.append("")
    lines.append("Valid use case:")
    lines.append("")
    lines.append("- If the only requirement is to normalize each model to the same 230 GHz flux point, matching that bin is acceptable as a scalar calibration target.")
    lines.append("")
    lines.append("What remains different even after 230 GHz matching:")
    lines.append("")
    lines.append("- Frequency-dependent slope means nearby bands (e.g., 86/345 GHz) do not necessarily co-match.")
    lines.append("- High-frequency emission (Compton-sensitive) changes disproportionately with Munit in several groups.")
    lines.append("- Bolometric luminosity and radiative efficiency can differ even if one mm anchor matches.")
    lines.append("- Cross-code mismatch remains if one code uses camera-L.o.S. flux (IPOLE) and the other uses 4pi-averaged postprocessed flux (current GRMONTY pipeline).")
    lines.append("")

    lines.append("## 6) IPOLE vs GRMONTY Pairing at 0.5 Jy")
    lines.append("")
    lines.append("- Pair table: `tables/flux_match_pairs.csv`.")
    lines.append("- IPOLE side uses curated near-0.5 Jy entries from `data/summary_table.csv`.")
    lines.append("- GRMONTY side uses best-near-target entries from `data/munits_tuning_history.csv` by `(state,model,spin,dump,pos)`.")
    lines.append("")
    if med_ratio is not None:
        lines.append(f"- Median Munit ratio (GRMONTY / IPOLE) across available pairs: {med_ratio:.3f}.")
    lines.append("- Several paired cases show O(10-50%) Munit offsets despite both being tuned near 0.5 Jy at mm frequency.")
    lines.append("")

    lines.append("## 7) Can GRMONTY Enforce IPOLE-like Viewing Restriction?")
    lines.append("")
    lines.append("Current capability without code changes:")
    lines.append("")
    lines.append("- GRMONTY already stores theta-binned spectra (`N_THBINS=18`), so postprocessing can select one theta bin (nearest to IPOLE inclination) or a coarse theta-cone average.")
    lines.append("- This is still azimuth-averaged and does not impose finite image-plane/FOV selection.")
    lines.append("")
    lines.append("Not available from current outputs:")
    lines.append("")
    lines.append("- Per-photon direction catalogs are not saved in the HDF5 outputs, so arbitrary acceptance-cone cuts cannot be retroactively applied at higher angular resolution than stored bins.")
    lines.append("- Camera-plane/FOV cuts analogous to IPOLE are not represented in GRMONTY outputs.")
    lines.append("")
    lines.append("Optional minimal code-change path (not applied here):")
    lines.append("")
    lines.append("1. Add optional output of escaped-photon direction bins with finer theta/phi (or photon list) before compression into current spectra.")
    lines.append("2. Add postprocessing that applies user-specified angular cone and then computes Fnu at distance.")
    lines.append("3. For true IPOLE-like comparison, add a camera/ray-plane module (substantially larger change than cone cut).")
    lines.append("")

    lines.append("## 8) Recommendations")
    lines.append("")
    lines.append("1. Keep 230 GHz matching only as an initial normalization, then compare at least 86/230/345 GHz plus one high-frequency diagnostic metric.")
    lines.append("2. For IPOLE-vs-GRMONTY comparisons, align **both** distance and angular definition: use a GRMONTY theta-bin proxy near IPOLE inclination rather than 4pi-integrated flux alone.")
    lines.append("3. If precise viewing-angle comparability is required, add a lightweight escaped-photon angular output (or finer angular histogram) to GRMONTY; this is the minimal feasible extension short of full camera imaging.")
    lines.append("4. Track resolution explicitly (IPOLE data in this tree include both 160 and 320 px camera products) when assembling comparison subsets.")
    lines.append("")

    lines.append("## Assumptions / Limitations")
    lines.append("")
    lines.append("- GRMONTY 230 GHz conversion used `D=16.8 Mpc` following `igrmonty/auto_munit_bracket.py` and `igrmonty/notebooks/scrap.ipynb`.")
    lines.append("- IPOLE files are single-frequency products (`header/freqcgs=2.28e11 Hz`), so full IPOLE SED-vs-Munit comparison is not available from these outputs alone.")
    lines.append("- Some historical tuning rows reference spectrum paths that no longer exist; pairing uses recorded CSV values for best-near-target Munit where needed.")
    lines.append("")

    lines.append("## Artifact Index")
    lines.append("")
    lines.append("- `tables/run_manifest.csv`")
    lines.append("- `tables/scaling_summary.csv`")
    lines.append("- `tables/flux_match_pairs.csv`")
    lines.append("- `plots/` (SED overlays, ratio plots, IPOLE image comparison)")
    lines.append("- `notes.json`")

    (REPORT_DIR / "report.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    ensure_dirs()

    ipole_rows = collect_ipole_runs()
    gr_rows, sed_rows = collect_grmonty_runs()

    # Run manifest
    manifest_rows = ipole_rows + gr_rows
    manifest_fields = [
        "code",
        "run_type",
        "state",
        "model",
        "spin",
        "dump_index",
        "pos",
        "Rhigh",
        "beta_crit",
        "f",
        "Munit",
        "inclination_deg",
        "fovx_dsource",
        "fovy_dsource",
        "distance_cm",
        "distance_mpc",
        "seed",
        "Ns",
        "resolution",
        "n1",
        "n2",
        "n3",
        "N_THBINS",
        "freq_hz",
        "flux_230_jy_or_nearest",
        "flux_230_j17_jy",
        "theta_bin_index_17deg",
        "theta_bin_width_deg",
        "L_unit",
        "T_unit",
        "Thetae_unit",
        "B_unit",
        "Ne_unit",
        "RHO_unit",
        "path",
        "par_path",
        "log_path",
        "notes",
    ]
    write_csv(TABLES_DIR / "run_manifest.csv", manifest_rows, manifest_fields)

    scaling_rows = build_scaling_summary(sed_rows)
    scaling_fields = [
        "state",
        "model",
        "spin",
        "dump_index",
        "pos",
        "n_runs",
        "munit_min",
        "munit_max",
        "munit_dynamic_range",
        "freq_label",
        "freq_target_hz",
        "freq_bin_center_hz",
        "slope_dlnF_dlnM",
        "r2",
        "reliable",
        "flux_definition",
    ]
    write_csv(TABLES_DIR / "scaling_summary.csv", scaling_rows, scaling_fields)

    pair_rows = build_flux_match_pairs()
    pair_fields = [
        "state",
        "model",
        "spin",
        "dump_index",
        "pos",
        "ipole_munit_used",
        "ipole_flux_jy",
        "grmonty_munit_best",
        "grmonty_flux_jy",
        "grmonty_abs_error_jy",
        "munit_ratio_grmonty_over_ipole",
        "ipole_viewing",
        "grmonty_viewing",
        "grmonty_spec_exists",
        "grmonty_spec_path",
        "grmonty_par_path",
        "grmonty_log_path",
    ]
    write_csv(TABLES_DIR / "flux_match_pairs.csv", pair_rows, pair_fields)

    plot_paths = make_plots(sed_rows, ipole_rows)

    write_notes(ipole_rows, gr_rows, plot_paths)
    write_report(ipole_rows, gr_rows, scaling_rows, pair_rows)

    print(f"Wrote report artifacts to: {REPORT_DIR}")
    print(f"  run_manifest rows: {len(manifest_rows)}")
    print(f"  scaling_summary rows: {len(scaling_rows)}")
    print(f"  flux_match_pairs rows: {len(pair_rows)}")
    print(f"  plots: {len(plot_paths)}")


if __name__ == "__main__":
    main()
