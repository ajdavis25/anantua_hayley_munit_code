#!/usr/bin/env python3
"""QA pipeline for finished GRMONTY outputs in /work/vmo703/igrmonty_outputs/m87.

Outputs written under:
  - /work/vmo703/igrmonty_outputs/m87/_qa/grmonty_output_inventory.csv
  - /work/vmo703/igrmonty_outputs/m87/_qa/grmonty_h5_sanity_report.csv
  - /work/vmo703/igrmonty_outputs/m87/_qa/grmonty_pairwise_report.csv
  - /work/vmo703/igrmonty_outputs/m87/_qa/plots/
  - /work/vmo703/igrmonty_outputs/m87/_qa/QA_REPORT.md
"""

import csv
import json
import math
import re
import argparse
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402  (after Agg backend)


BASE_DIR = Path("/work/vmo703/igrmonty_outputs/m87")
QA_DIR = BASE_DIR / "_qa"
PLOTS_DIR = QA_DIR / "plots"
PAIR_PLOTS_DIR = PLOTS_DIR / "pairs"

INVENTORY_CSV = QA_DIR / "grmonty_output_inventory.csv"
SANITY_CSV = QA_DIR / "grmonty_h5_sanity_report.csv"
PAIRWISE_CSV = QA_DIR / "grmonty_pairwise_report.csv"
REPORT_MD = QA_DIR / "QA_REPORT.md"


def safe_slug(text: str) -> str:
    return re.sub(r"[^A-Za-z0-9._+-]+", "_", text).strip("_")


def as_scalar(value: Any) -> Any:
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8", "ignore")
    if isinstance(value, np.ndarray):
        if value.shape == ():
            return as_scalar(value.item())
        return value
    if isinstance(value, np.generic):
        return value.item()
    return value


def float_or_none(value):
    try:
        if value is None:
            return None
        out = float(value)
        if math.isnan(out) or math.isinf(out):
            return None
        return out
    except Exception:
        return None


def int_or_none(value):
    try:
        if value is None:
            return None
        return int(value)
    except Exception:
        return None


def parse_filename(name):
    stem = Path(name).stem
    up = name.upper()
    out = {
        "filename": name,
        "run_id": stem.replace("spectrum_", "", 1),
        "state_letter": None,
        "state": "UNKNOWN",
        "spin_from_name": None,
        "dump_time": None,
        "pos": None,
        "heating_model": "UNKNOWN",
        "wjet_flag": "UNKNOWN",
        "wjet_label": None,
        "trial_tag": None,
    }

    m = re.search(r"_([SM])a([+-]?\d+(?:\.\d+)?)_([0-9]+)", name)
    if m:
        out["state_letter"] = m.group(1)
        out["spin_from_name"] = float_or_none(m.group(2))
        out["dump_time"] = int_or_none(m.group(3))

    pos_m = re.search(r"_pos([0-9]+)", name, flags=re.IGNORECASE)
    if pos_m:
        out["pos"] = int_or_none(pos_m.group(1))

    trial_m = re.search(r"(trial[0-9]+|TEST)", name, flags=re.IGNORECASE)
    if trial_m:
        out["trial_tag"] = trial_m.group(1)

    if "CRITBETA" in up:
        out["heating_model"] = "CRITBETA"
    elif "RBETA" in up:
        out["heating_model"] = "RBETA"

    if "WJET" in up or "STRONGJET" in up or "WEAKJET" in up:
        out["wjet_flag"] = "ON"
        if "STRONGJET" in up:
            out["wjet_label"] = "strongJET"
        elif "WEAKJET" in up:
            out["wjet_label"] = "weakJET"
        else:
            out["wjet_label"] = "wJET"
    elif re.search(r"(^|_)JET(_|\.|$)", up):
        out["wjet_flag"] = "ON"
        out["wjet_label"] = "JET"
    else:
        out["wjet_flag"] = "OFF"
        out["wjet_label"] = "none"

    if out["state_letter"] == "S":
        out["state"] = "SANE"
    elif out["state_letter"] == "M":
        out["state"] = "MAD"
    elif "MAD" in up:
        out["state"] = "MAD"
    elif "SANE" in up:
        out["state"] = "SANE"

    return out


def parse_log_filename_config(name):
    # e.g. SANE_RBETAwJET_a-0.5_t4000_pos0_trial01.par
    m = re.match(
        r"^(SANE|MAD)_([A-Z]+(?:wJET)?)_a([+-]?\d+(?:\.\d+)?)_t(\d+)_pos(\d+)(?:_trial\d+)?\.(?:par|log)$",
        name,
    )
    if not m:
        return None

    state_name, model_name, spin_s, dump_s, pos_s = m.groups()
    if "CRITBETA" in model_name:
        heating = "CRITBETA"
    elif "RBETA" in model_name:
        heating = "RBETA"
    else:
        heating = "UNKNOWN"

    wjet = "ON" if "wJET" in model_name else "OFF"
    state_letter = state_name[0]  # S or M
    spin = float_or_none(spin_s)
    dump_time = int_or_none(dump_s)
    pos = int_or_none(pos_s)
    if None in (state_letter, spin, dump_time, heating, wjet, pos):
        return None
    return (state_letter, spin, dump_time, heating, wjet, pos)


def parse_h5_config_strict(name):
    # Strict canonical naming: spectrum_Sa-0.5_6000_RBETA_pos0.h5
    m = re.match(
        r"^spectrum_([SM])a([+-]?\d+(?:\.\d+)?)_(\d+)_(CRITBETA|RBETA)(wJET)?_pos(\d+)\.h5$",
        name,
    )
    if not m:
        return None

    state_letter, spin_s, dump_s, heating, wjet_tag, pos_s = m.groups()
    spin = float_or_none(spin_s)
    dump_time = int_or_none(dump_s)
    pos = int_or_none(pos_s)
    wjet = "ON" if wjet_tag else "OFF"
    if None in (state_letter, spin, dump_time, heating, wjet, pos):
        return None
    return (state_letter, spin, dump_time, heating, wjet, pos)


def config_from_meta(meta):
    state_letter = meta.get("state_letter")
    spin = float_or_none(meta.get("spin_from_name"))
    dump_time = int_or_none(meta.get("dump_time"))
    heating = meta.get("heating_model")
    wjet = meta.get("wjet_flag")
    pos = int_or_none(meta.get("pos"))
    if None in (state_letter, spin, dump_time, heating, wjet, pos):
        return None
    return (state_letter, spin, dump_time, heating, wjet, pos)


def get_group_scalar(group: h5py.Group, key: str) -> Any:
    try:
        return as_scalar(group[key][()])
    except Exception:
        return None


def find_scalar_by_regex(group: h5py.Group, pattern: str) -> Any:
    regex = re.compile(pattern, flags=re.IGNORECASE)
    for key in group.keys():
        if regex.search(key):
            value = get_group_scalar(group, key)
            if not isinstance(value, np.ndarray):
                return value
    return None


def fmt_spin_for_key(spin: Any) -> str:
    val = float_or_none(spin)
    if val is None:
        return "NA"
    return f"{val:+.4f}"


def append_note(rec, note):
    if not note:
        return
    rec.setdefault("_notes", [])
    if note not in rec["_notes"]:
        rec["_notes"].append(note)


def extract_spectrum_from_nulnu(nu_lnu):
    """Return 1D integrated Stokes-I spectrum and optional polarization stats."""
    meta = {
        "pol_max_abs_frac": None,
        "pol_nonfinite_count": 0,
    }

    if nu_lnu.ndim == 3:
        # Typical shape (8, n_nu, n_theta): channel 0 is Stokes I.
        stokes_i = np.asarray(nu_lnu[0, :, :], dtype=float)
        spectrum = np.nansum(stokes_i, axis=-1)
        return spectrum, meta

    if nu_lnu.ndim == 2:
        # Could be (n_nu, n_theta) or (channels, n_nu); integrate over last axis.
        spectrum = np.nansum(np.asarray(nu_lnu, dtype=float), axis=-1)
        return spectrum, meta

    if nu_lnu.ndim == 1:
        return np.asarray(nu_lnu, dtype=float), meta

    return None, meta


def compute_f230(nu_hz, spectrum):
    target = 230e9
    valid = np.isfinite(nu_hz) & np.isfinite(spectrum) & (nu_hz > 0) & (spectrum > 0)
    if int(np.sum(valid)) < 2:
        return None

    x = np.log10(nu_hz[valid])
    y = np.log10(spectrum[valid])
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    x0 = math.log10(target)
    if x0 < float(x[0]) or x0 > float(x[-1]):
        return None
    return float(10.0 ** np.interp(x0, x, y))


def compute_roughness(nu_hz, spectrum):
    band = (nu_hz >= 1.0e10) & (nu_hz <= 1.0e15)
    valid = band & np.isfinite(spectrum) & (spectrum > 0)
    if int(np.sum(valid)) < 10:
        valid = np.isfinite(spectrum) & (spectrum > 0)
        if int(np.sum(valid)) < 10:
            return None
    logy = np.log10(spectrum[valid])
    d = np.diff(logy)
    if d.size == 0:
        return None
    return float(np.sqrt(np.mean(d**2)))


def compare_spectra(nu_a, y_a, nu_b, y_b):
    # Use A-grid; interpolate B if needed.
    nu = np.asarray(nu_a, dtype=float)
    a = np.asarray(y_a, dtype=float)
    b = np.asarray(y_b, dtype=float)

    if a.shape != b.shape or not np.allclose(nu_a, nu_b, rtol=1e-10, atol=0.0):
        valid_b = np.isfinite(nu_b) & np.isfinite(b) & (nu_b > 0)
        if int(np.sum(valid_b)) < 2:
            return {
                "median_ratio": None,
                "p05_ratio": None,
                "p95_ratio": None,
                "max_abs_log10_ratio": None,
                "wide_jump_frac": None,
                "valid_bins": 0,
                "suspicious_pair": True,
                "notes": "insufficient bins for interpolation",
            }
        xb = np.log10(nu_b[valid_b])
        yb = b[valid_b]
        order = np.argsort(xb)
        xb = xb[order]
        yb = yb[order]
        x = np.log10(nu)
        b = np.interp(x, xb, yb, left=np.nan, right=np.nan)

    finite_ab = np.isfinite(a) & np.isfinite(b) & (a > 0) & (b > 0)
    if not np.any(finite_ab):
        return {
            "median_ratio": None,
            "p05_ratio": None,
            "p95_ratio": None,
            "max_abs_log10_ratio": None,
            "wide_jump_frac": None,
            "valid_bins": 0,
            "suspicious_pair": True,
            "notes": "no positive overlap bins",
        }

    peak = float(np.nanmax(np.concatenate([a[finite_ab], b[finite_ab]])))
    floor = peak * 1.0e-12
    valid = finite_ab & (a > floor) & (b > floor)
    n_valid = int(np.sum(valid))
    if n_valid < 20:
        return {
            "median_ratio": None,
            "p05_ratio": None,
            "p95_ratio": None,
            "max_abs_log10_ratio": None,
            "wide_jump_frac": None,
            "valid_bins": n_valid,
            "suspicious_pair": True,
            "notes": "insufficient overlap bins",
        }

    ratio = b[valid] / a[valid]
    logr = np.log10(ratio)
    med = float(np.median(ratio))
    p05 = float(np.percentile(ratio, 5))
    p95 = float(np.percentile(ratio, 95))
    max_dex = float(np.max(np.abs(logr)))
    wide_jump_frac = float(np.mean(np.abs(logr) > 2.0))  # factor >100 changes

    suspicious = False
    notes = []
    if (med > 1e4 or med < 1e-4) and wide_jump_frac > 0.6:
        suspicious = True
        notes.append("broad >4 dex median shift in high-power band")
    if max_dex > 4.0 and wide_jump_frac > 0.7:
        suspicious = True
        notes.append("many bins differ by >4 dex in high-power band")

    return {
        "median_ratio": med,
        "p05_ratio": p05,
        "p95_ratio": p95,
        "max_abs_log10_ratio": max_dex,
        "wide_jump_frac": wide_jump_frac,
        "valid_bins": n_valid,
        "suspicious_pair": suspicious,
        "notes": "; ".join(notes),
    }


def choose_representative(recs):
    # Prefer passing sanity + non-trial outputs.
    return sorted(
        recs,
        key=lambda r: (
            not bool(r.get("pass_sanity", False)),
            bool(r.get("trial_tag")),
            r.get("filename", ""),
        ),
    )[0]


def write_csv(path, rows, fieldnames):
    with path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            out_row = {k: row.get(k, "") for k in fieldnames}
            writer.writerow(out_row)


def plot_quicklook(rec, nu_hz, spectrum, out_path):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    pos = np.isfinite(nu_hz) & np.isfinite(spectrum) & (nu_hz > 0) & (spectrum > 0)
    if np.any(pos):
        ax.loglog(nu_hz[pos], spectrum[pos], lw=1.2)
    else:
        idx = np.arange(spectrum.size)
        ax.plot(idx, spectrum, lw=1.0)
        ax.set_xlabel("bin index")

    ax.set_title(
        f"{rec.get('filename')} | {rec.get('state')} {rec.get('heating_model')} "
        f"WJET={rec.get('wjet_flag')} pass={rec.get('pass_sanity')}"
    )
    ax.set_xlabel("nu [Hz]")
    ax.set_ylabel("nuLnu (integrated, arb.)")
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def plot_pair(
    title,
    label_a,
    nu_a,
    y_a,
    label_b,
    nu_b,
    y_b,
    out_path,
):
    fig, ax = plt.subplots(figsize=(7.2, 4.4))

    pos_a = np.isfinite(nu_a) & np.isfinite(y_a) & (nu_a > 0) & (y_a > 0)
    pos_b = np.isfinite(nu_b) & np.isfinite(y_b) & (nu_b > 0) & (y_b > 0)

    if np.any(pos_a):
        ax.loglog(nu_a[pos_a], y_a[pos_a], lw=1.2, label=label_a)
    if np.any(pos_b):
        ax.loglog(nu_b[pos_b], y_b[pos_b], lw=1.2, label=label_b)

    ax.set_title(title)
    ax.set_xlabel("nu [Hz]")
    ax.set_ylabel("nuLnu (integrated, arb.)")
    ax.grid(True, alpha=0.25)
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(out_path, dpi=140)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="QA GRMONTY outputs")
    parser.add_argument(
        "--crossref-logs",
        action="store_true",
        help="Filter analyzed .h5 files to configs present in log/par filenames.",
    )
    parser.add_argument(
        "--logs-dir",
        default="/work/vmo703/igrmonty/logs",
        help="Directory containing .par/.log files for cross-reference.",
    )
    parser.add_argument(
        "--canonical-only",
        action="store_true",
        help="When cross-referencing, keep only canonical output names without trial/TEST suffixes.",
    )
    args = parser.parse_args()

    QA_DIR.mkdir(parents=True, exist_ok=True)
    PLOTS_DIR.mkdir(parents=True, exist_ok=True)
    PAIR_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

    all_h5_files = sorted([p for p in BASE_DIR.rglob("*.h5") if "_qa" not in str(p)])
    h5_files = list(all_h5_files)
    filtered_out_noncanonical = []
    filtered_out_not_in_logs = []
    log_cfg_set = set()
    selected_cfg_set = set()
    missing_log_cfgs = []
    crossref_enabled = bool(args.crossref_logs)

    if crossref_enabled:
        logs_dir = Path(args.logs_dir)
        log_paths = sorted(
            [p for p in logs_dir.glob("*") if p.suffix.lower() in (".par", ".log")]
        )
        for p in log_paths:
            cfg = parse_log_filename_config(p.name)
            if cfg is not None:
                log_cfg_set.add(cfg)

        selected = []
        for p in all_h5_files:
            strict_cfg = parse_h5_config_strict(p.name)
            if args.canonical_only and strict_cfg is None:
                filtered_out_noncanonical.append(p)
                continue

            cfg = strict_cfg
            if cfg is None:
                cfg = config_from_meta(parse_filename(p.name))

            if cfg is None or cfg not in log_cfg_set:
                filtered_out_not_in_logs.append(p)
                continue
            selected.append(p)
            if strict_cfg is not None:
                selected_cfg_set.add(strict_cfg)
        h5_files = sorted(selected)
        missing_log_cfgs = sorted(log_cfg_set - selected_cfg_set)

    print(f"[qa] scanning {BASE_DIR}")
    print(f"[qa] found {len(all_h5_files)} .h5 files total")
    if crossref_enabled:
        print(
            f"[qa] crossref mode: selected {len(h5_files)} files from logs ({len(log_cfg_set)} log configs)"
        )
        print(
            "[qa] filtered out "
            f"{len(filtered_out_noncanonical)} non-canonical + {len(filtered_out_not_in_logs)} non-log-matched files"
        )
        print(f"[qa] log configs missing canonical output file: {len(missing_log_cfgs)}")
    else:
        print(f"[qa] selected {len(h5_files)} .h5 files (no crossref filter)")

    if h5_files:
        sample = [p.name for p in h5_files[:8]]
        print("[qa] sample files:")
        for name in sample:
            print(f"  - {name}")

    inventory_rows = []
    sanity_rows = []
    spectra = {}

    for path in h5_files:
        meta = parse_filename(path.name)
        rec = {
            "filepath": str(path),
            "filename": path.name,
            "run_id": meta["run_id"],
            "state_letter": meta["state_letter"],
            "state": meta["state"],
            "spin": meta["spin_from_name"],
            "spin_source": "filename" if meta["spin_from_name"] is not None else None,
            "dump_time": meta["dump_time"],
            "pos": meta["pos"],
            "heating_model": meta["heating_model"],
            "wjet_flag": meta["wjet_flag"],
            "wjet_label": meta["wjet_label"],
            "trial_tag": meta["trial_tag"],
            "m_unit": None,
            "mbh": None,
            "tp_over_te": None,
            "bias": None,
            "rmax": None,
            "rhigh": None,
            "beta_crit": None,
            "param_keys_count": None,
            "header_keys_count": None,
            "top_level_keys": "",
            "parsed_params": "",
            "open_ok": False,
            "has_spectrum": False,
            "nan_count": 0,
            "inf_count": 0,
            "lnu_min": None,
            "lnu_max": None,
            "nu_min_hz": None,
            "nu_max_hz": None,
            "freq_negative_count": 0,
            "freq_nonmonotonic": False,
            "spectrum_bins": 0,
            "spectrum_finite_count": 0,
            "spectrum_positive_count": 0,
            "spectrum_negative_count": 0,
            "spectrum_neg_frac": None,
            "spectrum_all_zero": False,
            "spectrum_min": None,
            "spectrum_max": None,
            "spectrum_max_abs": None,
            "roughness_logdiff_rms": None,
            "f230ghz_nuLnu": None,
            "pol_products_present": False,
            "pol_max_abs_frac": None,
            "pol_nonfinite_count": 0,
            "hard_fail": False,
            "suspicious": False,
            "pass_sanity": False,
            "notes": "",
            "_notes": [],
        }

        parsed_params = {
            "state": rec["state"],
            "heating_model": rec["heating_model"],
            "wjet_flag": rec["wjet_flag"],
            "spin": rec["spin"],
            "dump_time": rec["dump_time"],
            "trial_tag": rec["trial_tag"],
        }

        try:
            with h5py.File(path, "r") as h5f:
                rec["open_ok"] = True
                keys = list(h5f.keys())
                rec["top_level_keys"] = "|".join(keys)

                params = h5f.get("params")
                header = h5f.get("fluid_header")
                output = h5f.get("output")

                if isinstance(params, h5py.Group):
                    rec["param_keys_count"] = len(params.keys())
                    rec["m_unit"] = float_or_none(get_group_scalar(params, "M_unit"))
                    rec["mbh"] = float_or_none(get_group_scalar(params, "MBH"))
                    rec["tp_over_te"] = float_or_none(get_group_scalar(params, "TP_OVER_TE"))
                    rec["bias"] = float_or_none(get_group_scalar(params, "bias"))
                    rec["rhigh"] = float_or_none(
                        find_scalar_by_regex(params, r"(rhigh|r_high)")
                    )
                    rec["beta_crit"] = float_or_none(
                        find_scalar_by_regex(params, r"(beta_crit|betacrit)")
                    )

                    if rec["spin"] is None:
                        spin_param = float_or_none(get_group_scalar(params, "a"))
                        if spin_param is not None:
                            rec["spin"] = spin_param
                            rec["spin_source"] = "params"

                if isinstance(header, h5py.Group):
                    rec["header_keys_count"] = len(header.keys())
                    rec["rmax"] = float_or_none(get_group_scalar(header, "rmax"))
                    if rec["spin"] is None:
                        spin_hdr = float_or_none(get_group_scalar(header, "a"))
                        if spin_hdr is not None:
                            rec["spin"] = spin_hdr
                            rec["spin_source"] = "fluid_header"

                parsed_params.update(
                    {
                        "m_unit": rec["m_unit"],
                        "mbh": rec["mbh"],
                        "tp_over_te": rec["tp_over_te"],
                        "bias": rec["bias"],
                        "rmax": rec["rmax"],
                        "rhigh": rec["rhigh"],
                        "beta_crit": rec["beta_crit"],
                    }
                )

                if isinstance(output, h5py.Group) and ("lnu" in output) and ("nuLnu" in output):
                    rec["has_spectrum"] = True

                    lnu = np.asarray(output["lnu"][:], dtype=float)
                    nu_lnu = np.asarray(output["nuLnu"][:], dtype=float)
                    spectrum, pol_meta = extract_spectrum_from_nulnu(nu_lnu)
                    rec["pol_max_abs_frac"] = pol_meta["pol_max_abs_frac"]
                    rec["pol_nonfinite_count"] = pol_meta["pol_nonfinite_count"]

                    # Only run polarization-specific checks if explicit polarization datasets exist.
                    pol_ratio_keys = {
                        "p_over_i",
                        "v_over_i",
                        "p_i",
                        "v_i",
                        "pi",
                        "vi",
                    }
                    pol_value_keys = {
                        "q",
                        "u",
                        "v",
                        "stokes_q",
                        "stokes_u",
                        "stokes_v",
                    }
                    for key in output.keys():
                        k = key.lower()
                        if k not in pol_ratio_keys and k not in pol_value_keys:
                            continue
                        rec["pol_products_present"] = True
                        arr = np.asarray(output[key][:], dtype=float)
                        rec["pol_nonfinite_count"] += int(np.sum(~np.isfinite(arr)))
                        if k in pol_ratio_keys:
                            finite = arr[np.isfinite(arr)]
                            if finite.size:
                                kmax = float(np.nanmax(np.abs(finite)))
                                if rec["pol_max_abs_frac"] is None:
                                    rec["pol_max_abs_frac"] = kmax
                                else:
                                    rec["pol_max_abs_frac"] = max(rec["pol_max_abs_frac"], kmax)

                    if spectrum is None:
                        append_note(rec, "could not extract 1D spectrum from nuLnu")
                    else:
                        nu_hz = np.exp(lnu)

                        rec["nan_count"] = int(np.isnan(lnu).sum() + np.isnan(nu_lnu).sum())
                        rec["inf_count"] = int(np.isinf(lnu).sum() + np.isinf(nu_lnu).sum())
                        rec["lnu_min"] = float(np.nanmin(lnu)) if lnu.size else None
                        rec["lnu_max"] = float(np.nanmax(lnu)) if lnu.size else None

                        finite_nu = np.isfinite(nu_hz)
                        if np.any(finite_nu):
                            rec["nu_min_hz"] = float(np.nanmin(nu_hz[finite_nu]))
                            rec["nu_max_hz"] = float(np.nanmax(nu_hz[finite_nu]))

                        rec["freq_negative_count"] = int(np.sum((~finite_nu) | (nu_hz <= 0)))
                        rec["freq_nonmonotonic"] = bool(np.any(np.diff(nu_hz) <= 0))

                        rec["spectrum_bins"] = int(spectrum.size)
                        finite_spec = np.isfinite(spectrum)
                        rec["spectrum_finite_count"] = int(np.sum(finite_spec))
                        rec["spectrum_positive_count"] = int(np.sum(spectrum > 0))
                        rec["spectrum_negative_count"] = int(np.sum(spectrum < 0))
                        rec["spectrum_neg_frac"] = float(
                            rec["spectrum_negative_count"] / max(rec["spectrum_bins"], 1)
                        )
                        rec["spectrum_all_zero"] = bool(
                            np.any(finite_spec) and np.allclose(spectrum[finite_spec], 0.0)
                        )

                        if np.any(finite_spec):
                            rec["spectrum_min"] = float(np.nanmin(spectrum[finite_spec]))
                            rec["spectrum_max"] = float(np.nanmax(spectrum[finite_spec]))
                            rec["spectrum_max_abs"] = float(
                                np.nanmax(np.abs(spectrum[finite_spec]))
                            )

                        rec["roughness_logdiff_rms"] = compute_roughness(nu_hz, spectrum)
                        rec["f230ghz_nuLnu"] = compute_f230(nu_hz, spectrum)

                        spectra[str(path)] = {
                            "nu_hz": nu_hz,
                            "spectrum": spectrum,
                        }
                else:
                    append_note(rec, "missing output/lnu or output/nuLnu")
        except Exception as exc:
            append_note(rec, f"open/read exception: {exc}")

        rec["parsed_params"] = json.dumps(parsed_params, sort_keys=True)
        inventory_rows.append(rec.copy())
        sanity_rows.append(rec)

    # Post-pass outlier checks requiring all rows.
    rough_vals = np.array(
        [r["roughness_logdiff_rms"] for r in sanity_rows if r["roughness_logdiff_rms"] is not None],
        dtype=float,
    )
    if rough_vals.size >= 5:
        med = float(np.median(rough_vals))
        mad = float(np.median(np.abs(rough_vals - med)))
        rough_thresh = med + (5.0 * mad if mad > 0 else max(0.15, 2.0 * med))
        for r in sanity_rows:
            rv = r["roughness_logdiff_rms"]
            if rv is not None and rv > rough_thresh:
                append_note(r, f"roughness outlier ({rv:.4g} > {rough_thresh:.4g})")

    max_abs_vals = np.array(
        [r["spectrum_max_abs"] for r in sanity_rows if r["spectrum_max_abs"] is not None and r["spectrum_max_abs"] > 0],
        dtype=float,
    )
    if max_abs_vals.size >= 5:
        med_abs = float(np.median(max_abs_vals))
        huge_thresh = med_abs * 1.0e4
        for r in sanity_rows:
            ma = r["spectrum_max_abs"]
            if ma is not None and med_abs > 0 and ma > huge_thresh:
                append_note(r, f"very large spectrum scale ({ma:.4g} > {huge_thresh:.4g})")

    # Finalize statuses.
    for r in sanity_rows:
        hard = False
        suspicious = False

        if not r["open_ok"]:
            hard = True
            append_note(r, "open failed")
        if r["open_ok"] and not r["has_spectrum"]:
            hard = True
            append_note(r, "missing spectrum datasets")
        if r["nan_count"] > 0:
            hard = True
            append_note(r, "contains NaN values")
        if r["inf_count"] > 0:
            hard = True
            append_note(r, "contains Inf values")
        if r["freq_negative_count"] > 0:
            hard = True
            append_note(r, "non-positive/invalid frequency bins")
        if r["freq_nonmonotonic"]:
            hard = True
            append_note(r, "frequency grid is non-monotonic")
        if r["has_spectrum"] and r["spectrum_finite_count"] == 0:
            hard = True
            append_note(r, "spectrum has no finite bins")
        if r["spectrum_all_zero"]:
            hard = True
            append_note(r, "spectrum is all zeros")
        if r["has_spectrum"] and (r["spectrum_max_abs"] is None or r["spectrum_max_abs"] <= 0):
            hard = True
            append_note(r, "non-positive spectrum scale")

        if r["spectrum_neg_frac"] is not None and r["spectrum_neg_frac"] > 0.2:
            suspicious = True
            append_note(r, f"high negative-bin fraction ({r['spectrum_neg_frac']:.2%})")
        if r["pol_products_present"] and r["pol_max_abs_frac"] is not None and r["pol_max_abs_frac"] > 1.5:
            suspicious = True
            append_note(r, f"polarization fraction exceeds 1.5 ({r['pol_max_abs_frac']:.3g})")
        if r["pol_products_present"] and r["pol_nonfinite_count"] > 0:
            suspicious = True
            append_note(r, "non-finite polarization fraction bins")
        if any("outlier" in n for n in r.get("_notes", [])):
            suspicious = True

        r["hard_fail"] = hard
        r["suspicious"] = suspicious
        r["pass_sanity"] = bool(r["open_ok"] and r["has_spectrum"] and not hard and not suspicious)
        r["notes"] = "; ".join(r.get("_notes", []))

    # Keep inventory in sync with finalized columns.
    inv_by_path = {r["filepath"]: r for r in sanity_rows}
    for idx, row in enumerate(inventory_rows):
        latest = inv_by_path[row["filepath"]]
        inventory_rows[idx] = latest.copy()

    # Quicklook plots (all readable spectra).
    for r in sanity_rows:
        spec = spectra.get(r["filepath"])
        if not spec:
            continue
        out_png = PLOTS_DIR / f"quicklook_{safe_slug(Path(r['filename']).stem)}.png"
        plot_quicklook(r, spec["nu_hz"], spec["spectrum"], out_png)

    # Pairwise checks.
    pair_rows = []

    # Heating pairs: same config except heating model.
    heating_groups = defaultdict(lambda: defaultdict(list))
    for r in sanity_rows:
        if r["heating_model"] not in {"RBETA", "CRITBETA"}:
            continue
        if r["filepath"] not in spectra:
            continue
        key = (
            r["state"],
            fmt_spin_for_key(r["spin"]),
            r["dump_time"],
            r["wjet_flag"],
            r["pos"],
        )
        heating_groups[key][r["heating_model"]].append(r)

    for key, hm in sorted(
        heating_groups.items(),
        key=lambda kv: tuple("" if v is None else str(v) for v in kv[0]),
    ):
        if "RBETA" not in hm or "CRITBETA" not in hm:
            continue
        rbeta = choose_representative(hm["RBETA"])
        crit = choose_representative(hm["CRITBETA"])
        a = spectra[rbeta["filepath"]]
        b = spectra[crit["filepath"]]
        cmp_stats = compare_spectra(a["nu_hz"], a["spectrum"], b["nu_hz"], b["spectrum"])

        state, spin_key, dump_time, wjet_flag, pos = key
        title = (
            f"Heating Pair: {state} a={spin_key} t={dump_time} "
            f"WJET={wjet_flag} pos={pos}"
        )
        pair_png = PAIR_PLOTS_DIR / (
            f"pair_heating_state-{state}_a-{spin_key}_t-{dump_time}_"
            f"wjet-{wjet_flag}_pos-{pos}.png"
        )
        plot_pair(
            title,
            f"RBETA: {Path(rbeta['filename']).stem}",
            a["nu_hz"],
            a["spectrum"],
            f"CRITBETA: {Path(crit['filename']).stem}",
            b["nu_hz"],
            b["spectrum"],
            pair_png,
        )

        pair_rows.append(
            {
                "pair_type": "heating",
                "state": state,
                "spin": spin_key,
                "dump_time": dump_time,
                "heating_model": "RBETA_vs_CRITBETA",
                "wjet_flag": wjet_flag,
                "pos": pos,
                "file_ref": rbeta["filepath"],
                "file_cmp": crit["filepath"],
                "median_ratio_cmp_over_ref": cmp_stats["median_ratio"],
                "p05_ratio": cmp_stats["p05_ratio"],
                "p95_ratio": cmp_stats["p95_ratio"],
                "max_abs_log10_ratio": cmp_stats["max_abs_log10_ratio"],
                "wide_jump_frac": cmp_stats["wide_jump_frac"],
                "valid_bins": cmp_stats["valid_bins"],
                "suspicious_pair": cmp_stats["suspicious_pair"],
                "notes": cmp_stats["notes"],
                "plot_path": str(pair_png),
            }
        )

    # WJET pairs: same config except WJET off/on.
    wjet_groups = defaultdict(lambda: defaultdict(list))
    for r in sanity_rows:
        if r["heating_model"] == "UNKNOWN":
            continue
        if r["filepath"] not in spectra:
            continue
        key = (
            r["state"],
            fmt_spin_for_key(r["spin"]),
            r["dump_time"],
            r["heating_model"],
            r["pos"],
        )
        wjet_groups[key][r["wjet_flag"]].append(r)

    for key, flags in sorted(
        wjet_groups.items(),
        key=lambda kv: tuple("" if v is None else str(v) for v in kv[0]),
    ):
        if "OFF" not in flags or "ON" not in flags:
            continue
        off = choose_representative(flags["OFF"])
        on = choose_representative(flags["ON"])
        a = spectra[off["filepath"]]
        b = spectra[on["filepath"]]
        cmp_stats = compare_spectra(a["nu_hz"], a["spectrum"], b["nu_hz"], b["spectrum"])

        state, spin_key, dump_time, heating_model, pos = key
        title = (
            f"WJET Pair: {state} a={spin_key} t={dump_time} "
            f"{heating_model} pos={pos}"
        )
        pair_png = PAIR_PLOTS_DIR / (
            f"pair_wjet_state-{state}_a-{spin_key}_t-{dump_time}_"
            f"heat-{heating_model}_pos-{pos}.png"
        )
        plot_pair(
            title,
            f"WJET OFF: {Path(off['filename']).stem}",
            a["nu_hz"],
            a["spectrum"],
            f"WJET ON: {Path(on['filename']).stem}",
            b["nu_hz"],
            b["spectrum"],
            pair_png,
        )

        pair_rows.append(
            {
                "pair_type": "wjet",
                "state": state,
                "spin": spin_key,
                "dump_time": dump_time,
                "heating_model": heating_model,
                "wjet_flag": "OFF_vs_ON",
                "pos": pos,
                "file_ref": off["filepath"],
                "file_cmp": on["filepath"],
                "median_ratio_cmp_over_ref": cmp_stats["median_ratio"],
                "p05_ratio": cmp_stats["p05_ratio"],
                "p95_ratio": cmp_stats["p95_ratio"],
                "max_abs_log10_ratio": cmp_stats["max_abs_log10_ratio"],
                "wide_jump_frac": cmp_stats["wide_jump_frac"],
                "valid_bins": cmp_stats["valid_bins"],
                "suspicious_pair": cmp_stats["suspicious_pair"],
                "notes": cmp_stats["notes"],
                "plot_path": str(pair_png),
            }
        )

    # Write CSV outputs.
    inventory_fields = [
        "filepath",
        "filename",
        "run_id",
        "state_letter",
        "state",
        "spin",
        "spin_source",
        "dump_time",
        "pos",
        "heating_model",
        "wjet_flag",
        "wjet_label",
        "trial_tag",
        "m_unit",
        "mbh",
        "tp_over_te",
        "bias",
        "rmax",
        "rhigh",
        "beta_crit",
        "param_keys_count",
        "header_keys_count",
        "top_level_keys",
        "parsed_params",
        "pass_sanity",
    ]
    write_csv(INVENTORY_CSV, inventory_rows, inventory_fields)

    sanity_fields = [
        "filepath",
        "filename",
        "parsed_params",
        "open_ok",
        "has_spectrum",
        "nan_count",
        "inf_count",
        "lnu_min",
        "lnu_max",
        "nu_min_hz",
        "nu_max_hz",
        "freq_negative_count",
        "freq_nonmonotonic",
        "spectrum_bins",
        "spectrum_finite_count",
        "spectrum_positive_count",
        "spectrum_negative_count",
        "spectrum_neg_frac",
        "spectrum_all_zero",
        "spectrum_min",
        "spectrum_max",
        "spectrum_max_abs",
        "roughness_logdiff_rms",
        "f230ghz_nuLnu",
        "pol_products_present",
        "pol_max_abs_frac",
        "pol_nonfinite_count",
        "hard_fail",
        "suspicious",
        "pass_sanity",
        "notes",
    ]
    write_csv(SANITY_CSV, sanity_rows, sanity_fields)

    pair_fields = [
        "pair_type",
        "state",
        "spin",
        "dump_time",
        "heating_model",
        "wjet_flag",
        "pos",
        "file_ref",
        "file_cmp",
        "median_ratio_cmp_over_ref",
        "p05_ratio",
        "p95_ratio",
        "max_abs_log10_ratio",
        "wide_jump_frac",
        "valid_bins",
        "suspicious_pair",
        "notes",
        "plot_path",
    ]
    write_csv(PAIRWISE_CSV, pair_rows, pair_fields)

    # Summaries for report.
    total = len(sanity_rows)
    open_fail = sum(not r["open_ok"] for r in sanity_rows)
    hard_fail = sum(r["hard_fail"] for r in sanity_rows)
    suspicious = sum(r["suspicious"] for r in sanity_rows)
    passed = sum(r["pass_sanity"] for r in sanity_rows)
    failed = total - passed

    def summarize_by(field):
        c_total: Counter[str] = Counter()
        c_pass: Counter[str] = Counter()
        for row in sanity_rows:
            key = str(row.get(field, "UNKNOWN"))
            c_total[key] += 1
            if row["pass_sanity"]:
                c_pass[key] += 1
        keys = sorted(c_total.keys())
        return [(k, c_total[k], c_pass[k], c_total[k] - c_pass[k]) for k in keys]

    by_heating = summarize_by("heating_model")
    by_wjet = summarize_by("wjet_flag")
    by_state = summarize_by("state")
    by_spin = summarize_by("spin")

    def offender_score(row):
        score = 0
        if not row["open_ok"]:
            score += 100
        if row["open_ok"] and not row["has_spectrum"]:
            score += 80
        if row["nan_count"] > 0 or row["inf_count"] > 0:
            score += 70
        if row["spectrum_all_zero"]:
            score += 60
        if row["freq_negative_count"] > 0 or row["freq_nonmonotonic"]:
            score += 50
        if row["suspicious"]:
            score += 20
        score += len(row.get("_notes", []))
        return score

    offenders = [r for r in sanity_rows if not r["pass_sanity"]]
    offenders = sorted(offenders, key=offender_score, reverse=True)[:10]

    heating_pairs = [p for p in pair_rows if p["pair_type"] == "heating"]
    wjet_pairs = [p for p in pair_rows if p["pair_type"] == "wjet"]
    bad_heating_pairs = [p for p in heating_pairs if p["suspicious_pair"]]
    bad_wjet_pairs = [p for p in wjet_pairs if p["suspicious_pair"]]

    wjet_ratio_vals = [
        p["median_ratio_cmp_over_ref"]
        for p in wjet_pairs
        if p["median_ratio_cmp_over_ref"] is not None
    ]
    wjet_boost_count = sum(v > 1.0 for v in wjet_ratio_vals)
    wjet_dim_count = sum(v < 1.0 for v in wjet_ratio_vals)

    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = []
    lines.append("# GRMONTY QA Report")
    lines.append("")
    lines.append(f"Generated: {ts}")
    lines.append(f"Base directory: `{BASE_DIR}`")
    if crossref_enabled:
        lines.append(
            "Selection mode: `crossref-logs`"
            + (" + `canonical-only`" if args.canonical_only else "")
        )
        lines.append(f"Logs directory: `{args.logs_dir}`")
        lines.append(f"Log configs discovered: **{len(log_cfg_set)}**")
        lines.append(
            f"Filtered out non-canonical filenames: **{len(filtered_out_noncanonical)}**"
        )
        lines.append(
            f"Filtered out not matched to log configs: **{len(filtered_out_not_in_logs)}**"
        )
        lines.append(
            f"Log configs with no canonical `.h5` found: **{len(missing_log_cfgs)}**"
        )
    else:
        lines.append("Selection mode: `all-h5` (no log cross-reference)")
    lines.append("")

    if crossref_enabled and missing_log_cfgs:
        lines.append("Missing canonical outputs for log configs:")
        for cfg in missing_log_cfgs:
            lines.append(f"- `{cfg}`")
        lines.append("")
    lines.append("## Overall Counts")
    lines.append("")
    lines.append(f"- Total `.h5` outputs found: **{total}**")
    lines.append(f"- Open failures: **{open_fail}**")
    lines.append(f"- Hard sanity failures: **{hard_fail}**")
    lines.append(f"- Suspicious (soft) failures: **{suspicious}**")
    lines.append(f"- Pass sanity checks (`looks good`): **{passed}**")
    lines.append(f"- Total non-pass files: **{failed}**")
    lines.append("")
    lines.append("Artifacts:")
    lines.append(f"- Inventory CSV: `{INVENTORY_CSV}`")
    lines.append(f"- Sanity CSV: `{SANITY_CSV}`")
    lines.append(f"- Pairwise CSV: `{PAIRWISE_CSV}`")
    lines.append(f"- Plot directory: `{PLOTS_DIR}`")
    lines.append("")

    lines.append("## Top 10 Worst Offenders")
    lines.append("")
    if offenders:
        for i, row in enumerate(offenders, start=1):
            why = row["notes"] if row["notes"] else "unspecified"
            lines.append(f"{i}. `{row['filepath']}`")
            lines.append(f"   - reason: {why}")
    else:
        lines.append("- No offenders. All files passed current checks.")
    lines.append("")

    def append_summary_table(title, rows):
        lines.append(f"## {title}")
        lines.append("")
        lines.append("| Group | Total | Pass | Fail |")
        lines.append("|---|---:|---:|---:|")
        for key, t, p, f in rows:
            lines.append(f"| {key} | {t} | {p} | {f} |")
        lines.append("")

    append_summary_table("Summary by Heating Model", by_heating)
    append_summary_table("Summary by WJET Flag", by_wjet)
    append_summary_table("Summary by State (MAD/SANE)", by_state)
    append_summary_table("Summary by Spin", by_spin)

    lines.append("## Pairwise Comparative Checks")
    lines.append("")
    lines.append(
        f"- Heating pairs (`RBETA` vs `CRITBETA`) evaluated: **{len(heating_pairs)}**"
    )
    lines.append(f"- Suspicious heating pairs: **{len(bad_heating_pairs)}**")
    lines.append(f"- WJET OFF/ON pairs evaluated: **{len(wjet_pairs)}**")
    lines.append(f"- Suspicious WJET pairs: **{len(bad_wjet_pairs)}**")
    if wjet_pairs:
        lines.append(
            "- WJET behavior summary (based on median ON/OFF ratio): "
            f"ON>OFF in **{wjet_boost_count}** pairs, ON<OFF in **{wjet_dim_count}** pairs."
        )
    lines.append("")

    if bad_heating_pairs or bad_wjet_pairs:
        lines.append("Flagged pair details:")
        for p in bad_heating_pairs[:10] + bad_wjet_pairs[:10]:
            lines.append(
                f"- `{p['pair_type']}` `{p['file_ref']}` vs `{p['file_cmp']}`: "
                f"{p['notes'] or 'ratio outlier'}"
            )
        lines.append("")

    lines.append("## Common Problem Patterns")
    lines.append("")
    if failed == 0:
        lines.append("- No failure patterns detected under current sanity thresholds.")
    else:
        fail_rows = [r for r in sanity_rows if not r["pass_sanity"]]
        token_counts: Counter[str] = Counter()
        for row in fail_rows:
            for token in row.get("_notes", []):
                token_counts[token] += 1
        for reason, cnt in token_counts.most_common(12):
            lines.append(f"- {reason}: {cnt} files")
    lines.append("")

    REPORT_MD.write_text("\n".join(lines) + "\n")

    print("[qa] complete")
    print(f"[qa] inventory: {INVENTORY_CSV}")
    print(f"[qa] sanity:    {SANITY_CSV}")
    print(f"[qa] pairwise:  {PAIRWISE_CSV}")
    print(f"[qa] report:    {REPORT_MD}")
    print(f"[qa] plots:     {PLOTS_DIR}")


if __name__ == "__main__":
    main()
