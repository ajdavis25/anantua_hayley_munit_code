#!/usr/bin/env python3
"""MWL SED check: existing production GRMONTY spectra vs observed M87 core anchors.

Motivation (2026-08-06): committee member noted Doppler boosting of the M87 jet
is verified in radio/optical/X-ray. The observational anchors from the same
literature constrain the wJET supplement directly:

  - 230 GHz: tuning target 0.5 Jy (auto_munit_bracket F_TARGET_JY, EHT 2017
    compact-flux scale).
  - 2-10 keV core: L_X = 4.4e40 erg/s (Chandra+NuSTAR quasi-simultaneous,
    EHT MWL WG 2021, ApJL 911 L11; core dominates HST-1 in the 2017 low state).

For each top-level spectrum in /work/vmo703/igrmonty_outputs/m87 this script
computes, using the exact conventions of tools/viewing_cone_postprocess.py
(theta folded to [0,90] deg, dOmega-weighted, L_sun -> cgs):

  - baseline_4pi        : tuner-style all-sky-average spectrum
  - cone_renorm_4pi_equiv: observer-at-17-deg (+-10 deg cone) isotropic equiv
  - F_nu at 86/230/345 GHz per mode, L_X(2-10 keV), L_bol, Compton split
  - M_unit provenance cross-check vs /work/vmo703/data/final_grmonty_paper.csv

Outputs:
  _qa/mwl_sed_check.csv        one row per file x mode
  _qa/mwl_sed_verdicts.csv     one row per file (gates + provenance)

Run:  /work/vmo703/ipole_venv/bin/python mwl_sed_check.py
"""

import argparse
import csv
import importlib.util
import math
import os
import re
import sys

import h5py
import numpy as np

TOOL_PATH = "/work/vmo703/igrmonty/tools/viewing_cone_postprocess.py"
BASE_DIR = "/work/vmo703/igrmonty_outputs/m87"
QA_DIR = os.path.join(BASE_DIR, "_qa")
FROZEN_CSV = "/work/vmo703/data/final_grmonty_paper.csv"

OUT_MODES_CSV = os.path.join(QA_DIR, "mwl_sed_check.csv")
OUT_VERDICT_CSV = os.path.join(QA_DIR, "mwl_sed_verdicts.csv")

THETACAM_DEG = 17.0
CONE_HALF_ANGLE_DEG = 10.0
DISTANCE_MPC = 16.8

F230_TARGET_JY = 0.5           # auto_munit_bracket F_TARGET_JY at 230 GHz
LX_OBS_ERG_S = 4.4e40          # 2-10 keV core, EHT MWL WG 2021
LX_OVERSHOOT_FACTOR = 1.5      # MC-noise-tolerant overshoot threshold

KEV_HZ = 2.417989242e17        # 1 keV in Hz
XBAND_LO_HZ = 2.0 * KEV_HZ
XBAND_HI_HZ = 10.0 * KEV_HZ
NIR_HZ = 1.36e14               # 2.2 micron (K band), secondary reference only


def load_tool():
    spec = importlib.util.spec_from_file_location("vcp", TOOL_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def parse_name(name):
    out = {
        "state": None, "spin": None, "dump": None,
        "heating": None, "wjet": None, "pos": None, "trial": None,
    }
    m = re.search(r"_([SM])a([+-]?\d+(?:\.\d+)?)_(\d+)_", name)
    if m:
        out["state"] = "SANE" if m.group(1) == "S" else "MAD"
        out["spin"] = float(m.group(2))
        out["dump"] = int(m.group(3))
    if "CRITBETA" in name:
        out["heating"] = "CRITBETA"
    elif "RBETA" in name:
        out["heating"] = "RBETA"
    out["wjet"] = "ON" if "wJET" in name else "OFF"
    pm = re.search(r"_pos(\d+)", name)
    if pm:
        out["pos"] = int(pm.group(1))
    tm = re.search(r"(trial\d+|TEST)", name)
    if tm:
        out["trial"] = tm.group(1)
    return out


def model_key(meta):
    if None in (meta["state"], meta["heating"], meta["wjet"]):
        return None
    return meta["heating"] + ("wJET" if meta["wjet"] == "ON" else "")


def load_frozen_munits(path):
    rows = {}
    if not os.path.exists(path):
        return rows
    with open(path) as fh:
        for row in csv.DictReader(fh):
            try:
                key = (
                    row["state"].strip(),
                    row["model"].strip(),
                    float(row["spin"]),
                    int(row["dump_index"]),
                )
            except Exception:
                continue
            rows[key] = row
    return rows


def band_luminosity(nu_hz, nuLnu_cgs, lo_hz, hi_hz):
    """L(band) = int Lnu dnu = int nuLnu dln(nu) over [lo, hi]."""
    lnnu = np.log(nu_hz)
    mask = (nu_hz >= lo_hz) & (nu_hz <= hi_hz)
    n_bins = int(np.count_nonzero(mask))
    n_pos = int(np.count_nonzero(mask & (nuLnu_cgs > 0.0)))
    if n_bins < 2:
        return float("nan"), n_bins, n_pos
    val = float(np.trapz(nuLnu_cgs[mask], lnnu[mask]))
    return val, n_bins, n_pos


def read_file_extras(path):
    out = {}
    with h5py.File(path, "r") as f:
        p = f["params"]
        for key in ("M_unit", "Ns", "THETAE_MAX", "TP_OVER_TE", "MBH", "bias"):
            try:
                out[key] = float(p[key][()])
            except Exception:
                out[key] = None
        e = p.get("electrons")
        if e is not None:
            for key in (
                "type", "sigma_transition", "constant_beta_e0",
                "constant_beta_e0_exponent", "jet_sigma_cut", "jet_beta_cut",
                "jet_thetae", "jet_ne_mult", "positron_ratio",
            ):
                try:
                    out["e_" + key] = float(e[key][()])
                except Exception:
                    out["e_" + key] = None
        o = f["output"]
        for key in ("L", "Mdot", "MdotEdd", "efficiency", "Nmade", "Nrecorded", "Nscattered"):
            try:
                out["out_" + key] = float(o[key][()])
            except Exception:
                out["out_" + key] = None
        try:
            out["Lcomponent"] = np.array(o["Lcomponent"], dtype=float)
        except Exception:
            out["Lcomponent"] = None
        try:
            rs = o["run_status"][()]
            out["run_status"] = rs.decode() if isinstance(rs, bytes) else str(rs)
        except Exception:
            out["run_status"] = None
    return out


def main():
    ap = argparse.ArgumentParser(description="MWL SED check on a GRMONTY spectrum corpus")
    ap.add_argument(
        "--dir",
        default=os.path.join(BASE_DIR, "pre_bugfix_2026-07"),
        help="Directory of spectrum_*.h5 files (default: archived pre-fix corpus; "
             "point at an era subdir, e.g. .../m87/postfix_2026-08, for new runs). "
             "NOTE: output CSV names are fixed, so a rerun overwrites the previous era's CSVs.",
    )
    args = ap.parse_args()
    corpus_dir = args.dir

    vcp = load_tool()

    files = []
    for name in sorted(os.listdir(corpus_dir)):
        full = os.path.join(corpus_dir, name)
        if not os.path.isfile(full):
            continue
        if not name.endswith(".h5") or "_cone_" in name:
            continue
        files.append(full)

    frozen = load_frozen_munits(FROZEN_CSV)
    print(f"[info] {len(files)} spectra, {len(frozen)} frozen-M_unit rows")

    mode_rows = []
    verdict_rows = []

    for path in files:
        name = os.path.basename(path)
        meta = parse_name(name)
        try:
            data = vcp.load_spectrum(path)
        except Exception as exc:
            print(f"[fail] {name} :: {exc}", file=sys.stderr)
            verdict_rows.append({"filename": name, "error": str(exc)})
            continue

        nu = data["nu_hz"]
        n_thbins = int(data["dOmega_sr"].shape[0])
        _, centers_deg, _ = vcp.build_theta_bins(n_thbins)
        theta_meta = vcp.select_theta_bins(centers_deg, THETACAM_DEG, CONE_HALF_ANGLE_DEG)
        modes = vcp.compute_spectra_modes(
            nu, data["nuLnu_theta_freq_cgs"], data["dOmega_sr"],
            theta_meta["theta_mask"], DISTANCE_MPC,
        )
        extras = read_file_extras(path)

        mode_map = [
            ("baseline_4pi", modes["nuLnu_4pi_cgs"], modes["fnu_4pi_jy"]),
            ("cone_renorm_4pi_equiv", modes["nuLnu_cone_renorm_cgs"], modes["fnu_cone_renorm_jy"]),
            ("cone_physical", modes["nuLnu_cone_phys_cgs"], modes["fnu_cone_phys_jy"]),
        ]

        per_mode = {}
        for mode_name, nuLnu_mode, fnu_mode in mode_map:
            lx, nb, npos = band_luminosity(nu, nuLnu_mode, XBAND_LO_HZ, XBAND_HI_HZ)
            lbol, _, _ = band_luminosity(nu, nuLnu_mode, nu[0], nu[-1])
            row = {
                "filename": name,
                "state": meta["state"], "spin": meta["spin"], "dump": meta["dump"],
                "heating": meta["heating"], "wjet": meta["wjet"], "pos": meta["pos"],
                "trial": meta["trial"], "mode": mode_name,
                "fnu_86_jy": vcp.interpolate_fnu(nu, fnu_mode, 86.0e9, "linlog"),
                "fnu_230_jy": vcp.interpolate_fnu(nu, fnu_mode, 230.0e9, "linlog"),
                "fnu_345_jy": vcp.interpolate_fnu(nu, fnu_mode, 345.0e9, "linlog"),
                "fnu_nir_2p2um_jy": vcp.interpolate_fnu(nu, fnu_mode, NIR_HZ, "linlog"),
                "L_x_2_10_keV": lx,
                "L_x_nbins": nb,
                "L_x_npos_bins": npos,
                "L_bol_mode": lbol,
            }
            per_mode[mode_name] = row
            mode_rows.append(row)

        # provenance: M_unit vs frozen table
        munit_file = extras.get("M_unit")
        munit_frozen = None
        frozen_conv = None
        mk = model_key(meta)
        if mk is not None and meta["pos"] is not None:
            row = frozen.get((meta["state"], mk, meta["spin"], meta["dump"]))
            if row is not None:
                col = f"MunitUsed_pos{meta['pos']}"
                try:
                    munit_frozen = float(row.get(col, "nan"))
                except Exception:
                    munit_frozen = None
                frozen_conv = row.get("converged")
        munit_match = None
        if munit_file and munit_frozen and munit_frozen > 0:
            munit_match = munit_file / munit_frozen

        # Compton split under layout [synch x (0,1,2,>2), brems x (0,1,2,>2)]
        lcomp = extras.get("Lcomponent")
        frac_scattered = None
        lcomp_str = ""
        if lcomp is not None and lcomp.size == 8 and np.sum(lcomp) > 0:
            tot = float(np.sum(lcomp))
            frac_scattered = 1.0 - (float(lcomp[0]) + float(lcomp[4])) / tot
            lcomp_str = ";".join(f"{v:.4e}" for v in lcomp)

        base = per_mode["baseline_4pi"]
        cone = per_mode["cone_renorm_4pi_equiv"]
        lx_ratio = cone["L_x_2_10_keV"] / LX_OBS_ERG_S if np.isfinite(cone["L_x_2_10_keV"]) else float("nan")
        gate_x = "OVERSHOOT" if (np.isfinite(lx_ratio) and lx_ratio > LX_OVERSHOOT_FACTOR) else "OK"

        verdict_rows.append({
            "filename": name,
            "state": meta["state"], "spin": meta["spin"], "dump": meta["dump"],
            "heating": meta["heating"], "wjet": meta["wjet"], "pos": meta["pos"],
            "trial": meta["trial"],
            "Ns": extras.get("Ns"),
            "Nrecorded": extras.get("out_Nrecorded"),
            "run_status": extras.get("run_status"),
            "thetae_max": extras.get("THETAE_MAX"),
            "e_type": extras.get("e_type"),
            "e_sigma_transition": extras.get("e_sigma_transition"),
            "e_jet_sigma_cut": extras.get("e_jet_sigma_cut"),
            "e_jet_thetae": extras.get("e_jet_thetae"),
            "e_positron_ratio": extras.get("e_positron_ratio"),
            "M_unit_file": munit_file,
            "M_unit_frozen": munit_frozen,
            "M_unit_file_over_frozen": munit_match,
            "frozen_converged": frozen_conv,
            "f230_base_jy": base["fnu_230_jy"],
            "f230_base_over_target": base["fnu_230_jy"] / F230_TARGET_JY if np.isfinite(base["fnu_230_jy"]) else float("nan"),
            "f230_cone_jy": cone["fnu_230_jy"],
            "f230_cone_over_target": cone["fnu_230_jy"] / F230_TARGET_JY if np.isfinite(cone["fnu_230_jy"]) else float("nan"),
            "f86_cone_jy": cone["fnu_86_jy"],
            "f345_cone_jy": cone["fnu_345_jy"],
            "fnu_nir_cone_jy": cone["fnu_nir_2p2um_jy"],
            "Lx_cone_erg_s": cone["L_x_2_10_keV"],
            "Lx_cone_over_obs": lx_ratio,
            "Lx_cone_npos_bins": cone["L_x_npos_bins"],
            "Lx_base_erg_s": base["L_x_2_10_keV"],
            "L_bol_cone": cone["L_bol_mode"],
            "L_bol_base": base["L_bol_mode"],
            "L_total_output": extras.get("out_L"),
            "frac_scattered_Lcomp": frac_scattered,
            "Lcomponent": lcomp_str,
            "gate_Lx": gate_x,
            "error": "",
        })
        print(
            f"[ok] {name}: F230(base)={base['fnu_230_jy']:.3g} Jy "
            f"F230(cone17)={cone['fnu_230_jy']:.3g} Jy "
            f"Lx(cone17)={cone['L_x_2_10_keV']:.3g} erg/s "
            f"(x{lx_ratio:.3g} obs) [{gate_x}]"
        )

    def write_rows(path, rows):
        fields = []
        for r in rows:
            for k in r:
                if k not in fields:
                    fields.append(k)
        with open(path, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow(r)

    write_rows(OUT_MODES_CSV, mode_rows)
    write_rows(OUT_VERDICT_CSV, verdict_rows)
    print(f"[done] {OUT_MODES_CSV}")
    print(f"[done] {OUT_VERDICT_CSV}")


if __name__ == "__main__":
    main()
