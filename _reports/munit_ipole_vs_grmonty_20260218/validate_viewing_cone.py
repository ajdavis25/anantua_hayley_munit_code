#!/usr/bin/env python3
"""Validate GRMONTY cone postprocessing outputs against baseline 4pi spectra.

This script:
- selects 3-5 GRMONTY runs from report `tables/run_manifest.csv`,
- runs `igrmonty/tools/viewing_cone_postprocess.py` on each,
- writes:
  * cone_validation_manifest.csv
  * cone_validation_plots/*.png

Usage:
  python3 /work/vmo703/_reports/munit_ipole_vs_grmonty_20260218/validate_viewing_cone.py
"""

from __future__ import print_function

import csv
import os
import subprocess
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


THETACAM_DEG = 17.0
CONE_HALF_ANGLE_DEG = 10.0
MAX_RUNS = 5
ANCHOR_HZ = 228.0e9


def _float_tag(val):
    s = "{:.1f}".format(float(val))
    s = s.replace("-", "m").replace(".", "p")
    return s


def root_dir():
    here = os.path.abspath(os.path.dirname(__file__))
    return os.path.abspath(os.path.join(here, "..", ".."))


def report_dir():
    return os.path.abspath(os.path.dirname(__file__))


def tool_path():
    return os.path.join(root_dir(), "igrmonty", "tools", "viewing_cone_postprocess.py")


def run_manifest_path():
    return os.path.join(report_dir(), "tables", "run_manifest.csv")


def pick_runs(csv_path, max_runs):
    rows = []
    with open(csv_path, "r", newline="") as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            if row.get("code", "") != "GRMONTY":
                continue
            p = row.get("path", "")
            if not p or not os.path.exists(p):
                continue
            rows.append(row)

    rows = sorted(
        rows,
        key=lambda r: (
            r.get("state", ""),
            r.get("model", ""),
            float(r.get("spin", "0") or "0"),
            int(r.get("dump_index", "0") or "0"),
            int(r.get("pos", "0") or "0"),
        ),
    )

    picked = []
    seen = set()
    # Prefer diversity by (state, model, spin)
    for row in rows:
        key = (row.get("state", ""), row.get("model", ""), row.get("spin", ""))
        if key in seen:
            continue
        picked.append(row)
        seen.add(key)
        if len(picked) >= max_runs:
            return picked

    # Fill remaining slots if needed.
    for row in rows:
        if row in picked:
            continue
        picked.append(row)
        if len(picked) >= max_runs:
            break

    return picked


def sidecar_paths(spec_path, thetacam_deg, cone_half_angle_deg):
    root, _ext = os.path.splitext(spec_path)
    tc_folded = float(thetacam_deg)
    if tc_folded > 180.0:
        tc_folded = 360.0 - tc_folded
    if tc_folded > 90.0:
        tc_folded = 180.0 - tc_folded
    tag = "tc{}_a{}".format(_float_tag(tc_folded), _float_tag(cone_half_angle_deg))
    npz_path = "{}_cone_{}.npz".format(root, tag)
    csv_path = "{}_cone_summary_{}.csv".format(root, tag)
    return npz_path, csv_path


def ensure_postprocess(spec_path):
    tool = tool_path()
    cmd = [
        "python3",
        tool,
        spec_path,
        "--thetacam-deg",
        str(THETACAM_DEG),
        "--cone-half-angle-deg",
        str(CONE_HALF_ANGLE_DEG),
        "--anchor-ghz",
        "228",
        "--overwrite",
    ]
    subprocess.check_call(cmd)
    return sidecar_paths(spec_path, THETACAM_DEG, CONE_HALF_ANGLE_DEG)


def read_summary_rows(csv_path):
    with open(csv_path, "r", newline="") as fh:
        rows = list(csv.DictReader(fh))
    out = {}
    for row in rows:
        out[row.get("mode", "")] = row
    return out


def safe_float(x):
    try:
        return float(x)
    except Exception:
        return float("nan")


def make_plots(npz_path, plot_dir):
    data = np.load(npz_path)
    nu = np.array(data["nu_hz"], dtype=np.float64)
    f0 = np.array(data["fnu_4pi_jy"], dtype=np.float64)
    f1 = np.array(data["fnu_cone_phys_jy"], dtype=np.float64)
    f2 = np.array(data["fnu_cone_renorm_jy"], dtype=np.float64)

    run_base = os.path.basename(npz_path).replace("_cone_", "_").replace(".npz", "")

    overlay_path = os.path.join(plot_dir, "{}_sed_overlay.png".format(run_base))
    ratio_path = os.path.join(plot_dir, "{}_ratio.png".format(run_base))

    plt.figure(figsize=(8.0, 5.5))
    plt.loglog(nu, f0, lw=2.0, label="baseline 4pi")
    plt.loglog(nu, f1, lw=1.8, label="cone physical (no renorm)")
    plt.loglog(nu, f2, lw=1.8, label="cone renorm 4pi-equiv")
    plt.axvline(ANCHOR_HZ, color="k", ls="--", lw=1.0, alpha=0.7, label="228 GHz anchor")
    plt.xlabel("frequency [Hz]")
    plt.ylabel("F_nu [Jy]")
    plt.title(run_base)
    plt.grid(True, which="both", ls=":", alpha=0.4)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(overlay_path, dpi=140)
    plt.close()

    with np.errstate(divide="ignore", invalid="ignore"):
        r1 = np.where(f0 > 0.0, f1 / f0, np.nan)
        r2 = np.where(f0 > 0.0, f2 / f0, np.nan)

    plt.figure(figsize=(8.0, 5.0))
    plt.semilogx(nu, r1, lw=1.8, label="cone physical / baseline")
    plt.semilogx(nu, r2, lw=1.8, label="cone renorm / baseline")
    plt.axhline(1.0, color="k", ls="--", lw=1.0, alpha=0.6)
    plt.axvline(ANCHOR_HZ, color="k", ls=":", lw=1.0, alpha=0.6)
    plt.xlabel("frequency [Hz]")
    plt.ylabel("ratio")
    plt.title(run_base)
    plt.grid(True, which="both", ls=":", alpha=0.4)
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(ratio_path, dpi=140)
    plt.close()

    return overlay_path, ratio_path


def main():
    rdir = report_dir()
    manifest_csv = run_manifest_path()
    out_manifest = os.path.join(rdir, "cone_validation_manifest.csv")
    plot_dir = os.path.join(rdir, "cone_validation_plots")

    if not os.path.exists(plot_dir):
        os.makedirs(plot_dir)

    runs = pick_runs(manifest_csv, MAX_RUNS)
    if len(runs) == 0:
        raise RuntimeError("No existing GRMONTY runs found in {}".format(manifest_csv))

    rows_out = []
    for row in runs:
        spec_path = row.get("path", "")
        npz_path, csv_path = ensure_postprocess(spec_path)
        by_mode = read_summary_rows(csv_path)
        overlay_path, ratio_path = make_plots(npz_path, plot_dir)

        base = by_mode.get("baseline_4pi", {})
        phys = by_mode.get("cone_physical", {})
        renm = by_mode.get("cone_renorm_4pi_equiv", {})

        b230 = safe_float(base.get("fnu_230_nearest_jy", "nan"))
        p230 = safe_float(phys.get("fnu_230_nearest_jy", "nan"))
        r230 = safe_float(renm.get("fnu_230_nearest_jy", "nan"))

        out = {
            "state": row.get("state", ""),
            "model": row.get("model", ""),
            "spin": row.get("spin", ""),
            "dump_index": row.get("dump_index", ""),
            "pos": row.get("pos", ""),
            "Munit": row.get("Munit", ""),
            "input_h5": spec_path,
            "cone_npz": npz_path,
            "cone_summary_csv": csv_path,
            "plot_sed_overlay": overlay_path,
            "plot_ratio": ratio_path,
            "thetacam_deg": THETACAM_DEG,
            "cone_half_angle_deg": CONE_HALF_ANGLE_DEG,
            "anchor_nearest_hz": base.get("anchor_nearest_hz", ""),
            "anchor_offset_pct": base.get("anchor_offset_pct", ""),
            "selected_theta_indices": base.get("selected_theta_indices", ""),
            "selected_omega_sr": base.get("selected_omega_sr", ""),
            "fnu_230_baseline_nearest_jy": b230,
            "fnu_230_cone_phys_nearest_jy": p230,
            "fnu_230_cone_renorm_nearest_jy": r230,
            "ratio_phys_to_baseline_230": (p230 / b230) if (b230 and np.isfinite(b230)) else float("nan"),
            "ratio_renorm_to_baseline_230": (r230 / b230) if (b230 and np.isfinite(b230)) else float("nan"),
            "fnu_228_baseline_interp_jy": base.get("fnu_anchor_interp_jy", ""),
            "fnu_228_cone_phys_interp_jy": phys.get("fnu_anchor_interp_jy", ""),
            "fnu_228_cone_renorm_interp_jy": renm.get("fnu_anchor_interp_jy", ""),
        }
        rows_out.append(out)

    # Write manifest
    fieldnames = []
    for r in rows_out:
        for k in r.keys():
            if k not in fieldnames:
                fieldnames.append(k)

    with open(out_manifest, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows_out:
            w.writerow(r)

    print("[done] runs={} manifest={} plots_dir={}".format(len(rows_out), out_manifest, plot_dir))


if __name__ == "__main__":
    sys.exit(main())
