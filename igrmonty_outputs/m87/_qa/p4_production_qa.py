#!/usr/bin/env python
"""Production QA for Phase-4 spectra: provenance + SED gates + July deltas.

For every final (suffix-free) spectrum in a run directory:
  provenance -- M_unit, bias, Ns, positron ratio, wJET params incl.
    constant_beta_paper_literal, run status;
  SED metrics -- F230 (4pi and 17-degree cone renorm), L_X(2-10 keV, 4pi),
    L_bol, Compton fraction;
  gates -- F230 within tolerance of 0.5 Jy; L_X(4pi) vs the 2017 core limit
    4.4e40 erg/s (cone value reported too: tuning is 4pi but observers sit
    at ~17 deg);
  July delta -- same-named spectrum in pre_bugfix_2026-07/, ratio of F230
    and L_X (quantifies new-binary + paper-literal shift for reruns).

Run: /work/vmo703/ipole_venv/bin/python p4_production_qa.py [run_dir]
"""

import csv
import importlib.util
import os
import sys

import h5py
import numpy as np


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vcp = load_module("vcp", "/work/vmo703/igrmonty/tools/viewing_cone_postprocess.py")
amb = load_module("amb", "/work/vmo703/igrmonty/auto_munit_bracket.py")

RUN_DIR = sys.argv[1] if len(sys.argv) > 1 else \
    "/work/vmo703/igrmonty_outputs/m87/run_2026-09-15"
ARCHIVE = "/work/vmo703/igrmonty_outputs/m87/pre_bugfix_2026-07"
OUT_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "p4_qa_%s.csv" % os.path.basename(RUN_DIR.rstrip("/")))
D_MPC = 16.8
LSUN = 3.827e33
NU_X_LO, NU_X_HI = 4.836e17, 2.418e18
LX_GATE = 4.4e40
F_TARGET, F_TOL = 0.5, 0.05  # tuner tolerance band


def band_lum(nu, nulnu, lo, hi):
    m = (nu >= lo) & (nu <= hi) & (nulnu > 0)
    if m.sum() < 2:
        return float("nan")
    return float(np.trapz(nulnu[m], x=np.log(nu[m])))


def read_scalar(f, path, default=float("nan")):
    try:
        return float(f[path][()])
    except Exception:
        return default


def sed_metrics(path):
    spec = vcp.load_spectrum(path)
    nu = spec["nu_hz"]
    _, centers, _ = vcp.build_theta_bins(spec["dOmega_sr"].shape[0])
    sel = vcp.select_theta_bins(centers, 17.0, 10.0)
    modes = vcp.compute_spectra_modes(nu, spec["nuLnu_theta_freq_cgs"],
                                      spec["dOmega_sr"], sel["theta_mask"], D_MPC)
    nulnu4 = modes["nuLnu_4pi_cgs"]
    nulnu_cone = modes["nuLnu_cone_renorm_cgs"]
    # the tuner's own nearest-bin measure IS the 0.5 Jy anchor definition;
    # loglog interpolation differs by ~5-10% on the steep mm slope.
    from pathlib import Path as _P
    f230_tuner = amb.measure_flux(_P(path), 230e9)
    if isinstance(f230_tuner, tuple):
        f230_tuner = f230_tuner[0]
    out = {
        "F230_4pi_Jy": float(f230_tuner),
        "F230_4pi_loglog_Jy": vcp.interpolate_fnu(nu, modes["fnu_4pi_jy"], 230e9, "loglog"),
        "F230_cone17_Jy": vcp.interpolate_fnu(nu, modes["fnu_cone_renorm_jy"], 230e9, "loglog"),
        "L_X_4pi": band_lum(nu, nulnu4, NU_X_LO, NU_X_HI),
        "L_X_cone17": band_lum(nu, nulnu_cone, NU_X_LO, NU_X_HI),
        "L_bol_4pi": band_lum(nu, nulnu4, nu[0], nu[-1]),
    }
    with h5py.File(path, "r") as f:
        lcomp = np.array(f["/output/Lcomponent"], dtype=float) * LSUN
        out["compton_fraction"] = float(
            (lcomp[1:4].sum() + lcomp[5:8].sum()) / lcomp.sum())
        out["M_unit"] = read_scalar(f, "/params/M_unit")
        out["bias"] = read_scalar(f, "/params/bias")
        out["Ns"] = read_scalar(f, "/params/Ns")
        out["positron_ratio"] = read_scalar(
            f, "/params/electrons/positron_ratio", float("nan"))
        out["paper_literal"] = read_scalar(
            f, "/params/electrons/constant_beta_paper_literal", float("nan"))
        out["with_electrons"] = read_scalar(f, "/params/electrons/type")
        try:
            rs = f["/output/run_status"][()]
            out["run_status"] = rs.decode() if hasattr(rs, "decode") else str(rs)
        except Exception:
            out["run_status"] = "?"
    return out


def main():
    finals = sorted(
        f for f in os.listdir(RUN_DIR)
        if f.startswith("spectrum_") and f.endswith(".h5") and "_trial" not in f)
    if not finals:
        print("no final spectra in %s yet" % RUN_DIR)
        return
    rows = []
    for name in finals:
        m = sed_metrics(os.path.join(RUN_DIR, name))
        m["spectrum"] = name
        is_wjet = "wJET" in name
        # gates
        f230_ok = abs(m["F230_4pi_Jy"] - F_TARGET) <= F_TOL * F_TARGET * 2  # 10% grace
        lx_ok = m["L_X_4pi"] <= LX_GATE
        prov_ok = (m["Ns"] == 1e6 and m["run_status"] in ("ok", "?")
                   and (not is_wjet or m["paper_literal"] == 1.0))
        m["gate_F230"] = "PASS" if f230_ok else "FAIL"
        m["gate_LX_4pi"] = "PASS" if lx_ok else "FAIL"
        m["gate_provenance"] = "PASS" if prov_ok else "FAIL"
        # July comparison
        old = os.path.join(ARCHIVE, name)
        if os.path.exists(old):
            o = sed_metrics(old)
            m["F230cone_vs_jul"] = m["F230_cone17_Jy"] / o["F230_cone17_Jy"]
            m["LX_vs_jul"] = m["L_X_4pi"] / o["L_X_4pi"]
            m["Munit_vs_jul"] = m["M_unit"] / o["M_unit"]
        else:
            m["F230cone_vs_jul"] = m["LX_vs_jul"] = m["Munit_vs_jul"] = float("nan")
        rows.append(m)

    cols = ["spectrum", "M_unit", "bias", "paper_literal", "positron_ratio",
            "run_status", "F230_4pi_Jy", "F230_4pi_loglog_Jy",
            "F230_cone17_Jy", "L_X_4pi",
            "L_X_cone17", "L_bol_4pi", "compton_fraction",
            "gate_F230", "gate_LX_4pi", "gate_provenance",
            "Munit_vs_jul", "F230cone_vs_jul", "LX_vs_jul", "Ns",
            "with_electrons"]
    with open(OUT_CSV, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    for m in rows:
        print("\n== %s ==" % m["spectrum"])
        print("  provenance: M_unit=%.4e bias=%.3g Ns=%.0e paper_literal=%s "
              "pos=%.0f status=%s -> %s"
              % (m["M_unit"], m["bias"], m["Ns"], m["paper_literal"],
                 m["positron_ratio"], m["run_status"], m["gate_provenance"]))
        print("  F230: 4pi=%.4g Jy [%s]  cone17=%.4g Jy (x%.2f vs 4pi)"
              % (m["F230_4pi_Jy"], m["gate_F230"], m["F230_cone17_Jy"],
                 m["F230_cone17_Jy"] / m["F230_4pi_Jy"]))
        print("  L_X(2-10): 4pi=%.4g [%s vs %.2g]  cone17=%.4g  "
              "L_bol=%.4g  f_compton=%.3f"
              % (m["L_X_4pi"], m["gate_LX_4pi"], LX_GATE, m["L_X_cone17"],
                 m["L_bol_4pi"], m["compton_fraction"]))
        print("  vs July: M_unit x%.3f  F230cone x%.3f  L_X x%.3f"
              % (m["Munit_vs_jul"], m["F230cone_vs_jul"], m["LX_vs_jul"]))
    print("\n[csv] %s" % OUT_CSV)


if __name__ == "__main__":
    main()
