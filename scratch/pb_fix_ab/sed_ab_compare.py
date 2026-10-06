#!/usr/bin/env python
"""Fixed-M_unit spectrum A/B: legacy vs paper-literal constant-beta P_B.

Both spectra: MAD Ma+0.94_4000 RBETAwJET rh80, M_unit 1.3488e25, seed 42,
Ns 1e6, bias 0.05 -- identical except constant_beta_paper_literal 0/1.
Computes F230 (4pi + 17-degree cone renorm), L_X(2-10 keV), L_bol, Compton
fraction; writes a summary CSV and the overlay figure.

Run: /work/vmo703/ipole_venv/bin/python sed_ab_compare.py
"""

import csv
import importlib.util

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


vcp = load_module("vcp", "/work/vmo703/igrmonty/tools/viewing_cone_postprocess.py")
stb = load_module("stb", "/work/vmo703/igrmonty/tools/sigma_thetae_breakdown.py")

BLUE, ORANGE = "#2a78d6", "#eb6834"   # paperlit, legacy (project two-series set)
SPECS = {
    "legacy": "/work/vmo703/scratch/pb_fix_ab/out/spec_mad_wjet_legacy_bhg.h5",
    "paperlit": "/work/vmo703/scratch/pb_fix_ab/out/spec_mad_wjet_paperlit.h5",
}
D_MPC = 16.8
LSUN = 3.827e33
NU_X_LO, NU_X_HI = 4.836e17, 2.418e18   # 2-10 keV
LX_GATE = 4.4e40                        # Chandra+NuSTAR 2017 core
F230_TARGET = 0.5
OUT_CSV = "/work/vmo703/scratch/pb_fix_ab/sed_ab_summary.csv"
OUT_PNG = ("/work/vmo703/igrmonty_outputs/m87/_qa/plots/zone_breakdown/"
           "pbfix_sed_ab.png")


def band_lum(nu, nulnu, lo, hi):
    m = (nu >= lo) & (nu <= hi) & (nulnu > 0)
    if m.sum() < 2:
        return float("nan")
    return float(np.trapz(nulnu[m], x=np.log(nu[m])))


def analyze(tag, path):
    spec = vcp.load_spectrum(path)
    nu = spec["nu_hz"]
    edges, centers, width = vcp.build_theta_bins(spec["dOmega_sr"].shape[0])
    sel = vcp.select_theta_bins(centers, 17.0, 10.0)
    modes = vcp.compute_spectra_modes(
        nu, spec["nuLnu_theta_freq_cgs"], spec["dOmega_sr"],
        sel["theta_mask"], D_MPC)
    with h5py.File(path, "r") as f:
        lcomp = np.array(f["/output/Lcomponent"], dtype=float) * LSUN
        l_tot = float(f["/output/L"][()])
        n_made = float(f["/output/Nmade"][()])
        n_rec = float(f["/output/Nrecorded"][()])
    nulnu4 = modes["nuLnu_4pi_cgs"]
    res = {
        "variant": tag,
        "F230_4pi_Jy": vcp.interpolate_fnu(nu, modes["fnu_4pi_jy"], 230e9, "loglog"),
        "F230_cone17_Jy": vcp.interpolate_fnu(nu, modes["fnu_cone_renorm_jy"], 230e9, "loglog"),
        "L_X_2_10keV": band_lum(nu, nulnu4, NU_X_LO, NU_X_HI),
        "L_bol_4pi": band_lum(nu, nulnu4, nu[0], nu[-1]),
        "L_total_output": l_tot,
        "compton_fraction": float(
            (lcomp[1:4].sum() + lcomp[5:8].sum()) / lcomp.sum()),
        "N_made": n_made,
        "N_recorded": n_rec,
    }
    return res, nu, nulnu4, lcomp


def main():
    rows, curves = [], {}
    for tag, path in SPECS.items():
        res, nu, nulnu4, lcomp = analyze(tag, path)
        rows.append(res)
        curves[tag] = (nu, nulnu4)
        print("[%s] Lcomponent/1e40 erg/s = %s" % (tag, np.round(lcomp / 1e40, 4)))

    L, P = rows[0], rows[1]
    print("\n%-22s %14s %14s %10s" % ("metric", "legacy", "paperlit", "leg/pl"))
    for k in ("F230_4pi_Jy", "F230_cone17_Jy", "L_X_2_10keV", "L_bol_4pi",
              "L_total_output", "compton_fraction"):
        r = L[k] / P[k] if P[k] else float("nan")
        print("%-22s %14.5g %14.5g %10.3f" % (k, L[k], P[k], r))
    print("\ngates: F230 target %.2g Jy (fixed M_unit, both un-retuned); "
          "L_X 2017 core gate %.2g erg/s" % (F230_TARGET, LX_GATE))
    print("M_unit retune factor implied (p=2): legacy x%.3f, paperlit x%.3f"
          % ((F230_TARGET / L["F230_4pi_Jy"]) ** 0.5,
             (F230_TARGET / P["F230_4pi_Jy"]) ** 0.5))

    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print("[csv] %s" % OUT_CSV)

    # ---- figure: SED overlay + ratio ------------------------------------
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.4, 5.0), constrained_layout=True,
        gridspec_kw=dict(width_ratios=[1.25, 1.0]))
    fig.patch.set_facecolor(stb.SURFACE)

    for tag, color in (("legacy", ORANGE), ("paperlit", BLUE)):
        nu, y = curves[tag]
        m = y > 0
        ax1.plot(nu[m], y[m], color=color, lw=2.0, label=tag)
    ax1.axvline(230e9, color=stb.INK2, lw=0.9, ls=(0, (4, 3)))
    ax1.text(230e9 * 1.4, 2e36, "230 GHz", color=stb.INK2, fontsize=8.5,
             rotation=90, va="bottom")
    ax1.axvspan(NU_X_LO, NU_X_HI, color=stb.INK2, alpha=0.10, lw=0)
    ax1.text(np.sqrt(NU_X_LO * NU_X_HI), 2e36, "2-10 keV", color=stb.INK2,
             fontsize=8.5, ha="center", va="bottom")
    ax1.axhline(LX_GATE, color=stb.INK2, lw=0.9, ls=(0, (1, 2)))
    ax1.text(2.2e9, LX_GATE * 0.5, r"2017 core X-ray gate $4.4{\times}10^{40}$",
             color=stb.INK2, fontsize=8.5, va="top")
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim(1e9, 3e20)
    ax1.set_ylim(1e36, 3e43)
    ax1.set_xlabel(r"$\nu$  [Hz]")
    ax1.set_ylabel(r"$\nu L_\nu$ (4$\pi$-avg)  [erg s$^{-1}$]")
    ax1.set_title("same $M_\\mathrm{unit}$, one flag flipped",
                  fontsize=10, color=stb.INK)
    ax1.legend(loc="upper right", frameon=False, fontsize=9)
    stb.style_ax(ax1)

    nuL, yL = curves["legacy"]
    nuP, yP = curves["paperlit"]
    m = (yL > 0) & (yP > 0)
    ax2.plot(nuL[m], yL[m] / yP[m], color=stb.INK, lw=2.0)
    ax2.axhline(1.0, color=stb.INK2, lw=0.9, ls=(0, (4, 3)))
    ax2.axvline(230e9, color=stb.INK2, lw=0.9, ls=(0, (4, 3)))
    ax2.axvspan(NU_X_LO, NU_X_HI, color=stb.INK2, alpha=0.10, lw=0)
    r230 = L["F230_4pi_Jy"] / P["F230_4pi_Jy"]
    rX = L["L_X_2_10keV"] / P["L_X_2_10keV"]
    ax2.annotate("x%.1f at 230 GHz" % r230, xy=(230e9, r230),
                 xytext=(1.1e10, r230 * 2.4), color=stb.INK, fontsize=9,
                 arrowprops=dict(arrowstyle="-", color=stb.INK2, lw=0.8))
    ax2.annotate("x%.0f in 2-10 keV" % rX,
                 xy=(np.sqrt(NU_X_LO * NU_X_HI), rX),
                 xytext=(3e14, rX * 2.2), color=stb.INK, fontsize=9,
                 arrowprops=dict(arrowstyle="-", color=stb.INK2, lw=0.8))
    ax2.set_xscale("log")
    ax2.set_yscale("log")
    ax2.set_xlim(1e9, 3e20)
    ax2.set_ylim(0.5, 300)  # clip the MC-noise spike at the legacy cutoff
    ax2.set_xlabel(r"$\nu$  [Hz]")
    ax2.set_ylabel(r"legacy / paper-literal")
    ax2.set_title("where the $12\\pi$ lives in the SED", fontsize=10, color=stb.INK)
    stb.style_ax(ax2)

    fig.suptitle("MAD Ma+0.94_4000 rh80 wJET -- legacy vs paper-literal $P_B$, "
                 "fixed $M_\\mathrm{unit} = 1.3488{\\times}10^{25}$",
                 fontsize=12, color=stb.INK)
    fig.savefig(OUT_PNG, dpi=150, facecolor=stb.SURFACE)
    print("[fig] %s" % OUT_PNG)


if __name__ == "__main__":
    main()
