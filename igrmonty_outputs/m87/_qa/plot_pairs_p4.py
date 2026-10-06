#!/usr/bin/env python
"""Phase-4 pair comparisons: pos0 vs pos1 SED overlay per branch pair.

One PNG per (family, spin, dump) -> plots/pairs_p4/: both positron states'
4pi-averaged SEDs, with the pair-driven shifts annotated (F230, L_X, and the
L_X pos1/pos0 ratio -- the campaign-wide pair signature is ~2x in X-ray).
Successor to the February pair_sweep 'pairs/' plots for the production era.

Run: /work/vmo703/ipole_venv/bin/python plot_pairs_p4.py
"""

import importlib.util
import os
import re

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

RUN_DIR = "/work/vmo703/igrmonty_outputs/m87/run_2026-09-15"
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       "plots", "pairs_p4")
BLUE, ORANGE = "#2a78d6", "#eb6834"   # pos0, pos1 (campaign convention)
NU_X_LO, NU_X_HI = 4.836e17, 2.418e18
LX_GATE = 4.4e40
D_MPC = 16.8


def metrics(path):
    s = vcp.load_spectrum(path)
    nu = s["nu_hz"]
    y4 = (s["nuLnu_theta_freq_cgs"] * s["dOmega_sr"].reshape(-1, 1)).sum(axis=0) / (4 * np.pi)
    m = (nu >= NU_X_LO) & (nu <= NU_X_HI) & (y4 > 0)
    lx = float(np.trapz(y4[m], x=np.log(nu[m])))
    modes = vcp.compute_spectra_modes(
        nu, s["nuLnu_theta_freq_cgs"], s["dOmega_sr"],
        np.ones(s["dOmega_sr"].shape[0], bool), D_MPC)
    f230 = vcp.interpolate_fnu(nu, modes["fnu_4pi_jy"], 230e9, "loglog")
    return nu, y4, lx, f230


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    finals = sorted(f for f in os.listdir(RUN_DIR)
                    if f.startswith("spectrum_") and f.endswith(".h5")
                    and "_trial" not in f)
    stems = {}
    for f in finals:
        m = re.match(r"(spectrum_.+)_pos([01])\.h5$", f)
        if m:
            stems.setdefault(m.group(1), {})[m.group(2)] = f
    n = 0
    for stem in sorted(stems):
        pair = stems[stem]
        if "0" not in pair or "1" not in pair:
            print("[warn] incomplete pair for %s" % stem)
            continue
        nu0, y0, lx0, f0 = metrics(os.path.join(RUN_DIR, pair["0"]))
        nu1, y1, lx1, f1 = metrics(os.path.join(RUN_DIR, pair["1"]))
        branch = stem.replace("spectrum_", "")

        fig, ax = plt.subplots(figsize=(7.4, 5.2))
        fig.patch.set_facecolor(stb.SURFACE)
        m0, m1 = y0 > 0, y1 > 0
        ax.plot(nu0[m0], y0[m0], color=BLUE, lw=2.0, label="pos0 (no pairs)")
        ax.plot(nu1[m1], y1[m1], color=ORANGE, lw=2.0, label="pos1 (pair-loaded)")
        ax.axvspan(NU_X_LO, NU_X_HI, color=stb.INK2, alpha=0.10, lw=0)
        ax.axhline(LX_GATE, color=stb.INK2, lw=0.9, ls=(0, (1, 2)))
        ax.axvline(230e9, color=stb.INK2, lw=0.9, ls=(0, (4, 3)))
        ax.text(0.02, 0.05,
                "F230: %.2f / %.2f Jy\nL_X: %.2g / %.2g erg/s\n"
                "pair X-ray boost: x%.2f"
                % (f0, f1, lx0, lx1, lx1 / lx0 if lx0 else float("nan")),
                transform=ax.transAxes, fontsize=9, color=stb.INK,
                va="bottom")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e9, 3e20)
        ax.set_ylim(1e36, 3e43)
        ax.set_xlabel(r"$\nu$  [Hz]")
        ax.set_ylabel(r"$\nu L_\nu$ (4$\pi$-avg)  [erg s$^{-1}$]")
        ax.set_title(branch + "  --  positron pair comparison",
                     fontsize=10, color=stb.INK)
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        stb.style_ax(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, "pair_%s.png" % branch), dpi=140,
                    facecolor=stb.SURFACE)
        plt.close(fig)
        n += 1
    print("%d pair plots -> %s" % (n, OUT_DIR))


if __name__ == "__main__":
    main()
