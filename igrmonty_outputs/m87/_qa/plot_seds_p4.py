#!/usr/bin/env python
"""Per-branch SED plots for every Phase-4 final spectrum.

One PNG per final in run_2026-09-15 -> plots/seds_p4/<branch>.png:
4pi-averaged and 17-deg-cone nuLnu, the 2-10 keV band, the 2017 core X-ray
gate, and the 230 GHz anchor. Same palette/axes style as the campaign
figures (pbfix_sed_ab et al.).

Run: /work/vmo703/ipole_venv/bin/python plot_seds_p4.py
"""

import importlib.util
import os

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
                       "plots", "seds_p4")
BLUE, ORANGE = "#2a78d6", "#eb6834"
NU_X_LO, NU_X_HI = 4.836e17, 2.418e18
LX_GATE = 4.4e40
D_MPC = 16.8


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    finals = sorted(f for f in os.listdir(RUN_DIR)
                    if f.startswith("spectrum_") and f.endswith(".h5")
                    and "_trial" not in f)
    for name in finals:
        spec = vcp.load_spectrum(os.path.join(RUN_DIR, name))
        nu = spec["nu_hz"]
        _, centers, _ = vcp.build_theta_bins(spec["dOmega_sr"].shape[0])
        sel = vcp.select_theta_bins(centers, 17.0, 10.0)
        modes = vcp.compute_spectra_modes(nu, spec["nuLnu_theta_freq_cgs"],
                                          spec["dOmega_sr"],
                                          sel["theta_mask"], D_MPC)
        y4 = modes["nuLnu_4pi_cgs"]
        yc = modes["nuLnu_cone_renorm_cgs"]
        f230 = vcp.interpolate_fnu(nu, modes["fnu_4pi_jy"], 230e9, "loglog")

        branch = name.replace("spectrum_", "").replace(".h5", "")
        fig, ax = plt.subplots(figsize=(7.4, 5.2))
        fig.patch.set_facecolor(stb.SURFACE)
        m = y4 > 0
        ax.plot(nu[m], y4[m], color=BLUE, lw=2.0, label=r"4$\pi$-averaged")
        m = yc > 0
        ax.plot(nu[m], yc[m], color=ORANGE, lw=1.6,
                label=r"17$^\circ$ cone (renorm)")
        ax.axvspan(NU_X_LO, NU_X_HI, color=stb.INK2, alpha=0.10, lw=0)
        ax.axhline(LX_GATE, color=stb.INK2, lw=0.9, ls=(0, (1, 2)))
        ax.text(2.2e9, LX_GATE * 0.45, r"2017 X-ray gate $4.4{\times}10^{40}$",
                color=stb.INK2, fontsize=8, va="top")
        ax.axvline(230e9, color=stb.INK2, lw=0.9, ls=(0, (4, 3)))
        ax.text(230e9 * 1.4, 3e36, "230 GHz\n%.2f Jy" % f230, color=stb.INK2,
                fontsize=8, rotation=90, va="bottom")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e9, 3e20)
        ax.set_ylim(1e36, 3e43)
        ax.set_xlabel(r"$\nu$  [Hz]")
        ax.set_ylabel(r"$\nu L_\nu$  [erg s$^{-1}$]")
        ax.set_title(branch, fontsize=10, color=stb.INK)
        ax.legend(loc="upper right", frameon=False, fontsize=9)
        stb.style_ax(ax)
        fig.tight_layout()
        fig.savefig(os.path.join(OUT_DIR, branch + ".png"), dpi=140,
                    facecolor=stb.SURFACE)
        plt.close(fig)
    print("%d SEDs -> %s" % (len(finals), OUT_DIR))


if __name__ == "__main__":
    main()
