#!/usr/bin/env python3
"""Figures for the MWL SED check (see mwl_sed_check.py for provenance).

Outputs under _qa/plots/mwl_sed/:
  fig1_lx_gate.png         L_X(2-10 keV, 17-deg cone) / observed, per family
  fig2_f230_frame_gap.png  F230 in the 17-deg cone vs the 0.5 Jy tuning target
  fig3_sed_families.png    cone-frame SEDs per family with observed anchors

Run:  /work/vmo703/ipole_venv/bin/python mwl_sed_figures.py
"""

import argparse
import importlib.util
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

BASE_DIR = "/work/vmo703/igrmonty_outputs/m87"
QA_DIR = os.path.join(BASE_DIR, "_qa")
OUT_DIR = os.path.join(QA_DIR, "plots", "mwl_sed")
TOOL_PATH = "/work/vmo703/igrmonty/tools/viewing_cone_postprocess.py"

# dataviz reference palette (light mode): slots 1-2 + text/surface tokens
C_POS0 = "#2a78d6"
C_POS1 = "#eb6834"
INK = "#0b0b0b"
INK2 = "#52514e"
SURFACE = "#fcfcfb"
GRID = dict(color=INK2, alpha=0.16, lw=0.7)

LX_OBS = 4.4e40
F230_TARGET = 0.5
D_CM = 16.8e6 * 3.085677581e18
FOURPI_D2 = 4.0 * np.pi * D_CM**2
NULNU_230_OBS = 0.5e-23 * 230.0e9 * FOURPI_D2  # 0.5 Jy at 230 GHz -> nuLnu
XLO, XHI = 2 * 2.417989242e17, 10 * 2.417989242e17
NULNU_X_OBS = LX_OBS / np.log(XHI / XLO)

FAMILIES = [
    ("MAD", "CRITBETA", "OFF", "MAD Crit-β (rh20 bc0.01 f0.5)"),
    ("MAD", "RBETA", "OFF", "MAD R-β (rh20)"),
    # Phase-4 additions (2026-09): the jet-supplement MAD families
    ("MAD", "CRITBETA", "ON", "MAD Crit-β wJET (rh20 bc0.01 f0.5)"),
    ("MAD", "RBETA", "ON", "MAD R-β wJET (rh80)"),
    ("SANE", "CRITBETA", "ON", "SANE Crit-β wJET (rh20 bc1 f0.5)"),
    ("SANE", "RBETA", "ON", "SANE R-β wJET (rh160)"),
]


def style_ax(ax):
    ax.set_facecolor(SURFACE)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(INK2)
        ax.spines[side].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.grid(True, which="major", axis="x", **GRID)
    ax.set_axisbelow(True)


def strip_plot(df, value_col, xlabel, refline, refline_label, out_png,
               shade_above_ref=False, annotate_max=False, extra_pts=None,
               xlim=None, title="", subtitle="", legend_loc="lower right"):
    fig, ax = plt.subplots(figsize=(8.4, 1.05 * len(FAMILIES) + 0.8), facecolor=SURFACE)
    style_ax(ax)

    rng = np.random.default_rng(11)
    for i, (state, heating, wjet, label) in enumerate(FAMILIES):
        sub = df[(df.state == state) & (df.heating == heating) & (df.wjet == wjet)]
        for pos, color, dy in ((0, C_POS0, 0.16), (1, C_POS1, -0.16)):
            vals = sub[sub.pos == pos][value_col].to_numpy()
            jit = rng.uniform(-0.05, 0.05, size=vals.size)
            ax.scatter(
                vals, np.full(vals.size, i + dy) + jit,
                s=64, color=color, edgecolors=SURFACE, linewidths=1.4, zorder=3,
            )
    if extra_pts is not None:
        for x, y, note in extra_pts:
            ax.scatter([x], [y], marker="x", s=56, color=INK2, zorder=3, linewidths=1.6)
            ax.annotate(note, (x, y), textcoords="offset points", xytext=(-6, 9),
                        fontsize=8, color=INK2, ha="right")

    ax.axvline(refline, color=INK, lw=1.1, ls=(0, (4, 3)), zorder=2)
    ax.text(refline * 1.08, -0.42, refline_label, fontsize=8.5, color=INK)
    if shade_above_ref:
        ax.axvspan(refline, ax.get_xlim()[1] if xlim is None else xlim[1],
                   color=INK2, alpha=0.07, zorder=1)

    if annotate_max:
        idx = df[value_col].idxmax()
        row = df.loc[idx]
        fam_i = next(i for i, f in enumerate(FAMILIES)
                     if f[0] == row.state and f[1] == row.heating and f[2] == row.wjet)
        ax.annotate(
            f"worst headroom: {row[value_col]:.2f}x",
            (row[value_col], fam_i - 0.16),
            textcoords="offset points", xytext=(10, 5), fontsize=8.5, color=INK,
        )

    ax.set_xscale("log")
    if xlim:
        ax.set_xlim(*xlim)
    ax.set_yticks(range(len(FAMILIES)))
    ax.set_yticklabels([f[3] for f in FAMILIES], fontsize=9.5, color=INK)
    ax.set_ylim(-0.6, len(FAMILIES) - 0.15)
    ax.invert_yaxis()
    ax.set_xlabel(xlabel, fontsize=10, color=INK)

    handles = [
        plt.Line2D([], [], marker="o", ls="", color=C_POS0, markeredgecolor=SURFACE, label="pos0 (no pairs)"),
        plt.Line2D([], [], marker="o", ls="", color=C_POS1, markeredgecolor=SURFACE, label="pos1 (pair-loaded)"),
    ]
    ax.legend(handles=handles, loc=legend_loc, frameon=False, fontsize=8.5,
              labelcolor=INK2)
    ax.set_title(title, fontsize=10.5, color=INK, loc="left", pad=24)
    ax.text(0, 1.012, subtitle, transform=ax.transAxes, fontsize=8.5, color=INK2)
    fig.tight_layout()
    fig.savefig(out_png, dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print("[fig]", out_png)


def main():
    ap = argparse.ArgumentParser(description="Figures for the MWL SED check")
    ap.add_argument(
        "--dir",
        default=os.path.join(BASE_DIR, "pre_bugfix_2026-07"),
        help="Corpus directory matching the CSVs from mwl_sed_check.py",
    )
    args = ap.parse_args()
    corpus_dir = args.dir

    os.makedirs(OUT_DIR, exist_ok=True)
    v = pd.read_csv(os.path.join(QA_DIR, "mwl_sed_verdicts.csv"))
    canon = v[v.trial.isna()].copy()

    trial = v[v.filename.str.contains("trial01") & (v.state == "SANE")]
    extra = []
    if len(trial):
        extra.append((float(trial.iloc[0].Lx_cone_over_obs), 3 - 0.16,
                      "stale trial01 (untuned M_unit)"))

    # data-driven title: never hardcode the run count or the verdict
    n_runs = len(canon)
    n_over = int((canon.Lx_cone_over_obs > 1.0).sum())
    if n_over:
        fig1_title = ("X-ray gate: %d of %d tuned runs exceed the observed "
                      "core in the 17° frame" % (n_over, n_runs))
    else:
        fig1_title = ("X-ray gate: all %d tuned runs sit below the observed "
                      "core" % n_runs)
    strip_plot(
        canon, "Lx_cone_over_obs",
        "L_X(2–10 keV) toward 17° ÷ observed core (4.4×10⁴⁰ erg s⁻¹)",
        1.0, "observed", os.path.join(OUT_DIR, "fig1_lx_gate.png"),
        shade_above_ref=True, annotate_max=True, extra_pts=extra,
        xlim=(2e-4, 30),
        title=fig1_title,
        subtitle="2017 EHT MWL campaign anchor · cone 17°±10° (folded)",
    )

    strip_plot(
        canon, "f230_cone_jy",
        "Fν(230 GHz) toward 17°  [Jy]",
        F230_TARGET, "0.5 Jy target", os.path.join(OUT_DIR, "fig2_f230_frame_gap.png"),
        xlim=(0.012, 1.3),
        title="Tuning-frame gap: the 17° observer sees less than the 4π-average target",
        subtitle="M_unit tuned on the 4π-average; cone/4π gap computed per run "
                 "(era: %d runs)" % n_runs,
        legend_loc="lower left",
    )

    # fig 3: per-family cone SEDs with observed anchors
    spec = importlib.util.spec_from_file_location("vcp", TOOL_PATH)
    vcp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(vcp)

    fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.2), facecolor=SURFACE,
                             sharex=True, sharey=True)
    for ax, (state, heating, wjet, label) in zip(axes.ravel(), FAMILIES):
        style_ax(ax)
        ax.grid(True, which="major", axis="y", **GRID)
        sub = canon[(canon.state == state) & (canon.heating == heating) & (canon.wjet == wjet)]
        for _, row in sub.iterrows():
            path = os.path.join(corpus_dir, row.filename)
            data = vcp.load_spectrum(path)
            nu = data["nu_hz"]
            nb = int(data["dOmega_sr"].shape[0])
            _, centers, _ = vcp.build_theta_bins(nb)
            tm = vcp.select_theta_bins(centers, 17.0, 10.0)
            modes = vcp.compute_spectra_modes(
                nu, data["nuLnu_theta_freq_cgs"], data["dOmega_sr"],
                tm["theta_mask"], 16.8,
            )
            y = modes["nuLnu_cone_renorm_cgs"]
            ok = (y > 0) & (nu > 1e9) & (nu < 3e21)
            ax.plot(nu[ok], y[ok], lw=1.3, alpha=0.55,
                    color=C_POS0 if row.pos == 0 else C_POS1, zorder=2)

        ax.scatter([230e9], [NULNU_230_OBS], marker="D", s=46, color=INK,
                   zorder=4, edgecolors=SURFACE, linewidths=1.2)
        ax.plot([XLO, XHI], [NULNU_X_OBS, NULNU_X_OBS], color=INK, lw=2.2,
                zorder=4, solid_capstyle="butt")
        ax.set_title(label, fontsize=10, color=INK, loc="left")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(1e9, 3e21)
        ax.set_ylim(1e35, 2e42)

    axes[0, 0].annotate("EHT-scale 0.5 Jy", (230e9, NULNU_230_OBS),
                        textcoords="offset points", xytext=(6, 8), fontsize=8, color=INK)
    axes[0, 0].annotate("Chandra+NuSTAR core", (XLO, NULNU_X_OBS),
                        textcoords="offset points", xytext=(-2, 8), fontsize=8, color=INK)
    for ax in axes[1, :]:
        ax.set_xlabel("ν  [Hz]", fontsize=10, color=INK)
    for ax in axes[:, 0]:
        ax.set_ylabel("ν Lν toward 17°  [erg s⁻¹]", fontsize=10, color=INK)
    handles = [
        plt.Line2D([], [], color=C_POS0, lw=2, label="pos0 (no pairs)"),
        plt.Line2D([], [], color=C_POS1, lw=2, label="pos1 (pair-loaded)"),
        plt.Line2D([], [], color=INK, lw=0, marker="D", label="observed anchors"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=9, labelcolor=INK2, bbox_to_anchor=(0.5, -0.005))
    fig.suptitle("Cone-frame SEDs (17°±10°) vs 2017 M87 core anchors — 6 dumps × 2 spins per panel",
                 fontsize=12, color=INK, x=0.02, ha="left")
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    fig.savefig(os.path.join(OUT_DIR, "fig3_sed_families.png"), dpi=170, facecolor=SURFACE)
    plt.close(fig)
    print("[fig]", os.path.join(OUT_DIR, "fig3_sed_families.png"))


if __name__ == "__main__":
    main()
