#!/usr/bin/env python
"""Poloidal slice of the assigned wJET Thetae: does it ever drop below the cap?

Answers ashton's 2026-08-07 question directly from the dump + numpy port
(tools/sigma_thetae_breakdown.py): map the post-clamp Thetae on a fixed-phi
plane and cut theta-profiles at a few radii, marking where the field is pinned
at THETAE_HARD_MAX = 1e3 versus where it drops (R-beta disk base, jet_thetae
override). Uses the same two presets as the zone breakdown (MAD test config,
SANE production config).

Run: /work/vmo703/ipole_venv/bin/python _qa/thetae_slice_check.py
"""

import importlib.util
import os

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

SPEC = importlib.util.spec_from_file_location(
    "stb", "/work/vmo703/igrmonty/tools/sigma_thetae_breakdown.py")
stb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stb)

OUT_DIR = os.path.join(stb.QA_DIR, "plots", "zone_breakdown")
ACCENT = "#eb6834"   # cap-pinned outline / reflines (palette slot 2)
BLUES = ["#79aede", "#3a7bc4", "#0d4a8f"]  # ordered radii, light -> dark
R_MAX_PLOT = 30.0
R_PROFILES = [5.0, 10.0, 20.0]
K_SLICE = 0


def run_preset(cfg, save_npy=None):
    hdr, prims = stb.load_dump(cfg["dump"])
    cfg = dict(cfg, beta_crit=1.0)

    Mbh_cgs = cfg["MBH"] * stb.MSUN
    L_unit = stb.GNEWT * Mbh_cgs / stb.CL**2
    RHO_unit = cfg["M_unit"] / L_unit**3
    B_unit = stb.CL * np.sqrt(4.0 * np.pi * RHO_unit)
    Ne_unit = RHO_unit / (stb.MP + stb.ME)

    r, th, gcov, gcon, sqrtg = stb.fmks_geometry(hdr)
    bsq = stb.fluid_quantities(hdr, prims, gcov, gcon)
    rho, uu = prims[..., 0], prims[..., 1]
    q = stb.thetae_wjet(cfg, rho, uu, np.sqrt(bsq), gam=hdr["gam"],
                        Ne_unit=Ne_unit, B_unit=B_unit)

    a = q["assigned"]
    if save_npy is not None:
        np.save(save_npy, a)
    pinned = a >= stb.THETAE_HARD_MAX * (1.0 - 1e-9)

    # summary numbers (full 3D grid)
    sup = q["in_supplement"]
    jet = q["in_jet"]
    base = ~sup & ~jet
    print(f"\n=== {cfg['tag']} ===")
    print(f"supplement zones: {sup.sum()} ({100*sup.mean():.3f}% of grid)")
    if sup.any():
        print(f"  pinned at cap: {100*pinned[sup].mean():.4f}%   "
              f"min/max assigned in supplement: {a[sup].min():.6g} / {a[sup].max():.6g}")
    print(f"override zones (Thetae={cfg['jet_thetae']:g}): {jet.sum()} "
          f"({100*jet.mean():.3f}%)  assigned min/max: {a[jet].min():.4g}/{a[jet].max():.4g}")
    print(f"base zones: {100*base.mean():.2f}%  assigned median/p95: "
          f"{np.median(a[base]):.4g} / {np.percentile(a[base], 95):.4g}"
          f"   pinned: {100*pinned[base].mean():.4f}%")
    print(f"whole grid pinned at cap: {100*pinned.mean():.3f}%")

    # ---- figure: poloidal slice + theta profiles --------------------------
    k = K_SLICE
    ak = a[:, :, k]
    sig = q["sigma"][:, :, k]
    pin_k = pinned[:, :, k]
    x = r * np.sin(th)
    z = r * np.cos(th)

    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(12.6, 5.4), constrained_layout=True,
        gridspec_kw=dict(width_ratios=[1.05, 1.0]))
    fig.patch.set_facecolor(stb.SURFACE)

    pc = ax1.pcolormesh(x, z, ak,
                        norm=LogNorm(vmin=stb.THETAE_MIN, vmax=stb.THETAE_HARD_MAX),
                        cmap="Blues", shading="gouraud", rasterized=True)
    ax1.contour(x, z, sig, levels=[2.0], colors=stb.INK2, linewidths=1.0,
                linestyles="dashed")
    ax1.contour(x, z, sig, levels=[10.0], colors=stb.INK2, linewidths=1.0)
    # the pinned skin is 1-2 cells thick and patchy in phi; dots render where a
    # contour of the mask would vanish. Faint = pinned at ANY phi (skin
    # footprint), solid = pinned in THIS slice.
    pin_any = pinned.any(axis=2)
    if pin_any.any():
        ax1.scatter(x[pin_any], z[pin_any], s=4, c=ACCENT, alpha=0.22,
                    marker="o", linewidths=0, zorder=4)
    if pin_k.any():
        ax1.scatter(x[pin_k], z[pin_k], s=8, c=ACCENT, marker="o",
                    linewidths=0, zorder=5)
    ax1.set_xlim(0, R_MAX_PLOT)
    ax1.set_ylim(-R_MAX_PLOT, R_MAX_PLOT)
    ax1.set_aspect("equal")
    ax1.set_xlabel(r"$x = r\,\sin\theta$  [$r_g$]")
    ax1.set_ylabel(r"$z = r\,\cos\theta$  [$r_g$]")
    ax1.set_title(r"assigned $\Theta_e$, $\phi$-slice $k=0$"
                  "\n" r"orange dots pinned at cap $10^3$ (solid: this slice,"
                  r" faint: any $\phi$); gray: $\sigma{=}2$/$10$",
                  fontsize=10, color=stb.INK)
    cb = fig.colorbar(pc, ax=ax1, shrink=0.9, pad=0.02)
    cb.set_label(r"$\Theta_e$")
    stb.style_ax(ax1)

    rcol = r[:, 0]
    th_deg_all = []
    for rt, c in zip(R_PROFILES, BLUES):
        i = int(np.argmin(np.abs(rcol - rt)))
        th_deg = np.degrees(th[i, :])
        th_deg_all.append(th_deg)
        ax2.plot(th_deg, a[i, :, k], color=c, lw=2.0,
                 label=rf"$r \approx {rcol[i]:.1f}\,r_g$")
    ax2.axhline(stb.THETAE_HARD_MAX, color=ACCENT, lw=1.2, ls=(0, (4, 3)))
    ax2.text(91, stb.THETAE_HARD_MAX * 1.25, r"cap $10^3$", color=ACCENT,
             fontsize=9, ha="center")
    ax2.axhline(cfg["jet_thetae"], color=stb.INK2, lw=1.0, ls=(0, (4, 3)))
    ax2.text(91, cfg["jet_thetae"] * 1.25, rf"override {cfg['jet_thetae']:g}",
             color=stb.INK2, fontsize=9, ha="center")
    ax2.set_yscale("log")
    ax2.set_ylim(3e-4, 8e3)
    ax2.set_xlim(0, 180)
    ax2.set_xticks([0, 45, 90, 135, 180])
    ax2.set_xlabel(r"$\theta$  [deg]  (0 = north pole, 90 = midplane)")
    ax2.set_ylabel(r"assigned $\Theta_e$")
    ax2.set_title(r"$\Theta_e(\theta)$ cuts through the same slice",
                  fontsize=10, color=stb.INK)
    ax2.legend(loc="lower center", frameon=False, fontsize=9)
    stb.style_ax(ax2)

    fig.suptitle(f"{cfg['tag']} — where the assigned $\\Theta_e$ sits vs the cap",
                 fontsize=12, color=stb.INK)
    out = os.path.join(OUT_DIR, f"thetae_slice_{cfg['tag'].replace('+','p').replace('.','_')}.png")
    fig.savefig(out, dpi=150, facecolor=stb.SURFACE)
    plt.close(fig)
    print(f"[fig] {out}")
    return out


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    for cfg in stb.PRESETS:
        run_preset(cfg)


if __name__ == "__main__":
    main()
