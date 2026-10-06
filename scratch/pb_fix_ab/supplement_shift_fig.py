#!/usr/bin/env python
"""The money plot for the paper-literal P_B fix: what changes in the sheath.

Left: assigned Thetae vs sigma in the supplement band (sigma in [2,10), beta >
0.1) for legacy vs paper-literal, with the 92-sigma supplement-only reference.
Right: distribution of supplement-zone Thetae -- legacy is a delta at the cap;
paper-literal spreads over ~270-1000 with an 11% remnant at the cap.
Also prints the emissivity-proxy (ne Thetae^2 B^2 sqrt(g)) shift at fixed
M_unit, the first-order expectation for the spectrum A/B.
"""

import importlib.util

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPEC = importlib.util.spec_from_file_location(
    "stb", "/work/vmo703/igrmonty/tools/sigma_thetae_breakdown.py")
stb = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stb)

BLUE, ORANGE = "#2a78d6", "#eb6834"
OUT_PNG = "/work/vmo703/igrmonty_outputs/m87/_qa/plots/zone_breakdown/pbfix_supplement_shift.png"

cfg = dict(stb.PRESETS[0], beta_crit=1.0)
hdr, prims = stb.load_dump(cfg["dump"])
Mbh_cgs = cfg["MBH"] * stb.MSUN
L_unit = stb.GNEWT * Mbh_cgs / stb.CL**2
RHO_unit = cfg["M_unit"] / L_unit**3
B_unit = stb.CL * np.sqrt(4.0 * np.pi * RHO_unit)
Ne_unit = RHO_unit / (stb.MP + stb.ME)

r, th, gcov, gcon, sqrtg = stb.fmks_geometry(hdr)
bsq = stb.fluid_quantities(hdr, prims, gcov, gcon)
rho, uu, b_code = prims[..., 0], prims[..., 1], np.sqrt(bsq)

cfg_pl = dict(cfg, constant_beta_e0=cfg["constant_beta_e0"] / (12.0 * np.pi))
qL = stb.thetae_wjet(cfg, rho, uu, b_code, gam=hdr["gam"], Ne_unit=Ne_unit, B_unit=B_unit)
qP = stb.thetae_wjet(cfg_pl, rho, uu, b_code, gam=hdr["gam"], Ne_unit=Ne_unit, B_unit=B_unit)

sup = qL["in_supplement"]  # identical mask in both variants (sigma/beta unchanged)
ne = np.clip(rho, 1e-30, None) * Ne_unit
Bc = b_code * B_unit
w = ne * Bc**2 * sqrtg[:, :, None]
for name, q in (("legacy", qL), ("paper-literal", qP)):
    w_em = w * q["assigned"] ** 2
    print(f"{name:14s} emissivity proxy: supplement = {w_em[sup].sum():.4e}  "
          f"total = {w_em.sum():.4e}  supplement share = {w_em[sup].sum()/w_em.sum():.4f}")
print(f"supplement-band proxy ratio legacy/paper-literal = "
      f"{(w*qL['assigned']**2)[sup].sum() / (w*qP['assigned']**2)[sup].sum():.2f}x")
print(f"total proxy ratio legacy/paper-literal = "
      f"{(w*qL['assigned']**2).sum() / (w*qP['assigned']**2).sum():.2f}x")

sig_s = qL["sigma"][sup]
aL, aP = qL["assigned"][sup], qP["assigned"][sup]

rng = np.random.RandomState(42)
idx = rng.choice(sig_s.size, size=min(6000, sig_s.size), replace=False)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.2, 5.0), constrained_layout=True)
fig.patch.set_facecolor(stb.SURFACE)

ax1.scatter(sig_s[idx], aL[idx], s=5, c=ORANGE, linewidths=0, alpha=0.5,
            label="legacy (as-shipped)")
ax1.scatter(sig_s[idx], aP[idx], s=5, c=BLUE, linewidths=0, alpha=0.35,
            label="paper-literal $B^2/8\\pi$")
ss = np.linspace(2, 10, 100)
ax1.plot(ss, 91.8 * ss, color=stb.INK2, lw=1.2, ls=(0, (4, 3)))
ax1.text(7.6, 91.8 * 7.6 * 0.82, r"$91.8\,\sigma$ (supplement term alone)",
         color=stb.INK2, fontsize=8.5, rotation=13)
ax1.axhline(1000, color=stb.INK2, lw=0.8)
ax1.set_xlim(2, 10)
ax1.set_ylim(0, 1060)
ax1.set_xlabel(r"$\sigma = b^2/\rho$")
ax1.set_ylabel(r"assigned $\Theta_e$")
ax1.set_title("supplement band: flat at the cap vs real $\\sigma$-dependence",
              fontsize=10, color=stb.INK)
ax1.legend(loc="center right", frameon=False, fontsize=9)
stb.style_ax(ax1)

bins = np.arange(250, 1026, 25)
ax2.hist(aP, bins=bins, color=BLUE, alpha=0.85, label="paper-literal")
ax2.axvline(1000, color=ORANGE, lw=2.0)
nL = aL.size
ax2.text(0.86, 0.86,
         f"legacy: {nL:,}/{nL:,} zones\nat exactly 1000",
         transform=ax2.transAxes, ha="right", color=ORANGE, fontsize=9)
pinnedP = (aP >= 1000 * (1 - 1e-9)).sum()
ax2.text(0.86, 0.70,
         f"paper-literal: {pinnedP:,} zones\n({100*pinnedP/nL:.1f}%) still at cap",
         transform=ax2.transAxes, ha="right", color=BLUE, fontsize=9)
ax2.set_xlabel(r"assigned $\Theta_e$ in supplement zones")
ax2.set_ylabel("zones")
ax2.set_title("supplement-zone $\\Theta_e$ distribution", fontsize=10, color=stb.INK)
ax2.legend(loc="upper left", frameon=False, fontsize=9)
stb.style_ax(ax2)

fig.suptitle("MAD Ma+0.94_4000 rh80 — one-line $P_B$ fix, field level "
             "(fixed $M_\\mathrm{unit}$)", fontsize=12, color=stb.INK)
fig.savefig(OUT_PNG, dpi=150, facecolor=stb.SURFACE)
print(f"[fig] {OUT_PNG}")
