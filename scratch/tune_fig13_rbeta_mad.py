#!/usr/bin/env python3
# tune_fig13_rbeta_mad.py
# find (A,S) for fig. 13 (MAD, a=-0.5, RBETA, t=5000) by
# 1) secant on MunitUsed per positron ratio to hit target flux,
#    with exactly ONE secant step for r=0.000 and TWO for r=1.000,
# 2) solve A,S from the two tuned MunitUsed values,
# 3) verify by re-running with the wrapper’s A,S formula.


import os, math, h5py, numpy as np
from ipole_many_models import runIPOLE


# --- config (match wrapper defaults) ---
simFile   = "/work/vmo703/grmhd_dump_samples/Ma-0.5_5000.h5"
nameBase  = "/work/vmo703/ipole_outputs/fig13_MADa-0.5_RBETA.h5"
ipole     = "/work/vmo703/aricarte/run_ipole.sh"


# RBETA settings
electronModel = 2 # RBETA = 2; CRITBETA = 4
Rhigh         = 20
sigma_cut     = 2.0
sigma_transition = 2.0 # no-jet configuration
thetacam      = 163.0
rmax_geo      = 50
freq_Hz       = 228e9
fov           = 160.0
npixel        = 320
counterjet    = 0


# base M and desired flux targets
M_base  = 7.49e24
T0, T1  = 0.55, 0.47 # targets for r=0 and r=1


def read_flux(h5path):
    with h5py.File(h5path, "r") as H:
        F  = float(H["/Ftot"][()])
        Q  = float(H["/Qt"][()]) if "/Qt" in H else None
        U  = float(H["/Ut"][()]) if "/Ut" in H else None
        V  = float(H["/Vt"][()]) if "/Vt" in H else None
    return F, Q, U, V


def run_once(positron_ratio, Munit_used):
    out = nameBase.replace(".h5", f"_{positron_ratio:1.3f}.h5")
    # run ipole via wrapper with explicit Munit_used
    runIPOLE(simFile, out, Munit_used, ipoleExecutable=ipole,
             thetacam=thetacam, Rhigh=Rhigh, freq_Hz=freq_Hz,
             fov=fov, npixel=npixel, counterjet=counterjet,
             rmax_geo=rmax_geo, positronRatio=positron_ratio,
             electronModel=electronModel, sigma_cut=sigma_cut,
             sigma_transition=sigma_transition)
    return read_flux(out)[0], out


def secant_for_flux(ratio, target, mu0, mu1, n_steps):
    """
    run exactly n_steps secant iterations starting from (mu0, mu1).
    if mu1 is None, seed mu1 by simple proportional scaling from mu0.
    """
    F0, _ = run_once(ratio, mu0)
    if mu1 is None:
        # proportional seed if a second guess is not provided
        mu1 = mu0 * (target / max(F0, 1e-30))
    F1, _ = run_once(ratio, mu1)
    print(f"[r={ratio}] start: mu0={mu0:.3e} F0={F0:.6f}  | mu1={mu1:.3e} F1={F1:.6f}")

    Mprev, Fprev = mu0, F0
    Mcurr, Fcurr = mu1, F1

    for i in range(1, n_steps+1):
        denom = (Fcurr - Fprev)
        if abs(denom) < 1e-12:
            Mnext = Mcurr * (target / max(Fcurr, 1e-30))
        else:
            Mnext = Mcurr + (target - Fcurr) * (Mcurr - Mprev) / denom
        Fnext, _ = run_once(ratio, Mnext)
        print(f"[r={ratio}] step {i}: mu={Mnext:.3e}  F={Fnext:.6f}")
        Mprev, Fprev, Mcurr, Fcurr = Mcurr, Fcurr, Mnext, Fnext

    return Mcurr, Fcurr


def verify_with_AS(A, S):
    # use the wrapper’s formula U(r)=A + S*M_base/(1+2r) to produce both r=0 and r=1
    results = {}
    for r, T in [(0.0, T0), (1.0, T1)]:
        U_used = A + S * M_base / (1.0 + 2.0*r)
        F, out = run_once(r, U_used)
        _, Q, Uq, V = read_flux(out)
        PoverI = (math.sqrt(Q**2 + Uq**2)/F) if (Q is not None and Uq is not None and F>0) else None
        VoverI = (V/F) if (V is not None and F>0) else None
        results[r] = (F, PoverI, VoverI, out)
        extra = (f"  P/I={PoverI:.3e}  V/I={VoverI:.3e}" if PoverI is not None else "")
        print(f"[verify r={int(r)}] Munit_used={U_used:.3e}  F={F:.6f} (target {T:.2f}){extra}")
    return results


if __name__ == "__main__":
    # initial guesses straight from CSV (A0,S0), only to seed the secants:
    A0 = 1.30e25
    S0 = 5.00

    # corresponding U(r) from those (seed points)
    U0_guess = A0 + S0*M_base           # r=0
    U1_guess = A0 + (S0*M_base)/3.0     # r=1

    # tune r=0 → 0.55 Jy (EXACTLY ONE secant step)
    U0_star, F0_star = secant_for_flux(0.0, T0, U0_guess, None, n_steps=1)

    # tune r=1 → 0.47 Jy (EXACTLY TWO secant steps)
    U1_star, F1_star = secant_for_flux(1.0, T1, U1_guess, None, n_steps=2)

    print(f"\nTuned Munit_used: U0*={U0_star:.6e} (F~{F0_star:.6f}),  U1*={U1_star:.6e} (F~{F1_star:.6f})")

    # solve for A,S that reproduce both simultaneously via RBETA mapping:
    # U(r) = A + S*M_base/(1+2r)
    # from r=0 and r=1 tuned values:
    #   S = (3/(2*M)) * (U0* - U1*)
    #   A = U1* - (S*M)/3
    S = 3.0*(U0_star - U1_star)/(2.0*M_base)
    A = U1_star - S*M_base/3.0

    print(f"\nDerived parameters for Fig.13 (RBETA MAD a=-0.5, t=5000):")
    print(f"  MunitOffset (A) = {A:.6e}")
    print(f"  MunitSlope  (S) = {S:.6e}")

    # verify by running again using U(r)=A + S*M/(1+2r)
    print("\nVerifying with these (A,S)…")
    verify_with_AS(A, S)
