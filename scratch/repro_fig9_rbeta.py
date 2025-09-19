#!/usr/bin/env python3
import os, h5py

# make sure both the repo root and Positron scripts are importable
POS_ROOT = "/work/vmo703/aricarte-copy/aricarte"
import sys
sys.path.insert(0, f"{POS_ROOT}/PositronIPOLEScripts")
sys.path.insert(0, POS_ROOT)

from ipole_many_models import runIPOLE

# --- inputs matching Fig 9 (SANE a=-0.5, T=25000M -> dump 5000) ---
sim_file   = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"
ipole_bin  = f"{POS_ROOT}/ipole+e-/ipole"
out_base   = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000"
npixel     = 320
thetacam   = 163.0
fov        = 160.0
rmax_geo   = 50
Rhigh      = 20
freq_Hz    = 228e9
counterjet = 0

# RBETA scaling (from your CSV for SANE -0.5, 5000):
Munit       = 1.83e27
MunitOffset = 2.50e28
MunitSlope  = 2.50e01

def munit_used(pos):
    return MunitOffset + MunitSlope * Munit / (1.0 + 2.0*pos)

# --- run pos0 ---
out0 = f"{out_base}_pos0.h5"
runIPOLE(
    sim_file,                   # positional: input dump
    out0,                       # positional: output .h5
    munit_used(0),              # positional: M_unit
    ipoleExecutable=ipole_bin,
    thetacam=thetacam, Rhigh=Rhigh, freq_Hz=freq_Hz,
    fov=fov, npixel=npixel, counterjet=counterjet, rmax_geo=rmax_geo,
    # RBETA:
    electronModel=2,            # 2 is RBETA in this branch (CRITBETA uses beta_crit params)
    positronRatio=0,
    sigma_cut=2.0,              # matches unlocked script defaults
)

# --- run pos1 ---
out1 = f"{out_base}_pos1.h5"
runIPOLE(
    sim_file, out1, munit_used(1),
    ipoleExecutable=ipole_bin,
    thetacam=thetacam, Rhigh=Rhigh, freq_Hz=freq_Hz,
    fov=fov, npixel=npixel, counterjet=counterjet, rmax_geo=rmax_geo,
    electronModel=2,
    positronRatio=1,
    sigma_cut=2.0,
)

# quick sanity print
for p in (out0, out1):
    with h5py.File(p, "r") as f:
        print(os.path.basename(p), "Ftot_unpol[Jy] =", f["/Ftot_unpol"][()])
