#!/usr/bin/env python3
# Minimal, “no magic” Fig. 9 reproducer/debugger for Sa-0.5_5000 RBETA pos0/1

import csv, os, re, subprocess, h5py, numpy as np
from pathlib import Path

# --- edit these if your paths differ ---
CSV      = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/munits_table.csv"
IP       = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
DUMP     = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"
OUTDIR   = "/work/vmo703/ipole_outputs"
# ---------------------------------------

TARGETS  = {0: 0.53, 1: 0.52}   # Fig. 9 I (unpolarized) targets in Jy
FREQ     = "2.28e11"            # 228 GHz used in the figure
THETA    = "163"
NXNY     = "320"
FOV      = "160"
MBH      = "6.2e9"
DSOURCE  = "1.69e7"             # 16.9 Mpc
RCAM     = "1000"

Path(OUTDIR).mkdir(parents=True, exist_ok=True)

def read_munit_used(csv_path, pos_ratio):
    """
    Return M_unitUsed for: timestep=5000, SANE, RBETA, spin=-0.5, positron=pos_ratio (0 or 1)
    Data start at line 4 in your file; we ignore earlier lines.
    """
    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))
    rows = rows[3:]  # skip the 3 header lines exactly as you noted
    for r in rows:
        # columns (per your paste): [Timestep,MAD/SANE,Model,Spin,positron frac,MunitOffset,MunitSlope,Munit ,MunitUsed,...]
        try:
            ts = int(r[0].strip())
            fam = r[1].strip()
            model = r[2].strip()
            spin = r[3].strip()
            pos  = int(r[4].strip())
        except Exception:
            continue
        if ts==5000 and fam=="SANE" and model=="RBETA" and spin=="-0.5" and pos==pos_ratio:
            return r, float(r[8])  # 8th col is MunitUsed
    raise RuntimeError("Could not find the CSV row for pos={}.".format(pos_ratio))

def par_text(out_h5, m_unit, pos_ratio):
    # absolutely minimal .par for RBETA runs at the Fig. 9 settings
    return f"""# Fig 9: RBETA Sa-0.5_5000, positronRatio={pos_ratio}
dump {DUMP}
outfile {out_h5}

# image/camera
nx {NXNY}
ny {NXNY}
freqcgs {FREQ}
thetacam {THETA}
phicam 0
fovx_dsource {FOV}
fovy_dsource {FOV}
rcam {RCAM}
counterjet 0

# M87 geometry
MBH {MBH}
dsource {DSOURCE}

# emission & e- temperature (RBETA)
emission_type 4
electronModel RBETA
trat_small 1
trat_large 20
beta_crit 1
beta_crit_coefficient 0.5
sigma_cut 2

# pairs & scaling
positronRatio {pos_ratio}
M_unit {m_unit:.6E}
"""

def run_ipole(par_file):
    print(f"\n=== running: {IP} -par {par_file}")
    # Don’t rely on stdout parsing; we’ll read the HDF5 after
    subprocess.check_call([IP, "-par", par_file])

def read_flux_and_fractions(h5_path):
    with h5py.File(h5_path, "r") as H:
        # these two should match: sum(unpol)*scale = /Ftot_unpol
        f_unpol = float(H["/Ftot_unpol"][()])
        scale   = float(H["/header/scale"][()])
        unpol   = np.array(H["/unpol"][()]).T * scale
        pol     = np.array(H["/pol"][()]).T * scale  # shape (ny,nx,4)
        I,Q,U,V = pol[:,:,0], pol[:,:,1], pol[:,:,2], pol[:,:,3]
        I_sum = float(I.sum())    # already in Jy (scale applied)
        Q_sum = float(Q.sum())
        U_sum = float(U.sum())
        V_sum = float(V.sum())
        P_over_I = np.hypot(Q_sum, U_sum)/I_sum if I_sum!=0 else np.nan
        V_over_I = V_sum/I_sum if I_sum!=0 else np.nan
        return f_unpol, I_sum, P_over_I, V_over_I

def one_case(pos_ratio):
    csv_row, mu = read_munit_used(CSV, pos_ratio)
    out_h5  = f"{OUTDIR}/fig9_rbeta_Sa-0.5_5000_pos{pos_ratio}.h5"
    par     = f"{OUTDIR}/fig9_rbeta_Sa-0.5_5000_pos{pos_ratio}.par"

    print("\n--------------------------")
    print(f"POS{pos_ratio}  (row from CSV): {csv_row}")
    print(f"Using M_unit (CSV MunitUsed) = {mu:.3E}")

    Path(par).write_text(par_text(out_h5, mu, pos_ratio))
    run_ipole(par)
    F0, I0, P0, V0 = read_flux_and_fractions(out_h5)
    print(f"Produced:  Ftot_unpol={F0:.6f} Jy  I={I0:.6f} Jy  P/I={P0:.3e}  V/I={V0:.3e}")

    # Optional: single “brute-force” rescale to match Fig. 9 target flux
    target = TARGETS[pos_ratio]
    mu_new = mu * (target / F0)
    if abs(F0 - target)/target > 0.02:  # only rescale if off by >2%
        print(f"Rescaling M_unit once to hit target {target:.2f} Jy:")
        print(f"  M_unit: {mu:.3E}  ->  {mu_new:.3E}")
        Path(par).write_text(par_text(out_h5, mu_new, pos_ratio))
        run_ipole(par)
        F1, I1, P1, V1 = read_flux_and_fractions(out_h5)
        print(f"After rescale: Ftot_unpol={F1:.6f} Jy  I={I1:.6f} Jy  P/I={P1:.3e}  V/I={V1:.3e}")
    else:
        print("Within 2% of target; no rescale.")

if __name__ == "__main__":
    one_case(0)
    one_case(1)
