#!/usr/bin/env python3
# Lock Sa-0.5_5000 RBETA (pos0/pos1) to the Fig. 9 fluxes using your CSV M_unitUsed
# 1) run with CSV M_unitUsed
# 2) read Ftot_unpol
# 3) set M_unit := M_unit * (target / Ftot_unpol) and rerun
# 4) print I, P/I, V/I

import csv, subprocess, h5py, numpy as np
from pathlib import Path

# --- paths (edit if yours differ) ---
CSV    = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/munits_table.csv"
IP = "/work/vmo703/aricarte/run_ipole.sh"
DUMP   = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"
OUTDIR = "/work/vmo703/ipole_outputs"
# ------------------------------------

TARGET = {0: 0.53, 1: 0.52}  # Jy, Fig. 9
CFG = dict(
    freq="2.28e11", thetacam="163", nx="320", fov="160",
    MBH="6.2e9", dsource="1.69e7", rcam="1000",
)

Path(OUTDIR).mkdir(parents=True, exist_ok=True)

def csv_munit_used(csv_path, pos):
    # Data start at line 4 in your file
    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))[3:]
    for r in rows:
        try:
            if (int(r[0])==5000 and r[1].strip()=="SANE" and r[2].strip()=="RBETA"
                and r[3].strip()=="-0.5" and int(r[4])==pos):
                return r, float(r[8])
        except Exception:
            continue
    raise RuntimeError(f"CSV row not found for pos={pos}")

def make_par(out_h5, munit, pos):
    return f"""# Fig9 RBETA Sa-0.5_5000 pos{pos}
dump {DUMP}
outfile {out_h5}

nx {CFG['nx']}
ny {CFG['nx']}
freqcgs {CFG['freq']}
thetacam {CFG['thetacam']}
phicam 0
fovx_dsource {CFG['fov']}
fovy_dsource {CFG['fov']}
rcam {CFG['rcam']}
counterjet 0

MBH {CFG['MBH']}
dsource {CFG['dsource']}

emission_type 4
electronModel RBETA
trat_small 1
trat_large 20
# RBETA uses the same ‘beta_crit’ threshold internally; this is normal:
beta_crit 1
beta_crit_coefficient 0.5
sigma_cut 2

positronRatio {pos}
M_unit {munit:.6E}
"""

def run(par):
    print(f"\n=== {IP} -par {par}")
    subprocess.check_call([IP, "-par", par])

def read_metrics(h5):
    with h5py.File(h5, "r") as H:
        f_unpol = float(H["/Ftot_unpol"][()])
        scale   = float(H["/header/scale"][()])
        pol     = np.array(H["/pol"][()]).T * scale
        I,Q,U,V = pol[:,:,0], pol[:,:,1], pol[:,:,2], pol[:,:,3]
        Is,Qs,Us,Vs = I.sum(), Q.sum(), U.sum(), V.sum()
        PoverI = np.hypot(Qs,Us)/Is if Is!=0 else np.nan
        VoverI = Vs/Is if Is!=0 else np.nan
        return f_unpol, Is, PoverI, VoverI

def do_one(pos):
    row, mu = csv_munit_used(CSV, pos)
    out = f"{OUTDIR}/fig9_rbeta_Sa-0.5_5000_pos{pos}.h5"
    par = f"{OUTDIR}/fig9_rbeta_Sa-0.5_5000_pos{pos}.par"

    print("\n-----------------------------")
    print(f"POS{pos} CSV row: {row}")
    print(f"Initial M_unit (CSV M_unitUsed): {mu:.3E}")

    # pass 1: CSV value
    Path(par).write_text(make_par(out, mu, pos))
    run(par)
    F0, I0, P0, V0 = read_metrics(out)
    print(f"CSV pass:     Ftot_unpol={F0:.6f}  I={I0:.6f}  P/I={P0:.3e}  V/I={V0:.3e}")

    # pass 2: exact rescale to Fig. 9 target
    t = TARGET[pos]
    mu2 = mu * (t / F0)
    Path(par).write_text(make_par(out, mu2, pos))
    print(f"Rescale M_unit: {mu:.3E} -> {mu2:.3E}  (aiming for {t:.2f} Jy)")
    run(par)
    F1, I1, P1, V1 = read_metrics(out)
    print(f"LOCKED flux:  Ftot_unpol={F1:.6f}  I={I1:.6f}  P/I={P1:.3e}  V/I={V1:.3e}")

if __name__ == "__main__":
    do_one(0)  # pos0
    do_one(1)  # pos1
