#!/usr/bin/env python3
import re, subprocess, h5py, numpy as np, sys, pathlib

IP = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
# Tune to these (Fig. 9 uses TOTAL polarized flux I ≡ /Ftot)
TARGETS = {
    "pos0": 0.53,
    "pos1": 0.52,
}

# Empirical nonlinearity from your runs (pos0: 7.06e28 -> 5.608e28 gave 0.667->0.460)
P_EXP = 1.60   # F ≈ const * (M_unit)^p in this range

def read_f(par_path):
    """Read outfile from .par, then /Ftot from the H5."""
    out = None
    for line in pathlib.Path(par_path).read_text().splitlines():
        m = re.match(r"^\s*outfile\s+(\S+)", line)
        if m: out = m.group(1)
    if out is None:
        raise RuntimeError("No 'outfile' in " + par_path)
    with h5py.File(out, "r") as H:
        Ftot = float(H["/Ftot"][()])
    return out, Ftot

def read_mu(par_path):
    txt = pathlib.Path(par_path).read_text()
    m = re.search(r"(?m)^\s*M_unit\s+([Ee0-9.+-]+)", txt)
    if not m:
        raise RuntimeError("No 'M_unit' in " + par_path)
    return float(m.group(1))

def write_mu(par_path, mu_new):
    p = pathlib.Path(par_path)
    txt = p.read_text()
    txt = re.sub(r"(?m)^\s*M_unit\s+[Ee0-9.+-]+", f"M_unit {mu_new:.6E}", txt)
    p.write_text(txt)

def run_ipole(par_path):
    print(f"\n=== {IP} -par {par_path}")
    # RHEL python3.6: use universal_newlines instead of text=
    p = subprocess.run([IP, "-par", par_path],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    print(p.stdout)
    return p.returncode == 0

def paper_metrics(h5_path):
    with h5py.File(h5_path, "r") as H:
        scale = float(H['/header/scale'][()])
        pol   = np.array(H['/pol']).transpose(1,0,2) * scale  # Jy/pixel
        Ftot  = float(H['/Ftot'][()])
        I = pol[:,:,0].sum(); Q = pol[:,:,1].sum()
        U = pol[:,:,2].sum(); V = pol[:,:,3].sum()
        PoverI = (Q*Q + U*U)**0.5 / Ftot
        VoverI = V / Ftot
        return Ftot, PoverI, VoverI

def tune_once(tag, par_path):
    target = TARGETS[tag]
    # 1) ensure there is a current H5 and read its Ftot
    out_h5, F_now = read_f(par_path)
    mu_now = read_mu(par_path)
    print(f"{tag}: current M_unit={mu_now:.6E}  Ftot={F_now:.6f}  target={target:.3f}")

    # 2) non-linear correction
    factor = (target / F_now) ** (1.0 / P_EXP)
    mu_new = mu_now * factor
    print(f"{tag}: propose M_unit -> {mu_new:.6E}  (factor {factor:.3f})")

    # 3) rewrite .par and run ipole once
    write_mu(par_path, mu_new)
    ok = run_ipole(par_path)
    if not ok: sys.exit(2)

    # 4) report paper-style metrics
    F, PoverI, VoverI = paper_metrics(out_h5)
    print(f"{tag}: AFTER  Ftot={F:.6f}  P/I={PoverI:.3e}  V/I={VoverI:.3e}")

if __name__ == "__main__":
    POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
    POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"
    tune_once("pos0", POS0_PAR)
    tune_once("pos1", POS1_PAR)
