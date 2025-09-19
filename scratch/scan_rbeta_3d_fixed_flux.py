#!/usr/bin/env python3
# scan_rbeta_3d_fixed_flux.py
# Small 3D scan for RBETA at fixed flux:
#   (Rhigh, beta_crit, sigma_cut) and retune M_unit to target Ftot each time.
#
# Uses the same "paper metrics" definition: Ftot from /Ftot; P/I, V/I from summed Stokes (scaled).

import os, re, sys, tempfile, subprocess, shutil, h5py, numpy as np, pathlib as p
IP = "/work/vmo703/aricarte/run_ipole.sh"
POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"

# modest neighborhoods around what you found
Rhigh_grid   = [15, 20, 25]
bcrit_grid   = [0.3, 0.5, 0.8]
sigma_grid   = [0.5, 1.0, 1.5, 2.0]  # lower = more low-σ plasma -> generally more |V/I| and lower P/I

targets = {
  "pos0": 0.530,
  "pos1": 0.520,
}

def run_ipole(par_path):
    # Python 3.6 compatible capture
    p1 = subprocess.Popen([IP, "-par", par_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    out,_ = p1.communicate()
    return p1.returncode, out

def get_Ftot(h5):
    with h5py.File(h5, "r") as H:
        return float(H["/Ftot"][()])

def paper_metrics(h5):
    with h5py.File(h5, "r") as H:
        scale = float(H["/header/scale"][()])
        Ftot  = float(H["/Ftot"][()])
        pol   = np.array(H["/pol"]).transpose(1,0,2) * scale  # Jy/pixel
        I = pol[:,:,0].sum()
        Q = pol[:,:,1].sum()
        U = pol[:,:,2].sum()
        V = pol[:,:,3].sum()
        # sanity
        if Ftot != 0.0:
            rel = abs(I - Ftot)/abs(Ftot)
            if rel > 5e-6:
                print("  [warn] sum(I) and /Ftot differ by {:.2e}".format(rel))
        PoverI = (Q*Q + U*U)**0.5 / (Ftot if Ftot != 0.0 else 1.0)
        VoverI = V / (Ftot if Ftot != 0.0 else 1.0)
        return Ftot, PoverI, VoverI

def read_param(par_text, key, default=None):
    m = re.search(r"(?m)^\s*"+re.escape(key)+r"\s+(\S+)", par_text)
    return m.group(1) if m else default

def set_or_replace(par_text, key, value):
    line = f"{key} {value}"
    if re.search(r"(?m)^\s*"+re.escape(key)+r"\s+\S+", par_text):
        return re.sub(r"(?m)^\s*"+re.escape(key)+r"\s+\S+.*$", line, par_text)
    else:
        return par_text.rstrip()+"\n"+line+"\n"

def retune_once(par_text, target_flux):
    # 1) write temp par, run, measure F, compute new M_unit, 2) write final par with tuned M_unit
    tmp_par = p.Path(tempfile.mkstemp(suffix=".par")[1])
    out_h5  = read_param(par_text, "outfile")
    if not out_h5:
        raise RuntimeError("No 'outfile' in par text")
    tmp_txt = par_text
    # ensure quench_output 0 so we get HDF5
    tmp_txt = set_or_replace(tmp_txt, "quench_output", "0")
    tmp_par.write_text(tmp_txt)
    rc, log = run_ipole(str(tmp_par))
    if rc != 0 or not p.Path(out_h5).exists():
        return None, None, None, log

    F = get_Ftot(out_h5)
    mu_old = float(read_param(par_text, "M_unit"))
    mu_new = mu_old * (target_flux / (F if F != 0.0 else target_flux))
    tuned = set_or_replace(par_text, "M_unit", f"{mu_new:.6E}")

    # second run with tuned M_unit
    tmp_par2 = p.Path(tempfile.mkstemp(suffix=".par")[1])
    tmp_par2.write_text(tuned)
    rc2, log2 = run_ipole(str(tmp_par2))
    if rc2 != 0 or not p.Path(out_h5).exists():
        return None, None, None, log + "\n" + (log2 or "")
    F2, PoverI, VoverI = paper_metrics(out_h5)
    # cleanup temp pars (keep .h5 outputs for inspection)
    try:
        tmp_par.unlink()
        tmp_par2.unlink()
    except Exception:
        pass
    return F2, PoverI, VoverI, (log + "\n" + log2)

def scan_one(base_par_path, tag, R_set, B_set, S_set):
    base = p.Path(base_par_path).read_text()
    tgt  = targets[tag]
    # make sure required fields exist
    assert read_param(base, "electronModel") == "RBETA", "electronModel must be RBETA for this scan"
    print("\n=== {} @ fixed flux {:.3f} ===".format(tag, tgt))
    mu0 = float(read_param(base, "M_unit"))
    print("start M_unit =", "{:.6E}".format(mu0))
    print(" Rhigh  bcrit  sigma  |    M_unit    |   Ftot   |    P/I    |    V/I    |  |Δ| to target")
    best = None

    # target P/I & V/I from Fig. 9
    P_tgt, V_tgt = ((3.93e-3, -5.60e-3) if tag=="pos0" else (1.62e-3, -2.81e-3))

    for R in R_set:
        for B in B_set:
            for S in S_set:
                # modify base par
                par = base
                par = set_or_replace(par, "trat_small", "1")
                par = set_or_replace(par, "trat_large", str(R))
                par = set_or_replace(par, "beta_crit", str(B))
                par = set_or_replace(par, "sigma_cut", str(S))
                # IMPORTANT: ensure positronRatio is correct for tag (in case the base par has wrong value)
                par = set_or_replace(par, "positronRatio", "0" if tag=="pos0" else "1")

                F, PoverI, VoverI, _log = retune_once(par, tgt)
                if F is None:
                    print(f"  R={R:>3}  b={B:>3}  s={S:>3} : ERROR")
                    continue
                mu_new = read_param(par, "M_unit")
                d = abs(PoverI - P_tgt) + abs(VoverI - V_tgt)
                print(f"  {R:>5}  {B:>4}  {S:>5}  | {float(mu_new):.3E} |  {F:6.3f}  | {PoverI:8.3e} | {VoverI:8.3e} | {d:8.3e}")
                rec = (d, R, B, S, float(mu_new), F, PoverI, VoverI)
                if best is None or d < best[0]:
                    best = rec

    if best:
        _, R, B, S, MU, F, P, V = best
        print(f"\n>>> best {tag}: Rhigh={R}  beta_crit={B}  sigma_cut={S}  ->  F={F:.3f}  P/I={P:.3e}  V/I={V:.3e}  (M_unit={MU:.3E})")
    else:
        print(f"\n(no successful runs for {tag})")

def main():
    # shrink grids if you want quicker passes
    scan_one(POS0_PAR, "pos0", Rhigh_grid, bcrit_grid, sigma_grid)
    scan_one(POS1_PAR, "pos1", Rhigh_grid, bcrit_grid, sigma_grid)

if __name__ == "__main__":
    main()
