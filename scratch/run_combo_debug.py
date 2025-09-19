#!/usr/bin/env python3
# run_combo_debug.py
import pathlib as p, re, subprocess, tempfile, h5py, numpy as np

IP = "/work/vmo703/aricarte/run_ipole.sh"
PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"   # change to pos1 to test that too
TARGET_F = 0.530  # 0.520 for pos1

def read_param(txt, key, default=None):
    m = re.search(r"(?m)^\s*"+re.escape(key)+r"\s+(\S+)", txt)
    return m.group(1) if m else default

def set_or_replace(txt, key, val):
    line = f"{key} {val}"
    if re.search(r"(?m)^\s*"+re.escape(key)+r"\s+\S+", txt):
        return re.sub(r"(?m)^\s*"+re.escape(key)+r"\s+\S+.*$", line, txt)
    else:
        return txt.rstrip()+"\n"+line+"\n"

def run(par_path):
    p1 = subprocess.Popen([IP, "-par", str(par_path)],
                          stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                          universal_newlines=True)
    out,_ = p1.communicate()
    return p1.returncode, out

def get_Ftot(h5):
    with h5py.File(h5,"r") as H:
        return float(H["/Ftot"][()])

def paper_metrics(h5):
    with h5py.File(h5,"r") as H:
        scale = float(H["/header/scale"][()])
        F = float(H["/Ftot"][()])
        pol = np.array(H["/pol"]).transpose(1,0,2)*scale
        I,Q,U,V = [pol[:,:,i].sum() for i in range(4)]
        return F, (Q*Q+U*U)**0.5/F, V/F

base = p.Path(PAR).read_text()

# *** pick one combo to debug ***
Rhigh = 20
bcrit = 0.5
sigma = 1.0
pos   = "0"  # "1" to test pos1

par = base
par = set_or_replace(par, "electronModel", "RBETA")
par = set_or_replace(par, "trat_small", "1")
par = set_or_replace(par, "trat_large", str(Rhigh))
par = set_or_replace(par, "sigma_cut",  str(sigma))
par = set_or_replace(par, "beta_crit",  str(bcrit))   # comment this line if it’s the culprit
par = set_or_replace(par, "positronRatio", pos)
par = set_or_replace(par, "quench_output", "0")

# write debug par so you can see exactly what we gave ipole
dbg_par = p.Path("/work/vmo703/ipole_outputs/_debug_last.par")
dbg_par.write_text(par)

# 1st run: measure F, retune M_unit
tmp1 = p.Path(tempfile.mkstemp(suffix=".par")[1])
tmp1.write_text(par)
rc, log = run(tmp1)
print("=== first run rc =", rc)
if rc != 0:
    print("--- ipole log (first run) ---")
    print(log)
    raise SystemExit(1)

out_h5 = read_param(par, "outfile")
F1 = get_Ftot(out_h5)
mu_old = float(read_param(par, "M_unit"))
mu_new = mu_old*(TARGET_F/(F1 if F1!=0 else TARGET_F))
par2 = set_or_replace(par, "M_unit", f"{mu_new:.6E}")

# 2nd run: locked flux + metrics
tmp2 = p.Path(tempfile.mkstemp(suffix=".par")[1])
tmp2.write_text(par2)
rc2, log2 = run(tmp2)
print("=== second run rc =", rc2)
if rc2 != 0:
    print("--- ipole log (second run) ---")
    print(log2)
    raise SystemExit(2)

F, P, V = paper_metrics(out_h5)
print(f"Rhigh={Rhigh}  beta_crit={bcrit}  sigma_cut={sigma}  ->  F={F:.3f}  P/I={P:.3e}  V/I={V:.3e}")
print("debug par saved to:", str(dbg_par))
