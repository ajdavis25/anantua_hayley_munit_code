#!/usr/bin/env python3
# /work/vmo703/scratch/scan_te_beta_fixed_flux.py
import os, re, tempfile, pathlib, subprocess, math, argparse

IP = "/work/vmo703/aricarte/run_ipole.sh"
POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"

TARGETS = {
    "pos0": {"F":0.53, "PoverI":3.93e-3, "VoverI":-5.60e-3, "par":POS0_PAR},
    "pos1": {"F":0.52, "PoverI":1.62e-3, "VoverI":-2.81e-3, "par":POS1_PAR},
}

# modest grids; expand if needed
TRAT_LARGE_GRID = [15, 20, 25]     # Rhigh
BETA_CRIT_GRID  = [0.5, 1.0]       # beta_crit threshold
SIGMA_CUT       = 2.0              # hold fixed for now

FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([Ee][-+]?\d+)?"

def read_mu(par_path):
    txt = pathlib.Path(par_path).read_text()
    m = re.search(r"(?m)^\s*M_unit\s+({})\s*$".format(FLOAT), txt)
    if not m: raise RuntimeError("M_unit not found in "+par_path)
    return float(m.group(1))

def patch(par_path, **kv):
    txt = pathlib.Path(par_path).read_text()
    # ensure keys exist / are replaced
    def set_key(k,v):
        nonlocal txt
        pat = r"(?m)^\s*{}\s+{}\s*$".format(re.escape(k), FLOAT)
        line = "{} {}".format(k, v)
        if re.search(pat, txt):
            txt = re.sub(pat, line, txt)
        else:
            txt += "\n" + line + "\n"

    for k,v in kv.items(): set_key(k,str(v))

    # enforce our fixed knobs
    set_key("sigma_cut", SIGMA_CUT)

    # quench output (numbers only, no .h5)
    if re.search(r"(?m)^\s*quench_output\s+{}".format(FLOAT), txt):
        txt = re.sub(r"(?m)^\s*quench_output\s+{}".format(FLOAT), "quench_output 1", txt)
    else:
        txt += "\nquench_output 1\n"

    fd,tmp = tempfile.mkstemp(prefix="scanfix_", suffix=".par")
    os.close(fd)
    pathlib.Path(tmp).write_text(txt)
    return tmp

def run_ipole(par):
    p = subprocess.run([IP, "-par", par],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    out = p.stdout
    I=Q=U=V=None
    for line in out.splitlines():
        if line.strip().startswith("I,Q,U,V [Jy]:"):
            nums = [float(m.group(0)) for m in re.finditer(FLOAT, line)]
            if len(nums) >= 4:
                I,Q,U,V = nums[:4]
    return p.returncode, out, I, Q, U, V

def retune_mu_to_flux(base_par, target_flux, mu_start, trat_large, beta_crit,
                      tol_rel=0.005, max_iters=6):
    mu = mu_start
    last_I = None
    for _ in range(max_iters):
        tmp = patch(base_par, M_unit="{:.6E}".format(mu),
                    trat_large=str(trat_large),
                    beta_crit=str(beta_crit))
        rc,out,I,Q,U,V = run_ipole(tmp)
        try: os.remove(tmp)
        except: pass
        if rc != 0 or I is None:
            return None, ("ipole failure\n" + out)
        err = (I - target_flux)/target_flux
        if abs(err) <= tol_rel:
            P_I = math.hypot(Q,U)/I if I else float('nan')
            V_I = V/I if I else float('nan')
            return {"mu": mu, "F": I, "P_I": P_I, "V_I": V_I, "stdout": out}, None
        # multiplicative correction (assumes near-linear response)
        # damp a bit to avoid oscillation in non-linear regime
        mu *= (target_flux / I)**0.9
        last_I = I
    # final report even if not within tol
    P_I = math.hypot(Q,U)/I if I else float('nan')
    V_I = V/I if I else float('nan')
    return {"mu": mu, "F": I, "P_I": P_I, "V_I": V_I, "stdout": out, "note":"max iters"}, None

def score(P_I, V_I, target_P, target_V):
    # weighted distance in (P/I, V/I) space (equal weights)
    return math.hypot(P_I - target_P, V_I - target_V)

def run_case(tag):
    tgt = TARGETS[tag]
    par = tgt["par"]
    F_tgt = tgt["F"]
    P_tgt = tgt["PoverI"]
    V_tgt = tgt["VoverI"]
    mu0 = read_mu(par)

    print(f"\n=== {tag} @ fixed flux {F_tgt:.3f} ===")
    print(f"start M_unit = {mu0:.6E}")
    print("{:>7} {:>6} {:>12} {:>8} {:>12} {:>12} {:>12}".format(
        "Rhigh","bcrit","M_unit","Ftot","P/I","V/I","|Δ|"))
    best = None

    for R in TRAT_LARGE_GRID:
        for bc in BETA_CRIT_GRID:
            res,err = retune_mu_to_flux(par, F_tgt, mu0, R, bc)
            if err:
                print("{:>7} {:>6} {:>12} {:>8} {:>12} {:>12} {:>12}".format(
                    f"{R}", f"{bc}", "ERROR", "", "", "", ""))
                tail = "\n".join(err.splitlines()[-30:])
                print(tail)
                continue
            d = score(res["P_I"], res["V_I"], P_tgt, V_tgt)
            print("{:>7} {:>6} {:>12.3E} {:>8.3f} {:>12.3e} {:>12.3e} {:>12.3e}".format(
                f"{R}", f"{bc}", res["mu"], res["F"], res["P_I"], res["V_I"], d))
            cand = (d, R, bc, res)
            if best is None or d < best[0]:
                best = cand

    if best:
        d,R,bc,res = best
        print("\n>>> best {}: Rhigh={}  beta_crit={}  ->  F={:.3f}  P/I={:.3e}  V/I={:.3e}  (M_unit={:.3E})".format(
            tag, R, bc, res["F"], res["P_I"], res["V_I"], res["mu"]))
    else:
        print(f"\n(no valid runs for {tag})")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=["pos0","pos1","both"], default="both")
    args = ap.parse_args()
    if args.case in ("pos0","both"):
        run_case("pos0")
    if args.case in ("pos1","both"):
        run_case("pos1")
