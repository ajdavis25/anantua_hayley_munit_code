#!/usr/bin/env python3
# /work/vmo703/scratch/scan_phys_knob_at_fixed_flux.py
import os, re, tempfile, pathlib, subprocess

IP = "/work/vmo703/aricarte/run_ipole.sh"
POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"

TARGET = {"pos0": 0.53, "pos1": 0.52}
SIGMA_CUTS = [1.0, 1.5, 2.0, 3.0]   # expand if needed

FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([Ee][-+]?\d+)?"

def read_mu(par_path):
    txt = pathlib.Path(par_path).read_text()
    m = re.search(r"(?m)^\s*M_unit\s+({})\s*$".format(FLOAT), txt)
    if not m: raise RuntimeError("M_unit not found in "+par_path)
    return float(m.group(1))

def patch(par_path, **kv):
    txt = pathlib.Path(par_path).read_text()
    # set/replace simple scalar keys in the par
    for key,val in kv.items():
        pat = r"(?m)^\s*{}\s+{}\s*$".format(re.escape(key), FLOAT)
        line = "{} {}".format(key, val)
        if re.search(pat, txt):
            txt = re.sub(pat, line, txt)
        else:
            txt += "\n" + line + "\n"
    # force: quench_output=1 (no .h5 writes)
    if re.search(r"(?m)^\s*quench_output\s+{}".format(FLOAT), txt):
        txt = re.sub(r"(?m)^\s*quench_output\s+{}".format(FLOAT), "quench_output 1", txt)
    else:
        txt += "\nquench_output 1\n"
    fd,tmp = tempfile.mkstemp(prefix="scanphys_", suffix=".par")
    os.close(fd)
    pathlib.Path(tmp).write_text(txt)
    return tmp

def run_ipole(par):
    p = subprocess.run([IP, "-par", par],
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                       universal_newlines=True)
    out = p.stdout
    # parse last "I,Q,U,V [Jy]:" line
    I=Q=U=V=None
    for line in out.splitlines():
        if line.strip().startswith("I,Q,U,V [Jy]:"):
            nums = []
            for m in re.finditer(FLOAT, line):
                nums.append(float(m.group(0)))
            if len(nums) >= 4:
                I,Q,U,V = nums[:4]
    return p.returncode, out, I, Q, U, V

def retune_mu_to_flux(base_par, target_flux, mu_start, sigma_cut):
    # do up to 2 multiplicative corrections to hit target flux quickly
    mu = mu_start
    for _ in range(2):
        tmp = patch(base_par, M_unit="{:.6E}".format(mu), sigma_cut=str(sigma_cut))
        rc,out,I,Q,U,V = run_ipole(tmp)
        try: os.remove(tmp)
        except: pass
        if rc != 0 or I is None:
            return None, ("ipole failed while tuning\n" + out)
        if I == 0: break
        mu *= (target_flux / I)
    # one final run at the tuned mu to report numbers
    tmp = patch(base_par, M_unit="{:.6E}".format(mu), sigma_cut=str(sigma_cut))
    rc,out,I,Q,U,V = run_ipole(tmp)
    try: os.remove(tmp)
    except: pass
    if rc != 0 or I is None:
        return None, ("ipole failed on final run\n" + out)
    P_I = (Q*Q + U*U)**0.5 / I if I else float('nan')
    V_I = V / I if I else float('nan')
    return {"mu": mu, "Ftot": I, "P_I": P_I, "V_I": V_I, "stdout": out}, None

def do_case(tag, par):
    target = TARGET[tag]
    mu0 = read_mu(par)
    print("\n=== {} @ fixed flux {:.3f} ===".format(tag, target))
    print("start M_unit = {:.6E}".format(mu0))
    print("{:>8} {:>12} {:>10} {:>10} {:>10}".format("sigma", "M_unit", "Ftot", "P/I", "V/I"))
    for sc in SIGMA_CUTS:
        res,err = retune_mu_to_flux(par, target, mu0, sc)
        if err:
            print("{:>8} {:>12} {:>10} {:>10} {:>10}".format(
                "{:.2f}".format(sc), "ERROR", "", "", ""))
            # tail the error lines to help debug
            tail = "\n".join(err.splitlines()[-40:])
            print(tail)
            continue
        print("{:>8} {:>12.3E} {:>10.3f} {:>10.3e} {:>10.3e}".format(
            "{:.2f}".format(sc), res["mu"], res["Ftot"], res["P_I"], res["V_I"]))

if __name__ == "__main__":
    do_case("pos0", POS0_PAR)
    do_case("pos1", POS1_PAR)
