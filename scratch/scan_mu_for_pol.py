#!/usr/bin/env python3
# /work/vmo703/scratch/scan_mu_for_pol.py
import os, re, tempfile, shutil, subprocess, pathlib

IP = "/work/vmo703/aricarte/run_ipole.sh"
POS0_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"

TARGETS = {"pos0": 0.53, "pos1": 0.52}
FACTORS = {
    "pos0": [0.96, 1.00, 1.04],
    "pos1": [1.05, 1.10, 1.20, 1.30],
}

FLOAT = r"[-+]?(\d+(\.\d*)?|\.\d+)([Ee][-+]?\d+)?"

def read_base_mu(par_path):
    text = pathlib.Path(par_path).read_text()
    m = re.search(r"(?m)^\s*M_unit\s+({})\s*$".format(FLOAT), text)
    if not m: raise RuntimeError("M_unit not found in {}".format(par_path))
    return float(m.group(1))

def patch_par_to_tmp(par_path, new_mu):
    # create a temp copy of the .par with updated M_unit and quench_output 1
    text = pathlib.Path(par_path).read_text()
    text = re.sub(r"(?m)^\s*M_unit\s+({})\s*$".format(FLOAT),
                  "M_unit {:.6E}".format(new_mu), text)
    if re.search(r"(?m)^\s*quench_output\s+({})\s*$".format(FLOAT), text):
        text = re.sub(r"(?m)^\s*quench_output\s+({})\s*$".format(FLOAT),
                      "quench_output 1", text)
    else:
        text += "\nquench_output 1\n"
    # (optional) ensure outfile exists but won’t be used
    # leave as-is since quench_output=1 skips writing .h5

    fd, tmp_path = tempfile.mkstemp(prefix="scanmu_", suffix=".par")
    os.close(fd)
    pathlib.Path(tmp_path).write_text(text)
    return tmp_path

def run_ipole_and_parse(par_path):
    p = subprocess.run([IP, "-par", par_path],
                       stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT,
                       universal_newlines=True)
    out = p.stdout
    # Always show helpful lines
    for line in out.splitlines():
        if line.strip().startswith("I,Q,U,V [Jy]:") or line.strip().startswith("LP,CP [%]:"):
            print("     " + line.strip())
    if p.returncode != 0:
        print("     [ipole exit code {}] Showing last 60 lines:".format(p.returncode))
        tail = "\n".join(out.splitlines()[-60:])
        print(tail)
        return None
    # Parse the last "I,Q,U,V [Jy]:" line
    iq_uv = None
    for line in out.splitlines():
        if line.strip().startswith("I,Q,U,V [Jy]:"):
            iq_uv = line
    if iq_uv is None:
        print("     Could not find 'I,Q,U,V [Jy]:' in output")
        return None
    nums = re.findall(FLOAT, iq_uv)
    # re.findall with groups returns tuples; reconstruct floats from the matched substrings
    floats = []
    for m in re.finditer(FLOAT, iq_uv):
        floats.append(float(m.group(0)))
    if len(floats) < 4:
        print("     Parse failed for I,Q,U,V")
        return None
    I, Q, U, V = floats[:4]
    P_over_I = (Q*Q + U*U) ** 0.5 / I if I != 0 else float('nan')
    V_over_I = V / I if I != 0 else float('nan')
    return I, P_over_I, V_over_I

def scan_case(tag, par_path, factors):
    target = TARGETS[tag]
    base_mu = read_base_mu(par_path)
    print("\n=== {} (target Ftot = {:.3f}) ===".format(tag, target))
    print("base M_unit = {:.6E}".format(base_mu))
    print("{:<8} {:>12} {:>10} {:>12} {:>12}".format("factor","M_unit","Ftot","P/I","V/I"))
    for f in factors:
        mu = base_mu * f
        tmp_par = patch_par_to_tmp(par_path, mu)
        print("  -> running", pathlib.Path(tmp_par).name)
        res = run_ipole_and_parse(tmp_par)
        try:
            os.remove(tmp_par)
        except Exception:
            pass
        if res is None:
            print("{:<8} {:>12} {:>10} {:>12} {:>12}".format("{:.3f}".format(f),
                  "{:.3E}".format(mu), "ERROR", "", ""))
            continue
        I, P_I, V_I = res
        print("{:<8} {:>12} {:>10.6f} {:>12.3e} {:>12.3e}".format(
            "{:.3f}".format(f), "{:.3E}".format(mu), I, P_I, V_I))

if __name__ == "__main__":
    scan_case("pos0", POS0_PAR, FACTORS["pos0"])
    scan_case("pos1", POS1_PAR, FACTORS["pos1"])
