#!/usr/bin/env python3
# tune_munit_rbeta_fast.py
#
# Finds M_unit values that reproduce Fig. 9 (RBETA, Sa-0.5 @ dump 5000)
# and ONLY keeps the two final 320×320 outputs. All intermediate .h5 are deleted.

import os, re, tempfile, subprocess, pathlib

# --- paths (edit if your layout differs) ---
IPOLE    = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
TEMPLATE = "/work/vmo703/scratch/runIPOLE_fig9_.par"   # your Fig.9 par template
OUTDIR   = "/work/vmo703/ipole_outputs"
pathlib.Path(OUTDIR).mkdir(parents=True, exist_ok=True)

# -------- helpers --------
def _sub_line(txt, key, value):
    """Replace 'key <number/string>' line with new value (adds if missing)."""
    pat = re.compile(rf"(?m)^{re.escape(key)}\s+\S+")
    repl = f"{key} {value}"
    return pat.sub(repl, txt) if pat.search(txt) else (txt + "\n" + repl + "\n")

def write_par_from_template(out_par, outfile_h5, munit, pos, npix):
    """Render a .par from TEMPLATE with all knobs set to the Fig. 9 spec."""
    txt = pathlib.Path(TEMPLATE).read_text()

    # camera / image geometry
    txt = _sub_line(txt, "nx", npix)
    txt = _sub_line(txt, "ny", npix)
    txt = _sub_line(txt, "thetacam", 163)
    txt = _sub_line(txt, "phicam", 0)
    txt = _sub_line(txt, "fovx_dsource", 160)
    txt = _sub_line(txt, "fovy_dsource", 160)

    # frequency & model choices (RBETA @ 228 GHz)
    txt = _sub_line(txt, "freqcgs", "2.28e11")
    txt = _sub_line(txt, "electronModel", 4)   # RBETA
    txt = _sub_line(txt, "emission_type", 4)

    # positron & electron-temperature tuning used in the positron scripts
    txt = _sub_line(txt, "positronRatio", pos)
    txt = _sub_line(txt, "trat_small", 1)
    txt = _sub_line(txt, "trat_large", 20)
    txt = _sub_line(txt, "sigma_cut", 2)

    # scaling & output
    txt = _sub_line(txt, "M_unit", f"{munit:.6e}")
    txt = _sub_line(txt, "outfile", outfile_h5)

    pathlib.Path(out_par).write_text(txt)

def run_ipole_stream(par_path):
    """Run ipole and capture stdout (py3.6-safe)."""
    proc = subprocess.Popen([IPOLE, "-par", par_path], stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out, _ = proc.communicate()
    return proc.returncode, out

def flux_from_stdout_or_h5(stdout_text, h5_path):
    """Prefer flux parsed from ipole stdout; fall back to h5dump on /Ftot_unpol."""
    m = re.search(r"Ftot:\s*([0-9.]+)\s*Jy", stdout_text)
    if m:
        return float(m.group(1))
    # fallback via h5dump (no h5py dependency)
    try:
        p = subprocess.Popen(["h5dump", "-d", "/Ftot_unpol", h5_path], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out, _ = p.communicate()
        # find the first float in the DATA block
        m2 = re.search(r"DATA\s*{\s*\(0\):\s*([0-9.eE+-]+)", out)
        if m2:
            return float(m2.group(1))
    except Exception:
        pass
    raise RuntimeError("Could not determine flux from stdout or H5.")

def run_once(munit, pos, tag, npix, keep=False):
    """Run ipole once; optionally delete the output unless it's the final result."""
    outfile = f"{OUTDIR}/fig9_rbeta_pos{pos}_{tag}_{npix}px.h5"
    with tempfile.NamedTemporaryFile(delete=False, suffix=".par") as tf:
        par_path = tf.name
    write_par_from_template(par_path, outfile, munit, pos, npix)
    rc, out = run_ipole_stream(par_path)
    try:
        pathlib.Path(par_path).unlink(missing_ok=True)
    except TypeError:
        # py3.6 compatibility
        if pathlib.Path(par_path).exists():
            os.remove(par_path)
    if rc != 0:
        raise RuntimeError("ipole failed\n" + out)
    flux = flux_from_stdout_or_h5(out, outfile)
    # delete intermediates to avoid clutter
    if (not keep) and ("final" not in tag):
        try: os.remove(outfile)
        except OSError: pass
    return flux

def tune(target, pos, coarse_np=160, fine_np=320, iters=10):
    """Bracket + bisection at 160px, confirm at 320px. Returns (M_unit, Ftot)."""
    lo, hi = 1e27, 1e29
    f_lo = run_once(lo, pos, "lo",  coarse_np, keep=False)
    f_hi = run_once(hi, pos, "hi",  coarse_np, keep=False)
    # expand bracket if needed
    safety = 0
    while not (f_lo <= target <= f_hi):
        if target < f_lo:
            hi, f_hi = lo, f_lo
            lo /= 10.0
            f_lo = run_once(lo, pos, "loexp", coarse_np, keep=False)
        else:
            lo, f_lo = hi, f_hi
            hi *= 10.0
            f_hi = run_once(hi, pos, "hiexp", coarse_np, keep=False)
        safety += 1
        if safety > 12:
            raise RuntimeError("Failed to bracket target flux.")
    # bisection
    for k in range(iters):
        mid = 0.5*(lo + hi)
        f_mid = run_once(mid, pos, f"mid{k+1}", coarse_np, keep=False)
        if f_mid < target:
            lo, f_lo = mid, f_mid
        else:
            hi, f_hi = mid, f_mid
    m_final = 0.5*(lo + hi)
    f_final = run_once(m_final, pos, "final", fine_np, keep=True)
    return m_final, f_final

# -------- main --------
if __name__ == "__main__":
    # Targets from Fig. 9 captions (RBETA @ 228 GHz)
    m0, f0 = tune(target=0.53, pos=0)
    print(f"pos0: M_unit ≈ {m0:.3e}, Ftot ≈ {f0:.3f} Jy  ->  {OUTDIR}/fig9_rbeta_pos0_final_320px.h5")

    m1, f1 = tune(target=0.52, pos=1)
    print(f"pos1: M_unit ≈ {m1:.3e}, Ftot ≈ {f1:.3f} Jy  ->  {OUTDIR}/fig9_rbeta_pos1_final_320px.h5")
