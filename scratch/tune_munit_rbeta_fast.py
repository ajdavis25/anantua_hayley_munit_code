import os, re, subprocess, tempfile, pathlib

IPOLE    = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
TEMPLATE = "/work/vmo703/scratch/runIPOLE_fig9_.par"
OUTDIR   = "/work/vmo703/ipole_outputs"
DUMP     = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"

def edit_par(base_par_text, **kv):
    txt = base_par_text
    def set_k(key, val):
        nonlocal txt
        if re.search(rf"(?m)^{key}\s+\S+", txt):
            txt = re.sub(rf"(?m)^{key}\s+\S+", f"{key} {val}", txt)
        else:
            txt += f"\n{key} {val}\n"
    for k, v in kv.items():
        set_k(k, v)
    return txt

def run_ipole_stream(par_path):
    # Stream stdout so you can see progress
    p = subprocess.Popen([IPOLE, "-par", str(par_path)],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True, bufsize=1)
    lines = []
    for line in p.stdout:
        print(line, end="")
        lines.append(line)
    p.wait()
    return p.returncode, "".join(lines)

def flux_from_stdout_or_h5(stdout_text, h5_path):
    m = re.search(r"Ftot:\s*([0-9.]+)\s*Jy", stdout_text)
    if m:
        return float(m.group(1))
    # fallback to h5dump
    p = subprocess.Popen(["h5dump", "-d", "/Ftot_unpol", h5_path],
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                         universal_newlines=True)
    out, _ = p.communicate()
    m2 = re.search(r"\(0\):\s*([0-9.eE+\-]+)", out or "")
    if m2:
        return float(m2.group(1))
    raise RuntimeError("Could not parse flux")

def run_once(munit, pos, tag, npix):
    base = pathlib.Path(TEMPLATE).read_text()
    outf = f"{OUTDIR}/fig9_rbeta_pos{pos}_{tag}_{npix}px.h5"
    txt = edit_par(
        base,
        dump=DUMP,
        outfile=outf,
        M_unit=f"{munit:.6e}",
        positronRatio=str(pos),
        electronModel="RBETA",
        freqcgs="2.28e11",
        nx=str(npix), ny=str(npix),
        fovx_dsource="160", fovy_dsource="160",
        thetacam="163",
        sigma_cut="2",
        trat_small="1", trat_large="20",
    )
    par = pathlib.Path(tempfile.mkstemp(suffix=".par")[1])
    par.write_text(txt)
    print(f"\n=== ipole run: pos={pos} npix={npix} M_unit={munit:.3e} → {outf}")
    rc, out = run_ipole_stream(par)
    if rc != 0:
        raise RuntimeError("ipole failed")
    return flux_from_stdout_or_h5(out, outf)

def tune(target, pos, coarse_iters=6, refine_iters=4):
    # Bracket (broad) at coarse resolution
    lo, hi = 1e27, 1e29
    f_lo = run_once(lo, pos, "lo", 160)
    f_hi = run_once(hi, pos, "hi", 160)
    while not (f_lo <= target <= f_hi):
        if target < f_lo:
            hi, f_hi = lo, f_lo
            lo /= 10.0
            f_lo = run_once(lo, pos, "loexp", 160)
        else:
            lo, f_lo = hi, f_hi
            hi *= 10.0
            f_hi = run_once(hi, pos, "hiexp", 160)
    # Coarse bisection
    for i in range(coarse_iters):
        mid = (lo + hi) / 2.0
        f_mid = run_once(mid, pos, f"mid{i+1}", 160)
        if f_mid < target: lo, f_lo = mid, f_mid
        else: hi, f_hi = mid, f_mid
    # Refine at full resolution
    for i in range(refine_iters):
        mid = (lo + hi) / 2.0
        f_mid = run_once(mid, pos, f"ref{i+1}", 320)
        if f_mid < target: lo, f_lo = mid, f_mid
        else: hi, f_hi = mid, f_mid
    final = (lo + hi) / 2.0
    f_final = run_once(final, pos, "final", 320)
    return final, f_final

if __name__ == "__main__":
    m0, f0 = tune(target=0.53, pos=0)
    print("\n>>> pos0: M_unit ≈ %.3e, Ftot ≈ %.3f Jy" % (m0, f0))
    m1, f1 = tune(target=0.52, pos=1)
    print("\n>>> pos1: M_unit ≈ %.3e, Ftot ≈ %.3f Jy" % (m1, f1))
