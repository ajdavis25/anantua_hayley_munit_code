#!/usr/bin/env python3
# Recreate Fig. 9 (SANE a=-0.5 at T=25000M) for M87
# - Fits M_unit against a flux goal using UNPOL/QUENCHED runs (fast, consistent)
# - Then renders full polarized images for the four cases:
#     {RBETA, CRITBETA} x {positronRatio=0, 1}
#
# Run:
#   PYTHONPATH=/work/vmo703/aricarte-copy/aricarte \
#   python /work/vmo703/scratch/repro_fig9_full.py

import os, re, tempfile, subprocess, pathlib
from pathlib import Path

# ---- Paths you already have ----
IPOLE = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
DUMP  = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"   # t = 25000M
OUTDIR = Path("/work/vmo703/ipole_outputs"); OUTDIR.mkdir(parents=True, exist_ok=True)

# ---- Figure-9 camera/e- model knobs (M87) ----
FREQ = "2.30e11"         # paper lists 230 GHz; if your numbers come out ~1–2% off, try 2.28e11
THETACAM = "163"
PHICAM   = "0"
FOV      = "160"
NPIX     = "320"
RMAX     = "50"
RCAM     = "1000"
COUNTERJET = "0"
MBH      = "6.2e9"
DSOURCE  = "16.9e6"

# Common temperature tuning used in the scripts
TRAT_SMALL = "1"
TRAT_LARGE = "20"
SIGMA_CUT  = "2"
BETA_CRIT_COEFF = "0.5"
BETA_CRIT  = "1"

# ---- Flux goals for M87 at 230 GHz (Fig. 9) ----
# If your exact PDF readback is slightly different, tweak these +/-0.01
GOAL_RBETA = {0: 0.53, 1: 0.52}
GOAL_CRIT  = {0: 0.53, 1: 0.52}

# ---- Helpers ----
def write_par(outfile, electron_model, positron_ratio, munit):
    """Write a strict .par capturing exactly the Fig. 9 choices."""
    lines = [
        f"dump {DUMP}",
        f"outfile {outfile}",
        f"MBH {MBH}",
        f"dsource {DSOURCE}",
        f"freqcgs {FREQ}",
        f"thetacam {THETACAM}",
        f"phicam {PHICAM}",
        f"nx {NPIX}",
        f"ny {NPIX}",
        f"fovx_dsource {FOV}",
        f"fovy_dsource {FOV}",
        f"rcam {RCAM}",
        f"counterjet {COUNTERJET}",
        f"rmax_geo {RMAX}",
        "",
        # electron temperature family
        f"electronModel {electron_model}",   # RBETA or CRITBETA (this ipole branch accepts the token)
        f"trat_small {TRAT_SMALL}",
        f"trat_large {TRAT_LARGE}",
        f"sigma_cut {SIGMA_CUT}",
        f"beta_crit_coefficient {BETA_CRIT_COEFF}",
        f"beta_crit {BETA_CRIT}",
        "",
        # normalization + pairs
        f"M_unit {munit:.6e}",
        f"positronRatio {positron_ratio}",
        "",
        # emission / polarization defaults
        "emission_type 4",  # thermal
    ]
    p = Path(tempfile.mkstemp(suffix=".par")[1])
    p.write_text("\n".join(lines))
    return p

def _run(cmd_list):
    # py3.6-compatible capture
    proc = subprocess.Popen(cmd_list, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    out = proc.communicate()[0].decode("utf-8", "replace")
    return proc.returncode, out

def _parse_unpol_flux(ipole_stdout):
    # ipole prints: "Ftot: X Jy (Y Jy unpol xfer)"
    m = re.search(r"\(.*?([\d.]+)\s+Jy\s+unpol xfer\)", ipole_stdout)
    if not m:
        raise RuntimeError("Could not parse unpolarized flux from ipole output.")
    return float(m.group(1))

def _parse_pol_flux(ipole_stdout):
    m = re.search(r"Ftot:\s*([\d.]+)\s*Jy\b", ipole_stdout)
    if not m:
        raise RuntimeError("Could not parse polarized Ftot from ipole output.")
    return float(m.group(1))

def measure_unpol_flux(parfile):
    rc, out = _run([IPOLE, "-par", str(parfile), "-unpol", "-quench"])
    if rc != 0:
        raise RuntimeError("ipole (unpol) failed:\n"+out)
    return _parse_unpol_flux(out)

def render_pol_image(parfile):
    rc, out = _run([IPOLE, "-par", str(parfile)])
    if rc != 0:
        raise RuntimeError("ipole (pol) failed:\n"+out)
    return _parse_pol_flux(out)

def tune_munit(target, electron_model, positron_ratio,
               m_lo=1e27, m_hi=1e29, max_iter=10, tag="RBETA"):
    """Robust bracket + bisection on UNPOL flux."""
    # expand bracket until it contains the target
    def flux_at(m, label):
        par = write_par("/dev/null", electron_model, positron_ratio, m)
        # give each probe a small on-disk tag if you want to inspect later:
        return measure_unpol_flux(par)

    f_lo = flux_at(m_lo, "lo")
    f_hi = flux_at(m_hi, "hi")

    expand_guard = 0
    while not (min(f_lo, f_hi) <= target <= max(f_lo, f_hi)):
        if target < min(f_lo, f_hi):
            m_lo /= 10.0; f_lo = flux_at(m_lo, "loexp")
        else:
            m_hi *= 10.0; f_hi = flux_at(m_hi, "hiexp")
        expand_guard += 1
        if expand_guard > 12:
            raise RuntimeError(f"Could not bracket target {target} Jy for {tag}, pos={positron_ratio}")

    # ensure increasing bracket (monotonic in this band)
    if f_lo > f_hi:
        m_lo, m_hi = m_hi, m_lo
        f_lo, f_hi = f_hi, f_lo

    m_mid = None; f_mid = None
    for _ in range(max_iter):
        m_mid = 0.5*(m_lo+m_hi)
        f_mid = flux_at(m_mid, "mid")
        if f_mid < target:
            m_lo, f_lo = m_mid, f_mid
        else:
            m_hi, f_hi = m_mid, f_mid
    return m_mid, f_mid

def do_one(mode, pos, target):
    assert mode in ("RBETA","CRITBETA")
    tag = f"{mode}_pos{pos}"
    print(f"\n=== Fitting {tag} to {target:.3f} Jy (UNPOL) ===")
    m_fit, f_fit = tune_munit(target, mode, pos, tag=mode)
    print(f"→ tuned M_unit ≈ {m_fit:.3e}, unpol F ≈ {f_fit:.3f} Jy")

    # final image (POL) to disk with that M_unit
    outfile = str(OUTDIR / f"fig9_{tag}.h5")
    par = write_par(outfile, mode, pos, m_fit)
    f_pol = render_pol_image(par)
    print(f"Rendered {outfile}    POL Ftot ≈ {f_pol:.3f} Jy")
    return outfile

def main():
    # RBETA (top row in paper’s caption) then CRITBETA
    for mode, goals in (("RBETA", GOAL_RBETA), ("CRITBETA", GOAL_CRIT)):
        for pos in (0, 1):
            do_one(mode, pos, goals[pos])

if __name__ == "__main__":
    main()
