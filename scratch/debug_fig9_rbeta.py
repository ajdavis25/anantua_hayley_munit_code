#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv, re, os, sys, subprocess, textwrap
from pathlib import Path
import h5py
import numpy as np

# ---------- paths you already have ----------
CSV      = Path("/work/vmo703/aricarte-copy/aricarte/ipole+e-/munits_table.csv")
DUMP     = Path("/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5")
IP       = Path("/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole")
OUTDIR   = Path("/work/vmo703/ipole_outputs")
OUTDIR.mkdir(parents=True, exist_ok=True)

# figure settings (paper)
FREQ = 2.28e11     # 228 GHz
THETA = 163.0
FOV = 160.0
NPX = 320
MBH = 6.2e9
DSOURCE = 16.9e6   # kpc -> ipole expects "kpc" units via this field

# --------------------------------------------------------------------------------
def read_munit_from_csv(csv_path: Path, timestep: int, model: str, spin: str, posfrac: int, flow: str):
    """
    Your CSV has 3 header/preamble rows; data starts at row 4.
    We read the 'MunitUsed' column directly for an exact match.
    """
    rows = []
    with open(csv_path, newline="") as f:
        r = csv.reader(f)
        # Skip the first 3 lines explicitly
        for _ in range(3):
            next(r, None)
        for row in r:
            # Expect exactly the schema you pasted
            # Timestep, MAD/SANE, Model, Spin, positron frac, MunitOffset, MunitSlope, Munit, MunitUsed, ...
            try:
                ts   = int(row[0])
                fam  = row[1].strip()
                mdl  = row[2].strip()
                spn  = row[3].strip()
                pf   = int(row[4])
                mu_u = row[8].strip()
            except Exception:
                continue
            if ts==timestep and fam==flow and mdl==model and spn==spin and pf==posfrac:
                rows.append((mu_u, row))

    print(f"\n[CSV] matched {len(rows)} row(s) for ts={timestep}, flow={flow}, model={model}, spin={spin}, pos={posfrac}")
    for mu_u, raw in rows:
        print("      raw:", raw)

    if not rows:
        raise RuntimeError("No matching row found in CSV. Check filters.")

    # Take the first match (there should be exactly one)
    mu_used_str = rows[0][0]
    # CSV uses scientific notation like 7.06E+28 — accept E/e
    if not re.match(r'^[0-9.+-Ee]+$', mu_used_str):
        raise RuntimeError(f"CSV M_unit parse failed: {mu_used_str}")
    mu_used = float(mu_used_str)
    print(f"[CSV] Using M_unit = {mu_used:.3e} (from MunitUsed column)")
    return mu_used

# --------------------------------------------------------------------------------
def write_par(par_path: Path, outfile: Path, m_unit: float, positron_ratio: int):
    """
    Only the essentials for RBETA; electronModel **numeric code 4** on this build.
    """
    par = textwrap.dedent(f"""\
        # Auto-generated for Fig. 9 RBETA sanity run
        dump {DUMP}
        outfile {outfile}

        freqcgs {FREQ:.6e}
        thetacam {THETA}
        phicam 0
        nx {NPX}
        ny {NPX}
        fovx_dsource {FOV}
        fovy_dsource {FOV}
        rcam 1000
        counterjet 0

        MBH {MBH:.6e}
        dsource {DSOURCE:.6e}

        emission_type 4
        electronModel 4      # RBETA (numeric code)
        trat_small 1
        trat_large 20
        sigma_cut 2

        M_unit {m_unit:.6e}
        positronRatio {positron_ratio}
    """)
    par_path.write_text(par)
    print(f"\n[PAR] wrote {par_path}")
    print("------ .par contents ------")
    print(par.rstrip())
    print("---------------------------")

# --------------------------------------------------------------------------------
def run_ipole(par_path: Path):
    """
    Run ipole; capture stdout for diagnostics (compatible with Python 3.6+).
    """
    cmd = [str(IP), "-par", str(par_path)]
    print(f"\n[RUN] {' '.join(cmd)}\n")
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    out_lines = []
    # stream output live
    for line in iter(p.stdout.readline, ''):
        sys.stdout.write(line)      # print live
        out_lines.append(line)
    p.wait()
    out = ''.join(out_lines)
    if p.returncode != 0:
        raise RuntimeError(f"ipole exited with code {p.returncode}")
    return out

# --------------------------------------------------------------------------------
def parse_ftot_from_stdout(stdout: str):
    m = re.search(r"Ftot:\s*([0-9.]+)\s*Jy\s*\(([^)]*)\)\s*scale=([0-9.eE+-]+)", stdout)
    if not m:
        print("[WARN] could not parse Ftot line from stdout.")
        return None, None
    ftot = float(m.group(1))
    scale = float(m.group(3))
    print(f"[PARSE] from stdout: Ftot = {ftot:.6f} Jy | scale = {scale:.6e}")
    return ftot, scale

# --------------------------------------------------------------------------------
def summarize_h5(h5_path: Path, expect_pos=None):
    """
    Validate the output numerically, independently of stdout.
    """
    print(f"\n[H5] opening {h5_path}")
    with h5py.File(h5_path, "r") as H:
        # header sanity
        dx = H["header"]["camera"]["dx"][()]
        dy = H["header"]["camera"]["dy"][()]
        lunit = H["header"]["units"]["L_unit"][()]
        dsource = H["header"]["dsource"][()]
        scale = H["header"]["scale"][()]
        nx = H["header"]["camera"]["nx"][()]
        ny = H["header"]["camera"]["ny"][()]
        freq = H["header"]["freqcgs"][()]
        try:
            munit_header = H["header"]["units"]["M_unit"][()]
        except Exception:
            munit_header = np.nan

        print(f"[H5] freq={freq:.3e} Hz | nx=ny={nx} | dx={dx:.6e} dy={dy:.6e} | L_unit={lunit:.6e} | dsource={dsource:.6e}")
        print(f"[H5] scale = {scale:.6e} | header M_unit = {munit_header:.3e}")

        # image arrays (transpose to match notebook convention)
        unpol = np.array(H["unpol"]).T * scale
        pol   = np.array(H["pol"]).T * scale  # (nx, ny, 4)

        # pixel size in (μas)^2 — matches the notebook logic
        pixelSize = dx * dy * (lunit / dsource * 2.06265e11)**2 / (pol.shape[0] * pol.shape[1])

        I = pol[:,:,0] / pixelSize
        Q = pol[:,:,1] / pixelSize
        U = pol[:,:,2] / pixelSize
        V = pol[:,:,3] / pixelSize

        Itot = I.sum() * pixelSize
        Qtot = Q.sum() * pixelSize
        Utot = U.sum() * pixelSize
        Vtot = V.sum() * pixelSize

        PoverI = np.hypot(Qtot, Utot) / Itot if Itot > 0 else np.nan
        VoverI = Vtot / Itot if Itot > 0 else np.nan

        # stored fields (if present)
        ftot_file = None
        for key in ("/Ftot", "/Ftot_unpol"):
            try:
                ftot_file = H[key][()]
                print(f"[H5] {key} = {ftot_file:.6f} Jy")
            except Exception:
                pass

        # n_pairs label (taken from filename suffix or header, if available)
        pos = None
        stem = h5_path.stem
        m = re.search(r'(?:^|_)pos(?P<val>[\d.]+)$', stem)
        if m: pos = float(m.group('val'))
        if pos is None and "electrons" in H["header"] and "positronRatio" in H["header"]["electrons"]:
            try:
                pos = float(H["header"]["electrons"]["positronRatio"][()])
            except Exception:
                pos = None
        print(f"[H5] inferred positronRatio = {pos}")

        print(f"[H5] Recomputed totals -> I={Itot:.6f} Jy | P/I={PoverI:.2e} | V/I={VoverI:.2e}")
        return dict(I=Itot, PoverI=PoverI, VoverI=VoverI, Ftot_field=ftot_file, pos=pos)

# --------------------------------------------------------------------------------
def main():
    # pull M_unit straight from your CSV (SANE, RBETA, spin -0.5, dump 5000)
    mu_pos0 = read_munit_from_csv(CSV, 5000, "RBETA", "-0.5", 0, "SANE")
    mu_pos1 = read_munit_from_csv(CSV, 5000, "RBETA", "-0.5", 1, "SANE")

    par0 = OUTDIR / "fig9_rbeta_Sa-0.5_5000_pos0.par"
    par1 = OUTDIR / "fig9_rbeta_Sa-0.5_5000_pos1.par"
    out0 = OUTDIR / "fig9_rbeta_Sa-0.5_5000_pos0.h5"
    out1 = OUTDIR / "fig9_rbeta_Sa-0.5_5000_pos1.h5"

    write_par(par0, out0, mu_pos0, 0)
    write_par(par1, out1, mu_pos1, 1)

    # run ipole and parse stdout for each
    out_txt0 = run_ipole(par0)
    ftot0, scale0 = parse_ftot_from_stdout(out_txt0)
    out_txt1 = run_ipole(par1)
    ftot1, scale1 = parse_ftot_from_stdout(out_txt1)

    # independently summarize from HDF5
    s0 = summarize_h5(out0, expect_pos=0)
    s1 = summarize_h5(out1, expect_pos=1)

    print("\n================ SUMMARY (compare to Fig. 9) ================")
    print("Target (paper): pos0 -> I≈0.53 Jy,  P/I≈3.93e-03, V/I≈-5.60e-03")
    print("                 pos1 -> I≈0.52 Jy,  P/I≈1.62e-03, V/I≈-2.81e-03")
    print("--------------------------------------------------------------")
    print(f"pos0  (M_unit {mu_pos0:.3e})  ipole Ftot={ftot0}  h5 I={s0['I']:.6f}  P/I={s0['PoverI']:.3e}  V/I={s0['VoverI']:.3e}  pos={s0['pos']}")
    print(f"pos1  (M_unit {mu_pos1:.3e})  ipole Ftot={ftot1}  h5 I={s1['I']:.6f}  P/I={s1['PoverI']:.3e}  V/I={s1['VoverI']:.3e}  pos={s1['pos']}")
    print("==============================================================\n")

if __name__ == "__main__":
    main()
