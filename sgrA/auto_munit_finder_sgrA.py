#!/usr/bin/env python3
"""
auto_munit_finder_sgrA.py
-------------------------
automatically tunes Munit scaling so IPOLE runs produce the target flux (~2.4 Jy)
for Sgr A* models, given a GRMHD dump and electron-heating prescription.

outputs a CSV with:
dump_index,timestep,MAD/SANE,model,spin,pos,MunitOffset,MunitSlope,Munit,MunitUsed,notes
"""

import sys, io, os, re, argparse, subprocess, math
import numpy as np, pandas as pd
from ipole_many_models import runIPOLE

# ensure real-time SLURM output
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, line_buffering=True)


# --------------------------------------------------------------------
# helper: run IPOLE once and extract flux (Jy), retry on failure
# --------------------------------------------------------------------
def _run_ipole_once(simFile, nameBase, MunitUsed, fpos, ipole, inclination,
                    Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
                    electronModel, sigma_transition, sigma_cut, retries=1):
    """run IPOLE once and parse flux (Jy) from stdout, retrying if parsing fails."""
    for attempt in range(retries + 1):
        print(f"[bisect] runIPOLE(sim={simFile}, MunitUsed={MunitUsed:.3e}, fpos={fpos}) (try {attempt+1}/{retries+1})", flush=True)
        result = subprocess.run(
            [
                "python", "-c",
                (
                    "from ipole_many_models import runIPOLE; "
                    f"runIPOLE('{simFile}', "
                    f"'{nameBase.replace('.h5', f'_{fpos:.1f}.h5')}', "
                    f"{MunitUsed}, ipoleExecutable='{ipole}', "
                    f"thetacam={inclination}, Rhigh={Rhigh}, freq_Hz={freq_Hz}, "
                    f"fov={fov}, npixel={npixel}, counterjet={counterjet}, "
                    f"rmax_geo={rmax_geo}, positronRatio={fpos}, "
                    "beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, "
                    f"electronModel={electronModel}, sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"
                )
            ],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
        )

        text = result.stdout + "\n" + result.stderr
        for pattern in [r"I,Q,U,V \[Jy\]:\s+([\d.\-eE]+)", r"Ftot:\s+([\d.\-eE]+)\s+Jy"]:
            m = re.search(pattern, text)
            if m:
                flux = float(m.group(1))
                print(f"[bisect] parsed flux={flux:.3f} Jy", flush=True)
                return flux

        print("[bisect] WARN: could not parse flux; raw output follows\n" + text, flush=True)

    return None


# --------------------------------------------------------------------
# geometric bracketing + bisection to match target flux
# --------------------------------------------------------------------
def _bracket_and_bisect_to_flux(simFile, nameBase, fpos, target_flux, init_MunitUsed,
                                ipole, inclination, Rhigh, freq_Hz, fov, npixel,
                                counterjet, rmax_geo, electronModel, sigma_transition, sigma_cut,
                                max_bracket_steps=12, max_bisect_steps=20, tol=0.05):
    """find MunitUsed that produces the target flux using geometric bracketing + bisection."""
    f0 = _run_ipole_once(simFile, nameBase, init_MunitUsed, fpos, ipole, inclination,
                         Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
                         electronModel, sigma_transition, sigma_cut)
    if f0 is None:
        return None, None

    if abs(f0 - target_flux) <= tol:
        return init_MunitUsed, f0

    lo_M, hi_M = None, None
    lo_F, hi_F = None, None
    M, F = init_MunitUsed, f0

    # geometric bracketing
    if F < target_flux:
        for _ in range(max_bracket_steps):
            M_up = M * 10.0
            F_up = _run_ipole_once(simFile, nameBase, M_up, fpos, ipole, inclination,
                                   Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
                                   electronModel, sigma_transition, sigma_cut)
            if F_up is None:
                return None, None
            if F_up >= target_flux:
                lo_M, lo_F, hi_M, hi_F = M, F, M_up, F_up
                break
            M, F = M_up, F_up
    else:
        for _ in range(max_bracket_steps):
            M_dn = M / 10.0
            if M_dn <= 0:
                break
            F_dn = _run_ipole_once(simFile, nameBase, M_dn, fpos, ipole, inclination,
                                   Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
                                   electronModel, sigma_transition, sigma_cut)
            if F_dn is None:
                return None, None
            if F_dn <= target_flux:
                lo_M, lo_F, hi_M, hi_F = M_dn, F_dn, M, F
                break
            M, F = M_dn, F_dn

    if lo_M is None or hi_M is None:
        print("[bisect] FAIL: could not bracket target flux.", flush=True)
        return None, None

    # bisection phase
    for _ in range(max_bisect_steps):
        mid_M = math.sqrt(lo_M * hi_M)
        mid_F = _run_ipole_once(simFile, nameBase, mid_M, fpos, ipole, inclination,
                                Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
                                electronModel, sigma_transition, sigma_cut)
        if mid_F is None:
            return None, None
        if abs(mid_F - target_flux) <= tol:
            print(f"[bisect] success M*={mid_M:.3e}, F*={mid_F:.3f} Jy", flush=True)
            return mid_M, mid_F
        if mid_F < target_flux:
            lo_M, lo_F = mid_M, mid_F
        else:
            hi_M, hi_F = mid_M, mid_F

    return (lo_M, lo_F) if abs(lo_F - target_flux) < abs(hi_F - target_flux) else (hi_M, hi_F)


# --------------------------------------------------------------------
# main Munit calibration routine
# --------------------------------------------------------------------
def makeMunitForSgrA(simFile, dump_index, timestep, evolution, model, spin,
                     Munit, MunitOffset_guess, MunitSlope_guess,
                     nameBase, inclination=45.0, Rhigh=20, freq_Hz=228e9,
                     ipole="/work/vmo703/aricarte/run_ipole.sh", fov=160.0,
                     npixel=320, counterjet=0, rmax_geo=50, electronModel=2,
                     sigma_transition=2.0, sigma_cut=2.0, row=0):
    """run Munit calibration for one simulation row."""
    print(f"[row start] row={row} simFile={simFile}", flush=True)
    target_flux = 2.4  # Jy
    fpos = 0.0         # no positrons for Sgr A*

    # skip missing dumps gracefully
    if not os.path.exists(simFile):
        print(f"[ERROR] missing dump file: {simFile}", flush=True)
        out_row = {
            "dump_index": dump_index, "timestep": timestep,
            "MAD/SANE": evolution, "model": model, "spin": spin, "pos": fpos,
            "MunitOffset": MunitOffset_guess, "MunitSlope": MunitSlope_guess,
            "Munit": Munit, "MunitUsed": np.nan, "notes": "missing_dump_file"
        }
        pd.DataFrame([out_row]).to_csv(
            "/work/vmo703/data/munits_results_sgrA.csv",
            mode="a", header=not os.path.exists("/work/vmo703/data/munits_results_sgrA.csv"), index=False
        )
        return

    # initial guess
    denom = (1 + 2 * fpos)
    init_MunitUsed = MunitOffset_guess + MunitSlope_guess * Munit / denom
    print(f"[bisect] initial guess MunitUsed={init_MunitUsed:.3e}", flush=True)

    # solve
    M_star, F_star = _bracket_and_bisect_to_flux(
        simFile, nameBase, fpos, target_flux, init_MunitUsed, ipole,
        inclination, Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
        electronModel, sigma_transition, sigma_cut
    )

    out_csv = "/work/vmo703/data/munits_results_sgrA.csv"
    if M_star is not None:
        out_row = {
            "dump_index": dump_index,
            "timestep": timestep,
            "MAD/SANE": evolution,
            "model": model,
            "spin": spin,
            "pos": fpos,
            "MunitOffset": MunitOffset_guess,
            "MunitSlope": MunitSlope_guess,
            "Munit": Munit,
            "MunitUsed": M_star,
            "notes": f"Matched {F_star:.2f} Jy @ {freq_Hz/1e9:.0f} GHz",
        }
    else:
        out_row = {
            "dump_index": dump_index,
            "timestep": timestep,
            "MAD/SANE": evolution,
            "model": model,
            "spin": spin,
            "pos": fpos,
            "MunitOffset": MunitOffset_guess,
            "MunitSlope": MunitSlope_guess,
            "Munit": Munit,
            "MunitUsed": np.nan,
            "notes": "failed_to_converge",
        }

    pd.DataFrame([out_row]).to_csv(
        out_csv, mode="a", header=not os.path.exists(out_csv), index=False
    )
    print(f"[save] wrote row {row} -> {out_csv}", flush=True)


# --------------------------------------------------------------------
# entry point
# --------------------------------------------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Auto-tune Munit scaling for Sgr A* IPOLE runs")
    parser.add_argument("--start_row", "-start_row", type=int, default=None)
    args = parser.parse_args()

    start_row = int(os.getenv("SLURM_ARRAY_TASK_ID", args.start_row or 0))
    print(f"[main] processing row {start_row}", flush=True)

    MunitOffset_guess, MunitSlope_guess = 5e24, 3.0
    params = pd.read_csv("/work/vmo703/data/sgra_batch.csv", engine="python", encoding="utf-8-sig").rename(columns=str.strip)
    if start_row >= len(params):
        print(f"[warn] row {start_row} exceeds csv length {len(params)} — exiting", flush=True)
        sys.exit()

    row = params.iloc[start_row]
    dump_index = int(row["Timestep"])
    timestep = dump_index * 5
    evolution, model, spin, Rhigh = row["MAD/SANE"], row["Model"], row["Spin"], row["Rhigh"]

    model_dict = {"RBETA": [2, 2.0], "RBETAwJET": [2, 1.0],
                  "CRITBETA": [4, 2.0], "CRITBETAwJET": [4, 1.0]}
    electronModel, sigma_cut = model_dict[model]

    sim_prefix = "Ma" if evolution == "MAD" else "Sa"
    Munit_val = 7.5e24 if evolution == "MAD" else 1.8e27
    spin_str = f"{float(spin):+g}"
    simFile = f"/work/vmo703/grmhd_dump_samples/{sim_prefix}{spin_str}_{dump_index}.h5"
    nameBase = f"/work/vmo703/ipole_outputs/sgrA/output_{sim_prefix}{spin_str}_{dump_index}_{model}.h5"

    makeMunitForSgrA(simFile, dump_index, timestep, evolution, model, spin,
                     Munit_val, MunitOffset_guess, MunitSlope_guess,
                     nameBase, Rhigh=Rhigh, electronModel=electronModel,
                     sigma_cut=sigma_cut, sigma_transition=2.0,
                     row=start_row, ipole="/work/vmo703/aricarte/run_ipole.sh")
