import os, math, time, logging
from typing import Tuple, Optional, List, Dict, Union

from config import DUMPS, OUTPUT, IPATH, MODEL_SETTINGS, defaults
from flux_tools import read_flux
from ipole_many_models import runIPOLE  # assumes you retain this import path
from munit_loader import MunitParams

logger = logging.getLogger(__name__)


def parse_sim_metadata(simFile: str) -> Tuple[str, float, int]:
    """Extracts MAD/SANE, spin, timestep from filename like Ma+0.94_4000.h5"""
    basename = os.path.basename(simFile)
    import re
    match = re.match(r"([MS]a)([+-]?\d+\.\d+)_([0-9]+)", basename)
    if not match:
        raise ValueError(f"cannot parse MAD/SANE, spin, timestep from: {basename}")
    prefix, spin_str, timestep = match.groups()
    mad_sane = "MAD" if prefix == "Ma" else "SANE"
    return mad_sane, float(spin_str), int(timestep)


def run_once(
    ratio: float,
    Munit_used: float,
    simFile: str,
    nameBase: str,
    electronModel: int,
    sigma_transition: float,
) -> Tuple[float, str]:
    """Runs IPOLE once and returns total flux and output path"""
    out = nameBase.replace(".h5", f"_{int(ratio)}.h5")
    runIPOLE(
        simFile,
        out,
        Munit_used,
        ipoleExecutable=IPATH,
        thetacam=defaults.thetacam,
        Rhigh=defaults.Rhigh,
        freq_Hz=defaults.freq_Hz,
        fov=defaults.fov,
        npixel=defaults.npixel,
        counterjet=defaults.counterjet,
        rmax_geo=defaults.rmax_geo,
        positronRatio=ratio,
        electronModel=electronModel,
        sigma_cut=defaults.sigma_cut,
        sigma_transition=sigma_transition,
    )
    F, *_ = read_flux(out)
    return F, out


def secant_for_flux(
    ratio: float,
    target: float,
    mu0: float,
    mu1: Optional[float],
    n_steps: int,
    simFile: str,
    nameBase: str,
    electronModel: int,
    sigma_transition: float,
) -> Tuple[float, float]:
    """Finds Munit that yields target flux via secant method"""
    F0, _ = run_once(ratio, mu0, simFile, nameBase, electronModel, sigma_transition)
    if mu1 is None:
        mu1 = mu0 * (target / max(F0, 1e-30))

    F1, _ = run_once(ratio, mu1, simFile, nameBase, electronModel, sigma_transition)

    Mprev, Fprev = mu0, F0
    Mcurr, Fcurr = mu1, F1

    for _ in range(n_steps):
        denom = (Fcurr - Fprev)
        Mnext = (
            Mcurr * (target / max(Fcurr, 1e-30))
            if abs(denom) < 1e-12
            else Mcurr + (target - Fcurr) * (Mcurr - Mprev) / denom
        )
        Fnext, _ = run_once(ratio, Mnext, simFile, nameBase, electronModel, sigma_transition)
        Mprev, Fprev, Mcurr, Fcurr = Mcurr, Fcurr, Mnext, Fnext

    return Mcurr, Fcurr


def process_dump_file(dump_file: str, munit_params_dict: Dict[Tuple[int, str, str, float], MunitParams]) -> List[Dict[str, Union[int, str, float, None]]]:
    """Main driver for processing a dump file across models & positron ratios"""
    simFile = os.path.join(DUMPS, dump_file)
    mad_sane, spin, timestep = parse_sim_metadata(simFile)
    rows = []
    t_start = time.time()

    for model_name, settings in MODEL_SETTINGS.items():
        key = (timestep, mad_sane, model_name, spin)
        if key not in munit_params_dict:
            logger.warning(f"[{os.getpid()}] skipping {key} — not found in CSV")
            continue

        p = munit_params_dict[key]
        electronModel = settings["electronModel"]
        sigma_transition = settings["sigma_transition"]

        dump_subdir = os.path.join(OUTPUT, dump_file[:-3])
        os.makedirs(dump_subdir, exist_ok=True)
        os.makedirs(os.path.join(dump_subdir, "images"), exist_ok=True)
        nameBase = os.path.join(dump_subdir, f"{model_name}.h5")

        U0_guess = p.A0 + p.S0 * p.M_base
        U1_guess = p.A0 + p.S0 * p.M_base / 3.0

        U0_star, _ = secant_for_flux(0.0, defaults.T0, U0_guess, None, 1, simFile, nameBase, electronModel, sigma_transition)
        U1_star, _ = secant_for_flux(1.0, defaults.T1, U1_guess, None, 2, simFile, nameBase, electronModel, sigma_transition)

        S = 3.0 * (U0_star - U1_star) / (2.0 * p.M_base)
        A = U1_star - S * p.M_base / 3.0

        for r in [0.0, 1.0]:
            U_used = A + S * p.M_base / (1.0 + 2.0 * r)
            F, out = run_once(r, U_used, simFile, nameBase, electronModel, sigma_transition)
            _, Q, Uq, V = read_flux(out)

            PoverI = (math.sqrt(Q**2 + Uq**2) / F) if (Q is not None and Uq is not None and F != 0) else None
            VoverI = (V / F) if (V is not None and F != 0) else None

            rows.append({
                "Timestep": timestep,
                "MAD/SANE": mad_sane,
                "Model": model_name,
                "Spin": spin,
                "Positron Ratio": r,
                "MunitOffsetDerived": A,
                "MunitSlopeDerived": S,
                "MunitBase": p.M_base,
                "MunitUsed": U_used,
                "Flux I": F,
                "P/I": PoverI,
                "V/I": VoverI,
                "Timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                "Elapsed (min)": round((time.time() - t_start) / 60, 2),
            })

    return rows
