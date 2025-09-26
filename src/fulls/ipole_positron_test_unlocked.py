#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# script to run IPOLE radiative transfer code on GRMHD simulation snapshots
# with positron fraction consideration, based on parameters from a CSV file
import sys, logging, pandas as pd
from pathlib import Path
from typing import Optional
from ipole_many_models import runIPOLE


# setup logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(sys.stdout),
        # logging.FileHandler("ipole_batch.log") # optional file output
    ]
)

# global constants for directories and executable paths
SIM_DIR = "/work/vmo703/grmhd_dump_samples"
OUT_DIR = "/work/vmo703/ipole_outputs"
IPOLE_EXEC = "/work/vmo703/aricarte/run_ipole.sh"
CSV_PATH = "/work/vmo703/data/munits_table.csv"


# defines the settings for each model
model_settings = {
    "RBETA":        {"electronModel": 2, "sigma_transition": 2.0},
    "RBETAWJET":    {"electronModel": 2, "sigma_transition": 1.0},
    "CRITBETA":     {"electronModel": 4, "sigma_transition": 2.0},
    "CRITBETAWJET": {"electronModel": 4, "sigma_transition": 1.0},
}


def compute_munit_used(
        Munit: float,
        MunitOffset: float,
        MunitSlope: float,
        positron_frac: int
    ) -> float:
    """
    computes the 'MunitUsed' value, a scaling factor for the GRMHD simulation data
    this formula adjusts Munit based on an offset, a slope, and the positron fraction

    args:
        Munit (float): base mass unit from the input CSV
        MunitOffset (float): offset value for the Munit calculation
        MunitSlope (float): slope value for the Munit calculation
        positron_frac (int): the positron fraction, used to scale the slope term

    returns:
        float: the calculated MunitUsed value
    """
    return MunitOffset + (MunitSlope * Munit) / (1.0 + 2.0 * positron_frac)


def run_from_csv(
        csv_path: str,
        sim_dir: str = SIM_DIR,
        out_dir: str = OUT_DIR,
        ipole_exec: str = IPOLE_EXEC,
        inclination: float = 163.0,
        Rhigh: int = 20,
        freq_Hz: float = 228e9,
        fov: float = 160.0,
        npixel: int = 320,
        counterjet: int = 0,
        rmax_geo: int = 50,
        row_index: Optional[int] = None
    ) -> None:
    """
    main function to read a CSV file and run the IPOLE radiative transfer code for each
    specified GRMHD simulation snapshot

    args:
        csv_path (str): the path to the input CSV file containing simulation parameters
        sim_dir (str): the base directory for the GRMHD simulation data files
        out_dir (str): the directory where the IPOLE output files will be saved
        ipole_exec (str): the path to the IPOLE executable script
        inclination (float): the inclination angle of the observer in degrees
        Rhigh (float): the high radius value for the simulation
        freq_Hz (float): the observation frequency in Hertz
        fov (float): the field of view in microarcseconds
        npixel (int): the number of pixels for the output image
        counterjet (int): a flag (0 or 1) to indicate if the counterjet should be included
        rmax_geo (int): the maximum radius for the geodesic integration
        row_index (int, optional): if specified, only runs IPOLE for this single row
                                   of the CSV file. used for parallel jobs (e.g., SLURM)
    """
    df = pd.read_csv(csv_path, skiprows=2)
    df.columns = df.columns.str.strip() # clean up spaces in headers

    if row_index is not None:
        row = df.iloc[int(row_index)]
        df = pd.DataFrame([row]) # single-row DataFrame

    for _, row in df.iterrows():
        if pd.isna(row["Model"]) or pd.isna(row["Munit"]):
            continue

        timestep = int(row["Timestep"])
        model = str(row["Model"]).upper()
        spin = str(row["Spin"]).replace("(", "").replace(")", "")
        positron_frac = int(row["positron frac"])

        Munit = float(row["Munit"])
        MunitOffset = float(row["MunitOffset"])
        MunitSlope = float(row["MunitSlope"])
        MunitUsed = compute_munit_used(Munit, MunitOffset, MunitSlope, positron_frac)

        mad_sane = str(row["MAD/SANE"]).strip().upper()
        if mad_sane == "MAD":
            prefix = "Ma"
        elif mad_sane == "SANE":
            prefix = "Sa"
        else:
            raise ValueError(f"unknown MAD/SANE value: {mad_sane}")

        simFile = Path(sim_dir) / f"{prefix}{spin}_{timestep}.h5"
        if not simFile.exists():
            logging.warning(f"skipping missing file: {simFile}")
            continue

        nameBase = Path(out_dir) / f"{mad_sane}_spin{spin}_t{timestep}_{model}_pos{positron_frac}.h5"

        logging.info(f"running {model} | timestep={timestep}, spin={spin}, pos={positron_frac}, MunitUsed={MunitUsed:.3e}")

        try:
            runIPOLE(
                str(simFile),
                str(nameBase),
                MunitUsed,
                ipoleExecutable=ipole_exec,
                thetacam=inclination,
                Rhigh=Rhigh,
                freq_Hz=freq_Hz,
                fov=fov,
                npixel=npixel,
                counterjet=counterjet,
                rmax_geo=rmax_geo,
                positronRatio=positron_frac,
                beta_crit_coefficient=0.5,
                constant_beta_e0=0.01,
                beta_crit=1,
                electronModel=model_settings[model]["electronModel"],
                sigma_transition=model_settings[model]["sigma_transition"],
                sigma_cut=2.0,
            )
        except Exception as e:
            logging.error(f"error running IPOLE for {simFile.name}: {e}")


if __name__ == "__main__":
    # read SLURM array index if provided
    # if no argument is provided, it processes the entire CSV file
    row_index = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run_from_csv(CSV_PATH, row_index=row_index)
