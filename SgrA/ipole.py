import sys, pandas as pd
from pathlib import Path
from ipole_many_models import runIPOLE

"""
model presets: electron-heating parameterizations
RBETA -> continuous beta scaling
CRITBETA -> threshold beta critical
WJET variants = milder sigma-transition, allowing jet contribution
"""
model_settings = {
    "RBETA":        {"electronModel": 2, "sigma_transition": 2.0},
    "RBETAWJET":    {"electronModel": 2, "sigma_transition": 1.0},
    "CRITBETA":     {"electronModel": 4, "sigma_transition": 2.0},
    "CRITBETAWJET": {"electronModel": 4, "sigma_transition": 1.0},
}


def compute_munit_used(Munit, MunitOffset, MunitSlope, positron_frac):
    # new target flux for sgr a*: 2.4 Jy
    return MunitOffset + (MunitSlope * Munit) / (1.0 + 2.0 * positron_frac)


def run_from_csv(csv_path,
                 sim_dir="/work/vmo703/grmhd_dump_samples",
                 out_dir="/work/vmo703/ipole_outputs/sgrA",
                 ipole_exec="/work/vmo703/aricarte/run_ipole.sh",
                 inclination=20.0, # 20-90 degrees (sgr a* orientation uncertain; jetless)
                 Rhigh=20, # maybe test 10-40 later
                 freq_Hz=228e9, # or 86e9 for multiband test
                 fov=160.0, # in microarcseconds; ring ~50μas for sgr a*
                 npixel=320, # can i reduce to 256 for faster runs?
                 counterjet=0,
                 rmax_geo=50,
                 row_index=None):
    df = pd.read_csv(csv_path, skiprows=0)
    df.columns = df.columns.str.strip() # clean up spaces in headers

    if row_index is not None:
        row = df.iloc[int(row_index)]
        df = pd.DataFrame([row]) # single-row DataFrame

    for _, row in df.iterrows():
        if pd.isna(row["model"]) or pd.isna(row["Munit"]):
            continue

        # note: dump_index corresponds to timestep: dump_index = timestep / 5
        timestep = int(row["dump_index"])
        model = str(row["model"]).upper()
        spin = str(row["spin"])
        positron_frac = 0

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
            print(f"skipping missing file: {simFile}")
            continue

        nameBase = Path(out_dir) / f"{mad_sane}_spin{spin}_t{timestep}_{model}_pos{positron_frac}.h5"

        print(f"running {model} | timestep={timestep}, spin={spin}, pos={positron_frac}, MunitUsed={MunitUsed:.3e}")

        runIPOLE(
            str(simFile),
            str(nameBase),
            MunitUsed, # density scale -> return for 2.4 Jy
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
            sigma_cut=2.0, # maybe lower to 1.5 for sgr a* later (less jet)
        )


if __name__ == "__main__":
    # read SLURM array index if provided
    row_index = int(sys.argv[1]) if len(sys.argv) > 1 else None
    run_from_csv("/work/vmo703/data/munits_results.csv", row_index=row_index)
