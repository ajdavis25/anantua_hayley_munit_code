#!/usr/bin/env python3


import os, re, csv, math, h5py, time, sys, fcntl, numpy as np
from dataclasses import dataclass
from pathlib import Path
from ipole_many_models import runIPOLE
from concurrent.futures import ProcessPoolExecutor, as_completed


@dataclass
class MunitParams:
    """
    a dataclass to store a set of physical and simulation parameters for a single
    simulation setup as read from the CSV file

    each instance of this class represents a unique combination of simulation inputs
    """
    timestep: int
    mad_sane: str
    model: str
    spin: float
    positron_frac: int
    munit_offset: float
    munit_slope: float
    munit: float
    munit_used: float
    spectrum_flag: str

    # aliases for brevity downstream
    @property
    def A0(self): return self.munit_offset
    @property
    def S0(self): return self.munit_slope
    @property
    def M_base(self): return self.munit


def load_munit_params(csv_path):
    """
    loads munit parameters from csv file into a dictionary

    args:
        csv_path (str): path to the CSV file

    returns:
        dict: a dictionary where keys are tuples of (timestep, mad_sane, model, spin)
              and values are MunitParams objects, or an empty dictionary if
              the file can't be processed
    """
    params_dict = {}
    with open(csv_path, newline='') as csvfile:
        for _ in range(2): next(csvfile)
        reader = csv.DictReader(csvfile)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]

        for row in reader:
            row = {k.strip(): (v.strip() if v else "") for k, v in row.items()}
            if not row["Timestep"] or not row["MAD/SANE"] or not row["Model"] or not row["Spin"]:
                continue
            try:
                spin_str = row["Spin"].replace("(+)" , "").replace("+", "")
                key = (int(row["Timestep"]), row["MAD/SANE"].upper(), row["Model"].upper(), float(spin_str))
                if int(row["positron frac"]) != 0:
                    continue
                params = MunitParams(
                    timestep=int(row["Timestep"]),
                    mad_sane=row["MAD/SANE"].upper(),
                    model=row["Model"].upper(),
                    spin=float(spin_str),
                    positron_frac=int(row["positron frac"]),
                    munit_offset=float(row["MunitOffset"]),
                    munit_slope=float(row["MunitSlope"]),
                    munit=float(row["Munit"]),
                    munit_used=float(row["MunitUsed"]),
                    spectrum_flag=row["M87 GRMONTY spectrum?"]
                )
                params_dict[key] = params
            except (KeyError, ValueError) as e:
                print(f"skipping row due to error: {e}\nRow: {row}")
    return params_dict


# --- CONSTANT SETTINGS ---
ipole = "/work/vmo703/aricarte/run_ipole.sh"
dump_dir = "/work/vmo703/grmhd_dump_samples/"
output_dir = "/work/vmo703/ipole_outputs/"
csv_path = "/work/vmo703/data/munits_table.csv"
summary_path = Path("/work/vmo703/data/summary_table.csv")

Rhigh = 20
sigma_cut = 2.0
thetacam = 163.0
rmax_geo = 50
freq_Hz = 228e9
fov = 160.0
npixel = 320
counterjet = 0
T0, T1 = 0.50, 0.50

model_settings = {
    "RBETA":        {"electronModel": 2, "sigma_transition": 2.0},
    "RBETAWJET":    {"electronModel": 2, "sigma_transition": 1.0},
    "CRITBETA":     {"electronModel": 4, "sigma_transition": 2.0},
    "CRITBETAWJET": {"electronModel": 4, "sigma_transition": 1.0},
}


def parse_sim_metadata(simFile):
    """
    parses metadata from filename

    this function expects filenames of the form:
        Ma+0.94_4000.h5
        Sa-0.50_2000.h5
        etc.
    
    args:
        simFile (str): path to the simulation file
    
    returns:
        tuple: a tuple containing the parsed metadata:
                (mad_sane (str), spin (float), timestep (int))

    raises:
        ValueError: if the filename does not match the expected pattern
    """
    basename = os.path.basename(simFile)
    match = re.match(r"([MS]a)([+-]?\d+\.\d+)_([0-9]+)", basename)
    if not match:
        raise ValueError(f"cannot parse MAD/SANE, spin, timestep from: {basename}")
    prefix, spin_str, timestep = match.groups()
    mad_sane = "MAD" if prefix == "Ma" else "SANE"
    return mad_sane, float(spin_str), int(timestep)


def read_flux(h5path):
    """
    reads total and polarized flux data from a HDF5 file
    this function opens an HDF5 file and extracts the total flux (F), and
    the polarized flux components (Q, U, V) if they exist using the h5py
    library to access the data

    args:
        h5path (str): path to the HDF5 file

    returns:
        tuple: a tuple containing the extracted flux values:
                (F: float, Q: float or None, U: float or None, V: float or None)
                returns None for a polarized component if it is not present in the file
    """
    with h5py.File(h5path, "r") as H:
        F = float(H["/Ftot"][()])
        Q = float(H["/Qt"][()]) if "/Qt" in H else None
        U = float(H["/Ut"][()]) if "/Ut" in H else None
        V = float(H["/Vt"][()]) if "/Vt" in H else None
    return F, Q, U, V


def run_once(r, Munit_used, simFile, nameBase, electronModel, sigma_transition):
    """
    executes a single run of the IPOLE simulation and reads the resulting flux
    this function constructs an output filename based on the input parameters,
    and then parses the output file to extract the total flux

    args:
        r (float): positron ratio (also used in output filename)
        Munit_used (float): Munit used in the simulation
        simFile (str): path to the simulation file
        nameBase (str): base name for the output HDF5 file
        electronModel (int): electron model to use in the simulation
                             (e.g. 2 or 4 for RBETA or CRITBETA respectively)
        sigma_transition (float): sigma transition value for the simulation
                                  (e.g. 1.0 for WJET models, 2.0 for standard models)

    returns:
        tuple: a tuple containing the total flux (F) and the output file
              (F: float, out: str)
    """
    out = nameBase.replace(".h5", f"_{int(r)}.h5")
    runIPOLE(simFile, out, Munit_used, ipoleExecutable=ipole,
             thetacam=thetacam, Rhigh=Rhigh, freq_Hz=freq_Hz,
             fov=fov, npixel=npixel, counterjet=counterjet,
             rmax_geo=rmax_geo, positronRatio=r,
             electronModel=electronModel, sigma_cut=sigma_cut,
             sigma_transition=sigma_transition)
    return read_flux(out)[0], out


def secant_for_flux(ratio, target, mu0, mu1, n_steps, simFile, nameBase, electronModel, sigma_transition):
    """
    finds the Munit value that produces a target flux using an iterative secant method
    this function performs a search in a total flux equal to the specified target value
    it uses an iterative approach, making successive calls to 'run_once' to refine the estimate

    args:
        ratio (float): positron ratio to use in the simulation
        target (float): desired target flux value to converge on
        mu0 (float): initial guess for Munit
        mu1 (float or None): an optional second guess for Munit, if None,
                             it is calculated based on the initial flux
        n_steps (int): number of secant iterations to perform
        simFile (str): path to the input simulation file
        nameBase (str): base filename for the output HDF5 file
        electronModel (int): electron model to use in the simulation
                             (e.g. 2 or 4 for RBETA or CRITBETA respectively)
        sigma_transition (float): sigma transition value for the simulation
                                  (e.g. 1.0 for WJET models, 2.0 for standard models)

    returns:
        tuple: a tuple containing the best estimate for the converged Munit and the
               corresponding flux values after the final iteration
               (Mcurr: float, Fcurr: float)
    """
    F0, _ = run_once(ratio, mu0, simFile, nameBase, electronModel, sigma_transition)
    if mu1 is None:
        mu1 = mu0 * (target / max(F0, 1e-30))
    F1, _ = run_once(ratio, mu1, simFile, nameBase, electronModel, sigma_transition)
    Mprev, Fprev = mu0, F0
    Mcurr, Fcurr = mu1, F1
    for _ in range(n_steps):
        denom = (Fcurr - Fprev)
        Mnext = Mcurr * (target / max(Fcurr, 1e-30)) if abs(denom) < 1e-12 else \
                Mcurr + (target - Fcurr) * (Mcurr - Mprev) / denom
        Fnext, _ = run_once(ratio, Mnext, simFile, nameBase, electronModel, sigma_transition)
        Mprev, Fprev, Mcurr, Fcurr = Mcurr, Fcurr, Mnext, Fnext
    return Mcurr, Fcurr


def process_dump_file(dump_file):
    simFile = os.path.join(dump_dir, dump_file)
    mad_sane, spin, timestep = parse_sim_metadata(simFile)
    rows = []
    t_start = time.time()

    for model_name, settings in model_settings.items():
        key = (timestep, mad_sane, model_name, spin)
        if key not in munit_params_dict:
            print(f"[{os.getpid()}] skipping {key} — not found in CSV")
            continue

        p = munit_params_dict[key]
        electronModel = settings["electronModel"]
        sigma_transition = settings["sigma_transition"]

        dump_subdir = os.path.join(output_dir, dump_file[:-3])
        os.makedirs(dump_subdir, exist_ok=True)
        os.makedirs(os.path.join(dump_subdir, "images"), exist_ok=True)
        nameBase = os.path.join(dump_subdir, f"{model_name}.h5")

        U0_guess = p.A0 + p.S0 * p.M_base
        U1_guess = p.A0 + p.S0 * p.M_base / 3.0
        U0_star, _ = secant_for_flux(0.0, T0, U0_guess, None, 1, simFile, nameBase, electronModel, sigma_transition)
        U1_star, _ = secant_for_flux(1.0, T1, U1_guess, None, 2, simFile, nameBase, electronModel, sigma_transition)
        S = 3.0 * (U0_star - U1_star) / (2.0 * p.M_base)
        A = U1_star - S * p.M_base / 3.0

        for r in [0.0, 1.0]:
            U_used = A + S * p.M_base / (1.0 + 2.0*r)
            F, out = run_once(r, U_used, simFile, nameBase, electronModel, sigma_transition)
            _, Q, Uq, V = read_flux(out)
            PoverI = (math.sqrt(Q**2 + Uq**2)/F) if (Q is not None and Uq is not None and F > 0) else None
            VoverI = (V/F) if (V is not None and F > 0) else None

            row = {
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
                "Elapsed (min)": round((time.time() - t_start) / 60, 2)
            }
            rows.append(row)

    return rows    


def append_to_summary(rows, path):
    with open(path, "a", newline="") as f:
        fcntl.flock(f, fcntl.LOCK_EX)
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        if f.tell() == 0:
            writer.writeheader()
        writer.writerows(rows)
        fcntl.flock(f, fcntl.LOCK_UN)


if __name__ == "__main__":
    munit_params_dict = load_munit_params(csv_path)
    all_dumps = sorted(f for f in os.listdir(dump_dir) if f.endswith(".h5"))

    if "SLURM_ARRAY_TASK_ID" not in os.environ:
        print("SLURM_ARRAY_TASK_ID not set. This script is for SLURM array jobs only")
        sys.exit(1)

    index = int(os.environ["SLURM_ARRAY_TASK_ID"])
    if index >= len(all_dumps):
        print(f"index {index} out of range for available dumps")
        sys.exit(1)

    dump_file = all_dumps[index]
    result_rows = process_dump_file(dump_file)

    if result_rows:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        append_to_summary(result_rows, summary_path)
        print(f"appended results for {dump_file} to {summary_path}")
