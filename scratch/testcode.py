#!/usr/bin/env python3


import os, re, csv, math, h5py, time, numpy as np
from dataclasses import dataclass
from pathlib import Path
from ipole_many_models import runIPOLE


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

    if not os.path.exists(out):
        runIPOLE(simFile, out, Munit_used, ipoleExecutable=ipole,
                thetacam=thetacam, Rhigh=Rhigh, freq_Hz=freq_Hz,
                fov=fov, npixel=npixel, counterjet=counterjet,
                rmax_geo=rmax_geo, positronRatio=r,
                electronModel=electronModel, sigma_cut=sigma_cut,
                sigma_transition=sigma_transition)
    else:
        print(f"skipping {out}, already exists.")

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


if __name__ == "__main__":
    """
    step by step execution of the full analysis pipeline:

    1. an empty dictionary 'results' is created to store the results of the analysis
    2. a dictionary of parameters for the "Munit" model is loaded from a CSV file
    3. a list of all simulation dump files ending with '.h5' in a specified directory is created
    4. a loop begins that iterates through each simulation dump file
    5. inside the loop, it extracts metadata like MAD/SANE type, spin, and timestep from the current dump file
    6. it then prints a message to indicate which file is being processed
    7. a nested loop begins that iterates through different electron heating models to be analyzed
    8. for each model, a unique 'key' is created using the timestep, MAD/SANE type, model name, and spin
    9. it checks if the 'key' exists in the 'munit_params_dict'; if not, it skips the current model and moves to the next
    10. it retrieves the model-specific parameters from the 'munit_params_dict' using the 'key'
    11. two initial 'guesses' for the 'Munit' variable are calculated based on the loaded parameters
    12. the code then uses a 'secant' method to find the optimal 'Munit' values, 'U0_star' and 'U1_star', for two different positron ratios
    13. it derives new 'Munit' parameters, 'A' (offset) and 'S' (slope), from the optimal 'U' values and prints them
    14. a third nested loop begins to run the analysis for two specific positron ratios (0.0 and 1.0)
    15. inside this loop, it calculates the 'Munit' value to be used for the current positron ratio
    16. it runs a simulation with the derived 'Munit' value and the current parameters
    17. the output of the simulation is read to extract flux values (I, Q, U, V)
    18. it calculates the polarization ratios P/I and V/I
    19. the calculated values (A, S, flux, polarization, etc.) are stored in the 'results' dictionary, organized by timestep, MAD/SANE, model, spin, and positron ratio
    20. it prints a verification message showing the calculated flux and polarization values for the current run
    21. after all models have been processed for a dump file, it prints the elapsed time
    22. after the main loop finishes, the code begins to prepare the results for export to a CSV file
    23. it iterates through the multi-level 'results' dictionary to flatten the data into a list of dictionaries (rows)
    24. finally, it checks if there are any rows to write and, if so, creates the output directory, writes the data to a CSV file with a header, and prints a success message.
    """
    results = {}
    munit_params_dict = load_munit_params(csv_path)
    all_dumps = sorted(f for f in os.listdir(dump_dir) if f.endswith(".h5"))

    for dump_file in all_dumps:
        simFile = os.path.join(dump_dir, dump_file)
        mad_sane, spin, timestep = parse_sim_metadata(simFile)

        print(f"starting timestep {timestep} ({mad_sane}, spin={spin})")
        t_start = time.time()

        for model_name, settings in model_settings.items():
            print(f"\n=== Running {model_name} on {dump_file} ===")
            key = (timestep, mad_sane, model_name, spin)
            if key not in munit_params_dict:
                print(f"skipping {key} — not found in CSV.")
                continue
            p = munit_params_dict[key]
            electronModel = settings["electronModel"]
            sigma_transition = settings["sigma_transition"]

            dump_subdir = os.path.join(output_dir, dump_file[:-3])
            os.makedirs(dump_subdir, exist_ok=True)

            images_subdir = os.path.join(dump_subdir, "images")
            os.makedirs(images_subdir, exist_ok=True)

            nameBase = os.path.join(dump_subdir, f"{model_name}.h5")

            U0_guess = p.A0 + p.S0 * p.M_base
            U1_guess = p.A0 + p.S0 * p.M_base / 3.0
            U0_star, _ = secant_for_flux(0.0, T0, U0_guess, None, 1, simFile, nameBase, electronModel, sigma_transition)
            U1_star, _ = secant_for_flux(1.0, T1, U1_guess, None, 2, simFile, nameBase, electronModel, sigma_transition)
            S = 3.0 * (U0_star - U1_star) / (2.0 * p.M_base)
            A = U1_star - S * p.M_base / 3.0
            print(f"derived MunitOffsetDerived(A) = {A:.3e}, MunitSlopeDerived(S) = {S:.3e}")

            for r in [0.0, 1.0]:
                U_used = A + S * p.M_base / (1.0 + 2.0*r)
                F, out = run_once(r, U_used, simFile, nameBase, electronModel, sigma_transition)
                _, Q, Uq, V = read_flux(out)
                PoverI = (math.sqrt(Q**2 + Uq**2)/F) if (Q and Uq and F > 0) else None
                VoverI = (V/F) if (V and F > 0) else None
                results.setdefault(timestep, {}) \
                       .setdefault(mad_sane, {}) \
                       .setdefault(model_name, {}) \
                       .setdefault(spin, {})[r] = {
                           "A": A, "S": S,
                           "Munit_base": p.M_base,
                           "Munit_used": U_used,
                           "flux_I": F,
                           "P/I": PoverI, "V/I": VoverI
                       }
                print(f"[verify r={int(r)}] Munit_used={U_used:.3e}  F={F:.6f}" +
                      (f"  P/I={PoverI:.3e}  V/I={VoverI:.3e}" if PoverI else ""))
            
            elapsed = time.time() - t_start
            print(f"finished timestep {timestep} in {elapsed/60:.1f} min")

    # --- CSV SUMMARY EXPORT ---
    rows = []
    for timestep, mad_dict in results.items():
        for mad_sane, model_dict in mad_dict.items():
            for model, spin_dict in model_dict.items():
                for spin, r_dict in spin_dict.items():
                    for r, vals in r_dict.items():
                        row = {
                            "Timestep": timestep,
                            "MAD/SANE": mad_sane,
                            "Model": model,
                            "Spin": spin,
                            "Positron Ratio": r,
                            "MunitOffsetDerived": vals["A"],
                            "MunitSlopeDerived": vals["S"],
                            "MunitBase": vals["Munit_base"],
                            "MunitUsed": vals["Munit_used"],
                            "Flux I": vals["flux_I"],
                            "P/I": vals["P/I"],
                            "V/I": vals["V/I"],
                            "Timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
                            "Elapsed (min)": round((time.time() - t_start) / 60, 2)
                        }
                        rows.append(row)

    if rows:
        summary_path.parent.mkdir(parents=True, exist_ok=True)
        with open(summary_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        print(f"wrote summary CSV to {summary_path}")
