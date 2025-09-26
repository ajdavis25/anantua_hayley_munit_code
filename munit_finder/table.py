import os
import csv
import numpy as np
from typing import List, Tuple
from optimizer import optimize_munit
from runner import IPoleRunner

def initialize_table(filename: str, formatting: str = 'Rhigh'):
    #make table if it doesn't exist, add header based on formatting
    if os.path.exists(filename):
        print(f"Table file '{filename}' already exists. Initialization skipped.")
        return
    os.makedirs(os.path.dirname(filename) or '.', exist_ok=True)
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if formatting == 'Rhigh':
            writer.writerow(['dumpPath', 'Munit', 'Rhigh', 'thetacam'])
        elif formatting == 'criticalBeta':
            writer.writerow(['dumpPath', 'Munit', 'beta_crit_coefficient', 'beta_crit', 'thetacam'])
        else:
            writer.writerow(['dumpPath', 'Munit', 'Param1', 'Param2', 'thetacam'])

def add_to_munit_table(runner: IPoleRunner, folder: str, thetacam: float, table_path: str,
                       flux_goal: float = 2.0, fractional_tolerance: float = 0.05, Munit_guess: float = 1e19,
                       dump_samples: int = 100, dump_range: Tuple[int, int] = (1000, 2000),
                       MBH: float = 4.14e6, dsource: float = 8.127e3, freq_Hz: float = 230e9,
                       npixel_max: int = 129, npixel_min: int = 17, fov: float = 200.0,
                       rmax_geo: float = 50.0, phicam: float = 0.0, counterjet: int = 0, target_nturns: int = -1,
                       beta_crit_coefficient: float = 0.5, beta_crit: float = 1.0,
                       kappa: float = 3.5, emission_type: int = 4, Rhigh: float = 20.0,
                       formatting: str = 'Rhigh', overwrite: bool = False, flip_retrograde: bool = True,
                       verbose: bool = False):
    #make sure folder ends with /
    if not folder.endswith(os.sep):
        folder = folder + os.sep

    #get all .h5 files in folder
    all_files = [f for f in os.listdir(folder) if f.endswith('.h5')]
    if not all_files:
        print(f"No HDF5 dump files found in folder: {folder}")
        return

    #extract dump number from file names
    dump_numbers = []
    for fname in all_files:
        try:
            if 'dump' in fname:
                num = int(fname.split('_')[1].split('.')[0])
            elif 'torus' in fname:
                num = int(fname.split('.')[2])
            elif fname.endswith('.h5'):
                num = int(fname.split('_')[-1].split('.')[0])
            else:
                continue
            dump_numbers.append((num, fname))
        except Exception:
            continue

    if not dump_numbers:
        print(f"No valid dump files found in folder: {folder}")
        return

    #clamp start/end to available range
    start, end = dump_range
    if end < start:
        start, end = end, start
    dump_numbers.sort(key=lambda x: x[0])
    available_nums = [num for num, _ in dump_numbers]
    start = max(start, available_nums[0])
    end = min(end, available_nums[-1])
    if start > end:
        print(f"Specified dump_range {dump_range} does not overlap with available dumps in {folder}")
        return

    #sample evenly from within the range
    count = min(dump_samples, end - start + 1)
    desired_numbers = np.linspace(start, end, count).astype(int)
    desired_numbers = np.unique(desired_numbers)
    file_map = {num: fname for num, fname in dump_numbers}
    selected_files = [os.path.join(folder, file_map[num]) for num in desired_numbers if num in file_map]
    if not selected_files:
        print(f"No dump files selected for folder {folder} in range {start}-{end}")
        return

    #use last dump in list as identifier
    id_file = selected_files[-1]

    #check if we've already done this combo
    if os.path.exists(table_path):
        existing_entries = []
        with open(table_path, 'r', newline='') as csvfile:
            reader = csv.reader(csvfile)
            next(reader, None)
            for row in reader:
                if not row:
                    continue
                existing_entries.append(row)

        combo_exists = False
        for entry in existing_entries:
            try:
                existing_file = entry[0]
                if formatting == 'Rhigh':
                    existing_Rhigh = float(entry[2]); existing_inc = float(entry[3])
                    if existing_file == id_file and abs(existing_Rhigh - Rhigh) < 1e-6 and abs(existing_inc - thetacam) < 1e-6:
                        combo_exists = True
                        break
                elif formatting == 'criticalBeta':
                    existing_bcoeff = float(entry[2]); existing_beta = float(entry[3]); existing_inc = float(entry[4])
                    if existing_file == id_file and abs(existing_bcoeff - beta_crit_coefficient) < 1e-6 and abs(existing_beta - beta_crit) < 1e-6 and abs(existing_inc - thetacam) < 1e-6:
                        combo_exists = True
                        break
            except Exception:
                continue

        if combo_exists:
            if overwrite:
                print("Warning: already in table, overwriting...")
            else:
                print("Skipping existing combo. Use overwrite=True to rerun.")
                return

    #flip theta if retrograde sim and user wants that
    if flip_retrograde and 'a-' in os.path.basename(os.path.normpath(folder)):
        thetacam = 180.0 - thetacam

    #params for ipole call
    ipole_params = {
        'MBH': MBH, 'dsource': dsource, 'freq_Hz': freq_Hz,
        'npixel_max': npixel_max, 'npixel_min': npixel_min, 'fov': fov,
        'rmax_geo': rmax_geo, 'phicam': phicam, 'counterjet': counterjet,
        'target_nturns': target_nturns, 'beta_crit_coefficient': beta_crit_coefficient,
        'beta_crit': beta_crit, 'kappa': kappa, 'emission_type': emission_type,
        'Rhigh': Rhigh
    }
    if verbose:
        ipole_params['print_command'] = True

    #run optimizer
    Munit_trials, flux_trials = optimize_munit(runner, selected_files, thetacam, flux_goal,
                                               fractional_tolerance=fractional_tolerance,
                                               Munit_guess=Munit_guess, **ipole_params)
    M_opt = Munit_trials[-1]

    #write result to csv
    with open(table_path, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        if formatting == 'Rhigh':
            row = [id_file, f"{M_opt}", f"{Rhigh:.2f}", f"{thetacam:.2f}"]
        elif formatting == 'criticalBeta':
            row = [id_file, f"{M_opt}", f"{beta_crit_coefficient:.2f}", f"{beta_crit:.2f}", f"{thetacam:.2f}"]
        else:
            row = [id_file, f"{M_opt}", "--", "--", f"{thetacam:.2f}"]
        writer.writerow(row)
        print(f"Added entry: {','.join(row)}")

def create_munit_table(runner: IPoleRunner, output_path: str, folders: List[str], inclination_list: List[float],
                       dump_samples: int = 100, dump_range: Tuple[int, int] = (1000, 2000),
                       MBH: float = 4.14e6, dsource: float = 8.127e3, kappa: float = 3.5, emission_type: int = 4,
                       flux_goal: float = 2.0, fractional_tolerance: float = 0.05,
                       formatting: str = 'Rhigh', Rhigh_list: List[float] = None,
                       beta_crit_coefficient: float = 0.5, beta_crit: float = 1.0,
                       npixel_max: int = 129, npixel_min: int = 17, fov: float = 200.0,
                       rmax_geo: float = 50.0, phicam: float = 0.0, counterjet: int = 0,
                       target_nturns: int = -1, flip_retrograde: bool = True,
                       overwrite: bool = False, freq_Hz: float = 230e9, verbose: bool = False):
    #make csv with header
    initialize_table(output_path, formatting=formatting)

    for folder in folders:
        #guess M_unit based on sim type
        is_mad = 'MAD' in folder
        default_guess = 1e17 if is_mad else 1e20

        for incl in inclination_list:
            thetacam = incl
            if formatting == 'Rhigh':
                for Rhigh in (Rhigh_list or [20.0]):
                    add_to_munit_table(runner=runner, folder=folder, thetacam=thetacam,
                                       table_path=output_path, flux_goal=flux_goal,
                                       fractional_tolerance=fractional_tolerance,
                                       Munit_guess=default_guess, dump_samples=dump_samples,
                                       dump_range=dump_range, MBH=MBH, dsource=dsource, freq_Hz=freq_Hz,
                                       npixel_max=npixel_max, npixel_min=npixel_min, fov=fov,
                                       rmax_geo=rmax_geo, phicam=phicam, counterjet=counterjet,
                                       target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient,
                                       beta_crit=beta_crit, kappa=kappa, emission_type=emission_type,
                                       Rhigh=Rhigh, formatting=formatting, overwrite=overwrite,
                                       flip_retrograde=flip_retrograde, verbose=verbose)
            elif formatting == 'criticalBeta':
                add_to_munit_table(runner=runner, folder=folder, thetacam=thetacam,
                                   table_path=output_path, flux_goal=flux_goal,
                                   fractional_tolerance=fractional_tolerance,
                                   Munit_guess=default_guess, dump_samples=dump_samples,
                                   dump_range=dump_range, MBH=MBH, dsource=dsource, freq_Hz=freq_Hz,
                                   npixel_max=npixel_max, npixel_min=npixel_min, fov=fov,
                                   rmax_geo=rmax_geo, phicam=phicam, counterjet=counterjet,
                                   target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient,
                                   beta_crit=beta_crit, kappa=kappa, emission_type=emission_type,
                                   Rhigh=Rhigh_list[0] if Rhigh_list else 20.0,
                                   formatting=formatting, overwrite=overwrite,
                                   flip_retrograde=flip_retrograde, verbose=verbose)
