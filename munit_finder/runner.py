#from __future__ import annotations
"""
basic ipole runner — finds the binary, runs it, parses output for flux
"""
import os
import subprocess
import re
from typing import List

class IPoleRunner:
    #wrapper for running ipole on dumps and grabbing flux
    def __init__(self, ipole_path: str = None, search_dir: str = "~/Desktop/ipole_versions", generate_images: bool = False):
        #use given path or auto-search if not provided
        if ipole_path:
            #make sure path is valid and executable
            if not os.path.isfile(ipole_path) or not os.access(ipole_path, os.X_OK):
                raise FileNotFoundError(f"Specified ipole executable not found or not executable: {ipole_path}")
            self.ipole_path = ipole_path
        else:
            #expand ~ and look for ipole in tree
            search_dir = os.path.expanduser(search_dir)
            self.ipole_path = self._find_ipole_executable(search_dir)
        #clean up path and save flag
        self.ipole_path = os.path.abspath(self.ipole_path)
        self.generate_images = generate_images

    def _find_ipole_executable(self, search_dir: str) -> str:
        #crawl directory and grab first valid ipole executable
        if not os.path.isdir(search_dir):
            raise FileNotFoundError(f"ipole search directory not found: {search_dir}")
        for root, _, files in os.walk(search_dir):
            if "ipole" in files:
                ipole_exec = os.path.join(root, "ipole")
                if os.access(ipole_exec, os.X_OK):
                    return ipole_exec
        raise FileNotFoundError(f"No ipole executable found in directory {search_dir} or its subdirectories")

    def compute_flux(self, dump_file: str, Munit: float, thetacam: float,
                     freq_Hz: float = 230e9, MBH: float = 4.14e6, dsource: float = 8.127e3,
                     npixel_max: int = 129, npixel_min: int = 17, fov: float = 200.0,
                     rmax_geo: float = 50.0, phicam: float = 0.0, counterjet: int = 0,
                     target_nturns: int = -1, Rhigh: float = 20.0,
                     beta_crit_coefficient: float = 0.5, beta_crit: float = 1.0,
                     kappa: float = 3.5, emission_type: int = 4,
                     unpolarized: bool = True, print_command: bool = False) -> float:
        #build ipole command from args
        cmd = [self.ipole_path]
        if unpolarized:
            cmd.append("-unpol")
        if not self.generate_images:
            cmd.append("-quench")  #skip images unless asked

        #main args passed into cmd
        cmd.append(f"--dump={dump_file}")
        cmd.append(f"--thetacam={thetacam}")
        cmd.append(f"--phicam={phicam}")
        cmd.append(f"--freqcgs={freq_Hz}")
        cmd.append(f"--MBH={MBH:.6e}")
        cmd.append(f"--dsource={dsource:.6e}")
        cmd.append(f"--M_unit={Munit}")
        cmd.append(f"--trat_large={Rhigh}")
        cmd.append(f"--nx={npixel_max}"); cmd.append(f"--ny={npixel_max}")
        cmd.append(f"--fov={fov}")
        cmd.append(f"--counterjet={counterjet}")
        cmd.append(f"--rmax_geo={rmax_geo}")
        cmd.append(f"--target_nturns={target_nturns}")
        cmd.append(f"--positronRatio=0")
        cmd.append(f"--beta_crit={beta_crit}")
        cmd.append(f"--beta_crit_coefficient={beta_crit_coefficient}")
        cmd.append(f"--emission_type={emission_type}")
        cmd.append(f"--kappa_kappa={kappa}")
        if npixel_min is not None and npixel_min < npixel_max:
            cmd.append(f"--nx_min={npixel_min}"); cmd.append(f"--ny_min={npixel_min}")

        if print_command:
            print("Executing:", " ".join(cmd))

        #run ipole and grab output
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = proc.communicate()
        stdout_str = stdout.decode('utf-8', errors='replace')
        stderr_str = stderr.decode('utf-8', errors='replace')

        if proc.returncode != 0:
            #ipole died, dump stderr
            raise RuntimeError(f"ipole execution failed (return code {proc.returncode}). Error output:\n{stderr_str}")

        #scan for flux in both stdout and stderr
        flux = None
        for stream_name, stream in [('stdout', stdout_str), ('stderr', stderr_str)]:
            for line in stream.splitlines():
                if 'unpol' in line.lower():
                    print(f"[debug] scanning {stream_name} line: {line.strip()}")
                    match = re.search(
                        r"\(?\s*([-+]?[0-9]*\.?[0-9]+(?:[eE][+-]?[0-9]+)?)\s*Jy\s+unpol\s+xfer",
                        line, flags=re.IGNORECASE
                    )
                    if match:
                        try:
                            flux = float(match.group(1))
                            print(f"[debug] regex matched flux from {stream_name}: {flux:.5e}")
                            break
                        except ValueError:
                            print(f"[debug] regex float conversion failed on {stream_name}")
                            continue
            if flux is not None:
                break

        #fallback if no match — pull any float < 100
        if flux is None:
            print("[debug] flux not matched via regex — using numeric fallback...")
            numbers = re.findall(r"[-+]?\d*\.?\d+(?:[eE][+-]?\d+)?", stdout_str + '\n' + stderr_str)
            candidates = [float(n) for n in numbers if 0 < float(n) < 100]
            if candidates:
                flux = candidates[-1]
                print(f"[debug] fallback flux matched: {flux:.5e}")
            else:
                print("[debug] fallback found no usable numeric candidates")

        #if still no match, raise error
        if flux is None:
            raise RuntimeError(
                f"Failed to parse flux from ipole output. Output was:\n{stdout_str}\nStderr:\n{stderr_str}"
            )

        print(f"Computed flux = {flux:.3e} Jy for M_unit = {Munit:.3e} (dump {os.path.basename(dump_file)})")
        return flux

    def compute_average_flux(self, dump_files: List[str], Munit: float, thetacam: float, **kwargs) -> float:
        #run flux calc on all dumps and average
        fluxes = []
        for df in dump_files:
            flux = self.compute_flux(df, Munit, thetacam, **kwargs)
            fluxes.append(flux)
        if len(fluxes) == 0:
            raise RuntimeError("No dump files provided for flux calculation")
        avg_flux = sum(fluxes) / len(fluxes)
        print(f"Average flux = {avg_flux:.2f} Jy for M_unit = {Munit:.3e} over {len(fluxes)} dumps")
        return avg_flux

    def find_munit(self, dump_file: str, model: str = "RBETA", pos: int = 0,
                   thetacam: float = 50.0, flux_goal: float = 2.0,
                   Munit_guess: float = 1e19, **kwargs) -> float:
        """
        iteratively solve for M_unit that matches the flux_goal,
        using the optimizer's secant method.
        """
        from optimizer import optimize_munit
        # optimizer expects a list of dump files
        dumps = [str(dump_file)]

        # run the optimization loop
        M_trials, F_trials = optimize_munit(
            self, dumps, thetacam=thetacam,
            flux_goal=flux_goal,
            Munit_guess=Munit_guess,
            **kwargs
        )

        M_final = M_trials[-1]
        print(f"find_munit: dump={dump_file}, model={model}, pos={pos} -> Munit={M_final:.3e}")
        return M_final
    