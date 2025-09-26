"""
this modulee contains the numerical routine to optimize the M_unit parameter using a secant-like method
sourced form angelo at EHT base script
"""
import numpy as np
from typing import Tuple, List
from runner import IPoleRunner

def optimize_munit(runner: IPoleRunner, dump_files: List[str], thetacam: float, flux_goal: float,
                   fractional_tolerance: float = 0.05, Munit_guess: float = 1e19, max_iterations: int = 50,
                   **ipole_params) -> Tuple[List[float], List[float]]:
    """
    find m_unit that gives desired flux using secant method in log-log space
    """
    #set up first 2 guesses below and above the initial guess
    M0 = Munit_guess / np.sqrt(10)
    M1 = Munit_guess * np.sqrt(10)

    #run sim and get flux at those guesses
    F0 = runner.compute_average_flux(dump_files, M0, thetacam, **ipole_params)
    F1 = runner.compute_average_flux(dump_files, M1, thetacam, **ipole_params)

    print(f"[debug] F0 = {F0:.3e}, F1 = {F1:.3e}")

    Munit_trials = [M0, M1]
    flux_trials = [F0, F1]

    #warn if we’re not bracketed around the target
    if (F0 - flux_goal) * (F1 - flux_goal) > 0:
        print("Warning: Initial fluxes do not bracket the target flux. Convergence is not guaranteed.")

    #secant loop begins
    iteration = 0
    while iteration < max_iterations:
        iteration += 1

        #make sure flux is positive for log math
        if F0 <= 0 or F1 <= 0:
            raise ValueError("Flux values must be positive for secant method in log-space.")

        #get slope in log-log space
        alpha = np.log(F1 / F0) / np.log(M1 / M0)
        if alpha == 0:
            raise RuntimeError("Secant method failed: computed power-law exponent is zero.")

        #solve for new m using the power law
        A0 = F0 / (M0 ** alpha)
        M_new = np.exp(np.log(flux_goal / A0) / alpha)

        #run sim again at new m
        F_new = runner.compute_average_flux(dump_files, M_new, thetacam, **ipole_params)
        print(f"[debug] iteration {iteration}: M_new = {M_new:.3e}, F_new = {F_new:.3e}")

        #track trials
        Munit_trials.append(M_new)
        flux_trials.append(F_new)

        #update bracket depending on if we overshot
        if F_new < flux_goal:
            M0, F0 = M_new, F_new
        else:
            M1, F1 = M_new, F_new

        #check if we’re close enough to be done
        if abs(F_new - flux_goal) / flux_goal <= fractional_tolerance:
            print(f"Converged after {iteration + 1} iterations: M_unit ≈ {M_new:.3e}")
            break
    else:
        #bail if still not converged
        raise RuntimeError("Optimization did not converge within the maximum number of iterations")

    return Munit_trials, flux_trials
