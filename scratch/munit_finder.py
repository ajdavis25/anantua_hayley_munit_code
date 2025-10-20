import sys, io, os, re, argparse, subprocess, math
import numpy as np, pandas as pd
from ipole_many_models import runIPOLE

# ensure real-time output in SLURM logs
try:
    sys.stdout.reconfigure(line_buffering=True)
except AttributeError:
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, line_buffering=True)

# helper functions for bracket + bisection flux solving
def _run_ipole_once(simFile, nameBase, MunitUsed, fpos, ipole, inclination, Rhigh, freq_Hz, fov, npixel,
                    counterjet, rmax_geo, electronModel, sigma_transition, sigma_cut):
    """run ipole once and parse flux (Jy)"""
    print(f"bisect: runIPOLE(sim={simFile}, MunitUsed={MunitUsed:.3e}, fpos={fpos})", flush=True)
    result = subprocess.run(
        [
            'python', '-c',
            (
                "from ipole_many_models import runIPOLE; "
                f"runIPOLE('{simFile}', "
                f"'{nameBase.replace('.h5', '_{0:1.3f}.h5'.format(fpos))}', "
                f"{MunitUsed}, ipoleExecutable='{ipole}', thetacam={inclination}, "
                f"Rhigh={Rhigh}, freq_Hz={freq_Hz}, fov={fov}, npixel={npixel}, "
                f"counterjet={counterjet}, rmax_geo={rmax_geo}, positronRatio={fpos}, "
                "beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, "
                f"electronModel={electronModel}, sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"
            )
        ],
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
    )
    text = result.stdout + "\n" + result.stderr
    m = re.search(r'I,Q,U,V \[Jy\]:\s+([\d.\-eE]+)', text)
    if m:
        flux = float(m.group(1))
        print(f"bisect: parsed flux={flux:.6g} Jy (old format)", flush=True)
        return flux
    m = re.search(r'Ftot:\s+([\d.\-eE]+)\s+Jy', text)
    if m:
        flux = float(m.group(1))
        print(f"bisect: parsed flux={flux:.6g} Jy (Ftot format)", flush=True)
        return flux
    print("bisect: WARN could not parse flux; raw output follows\n" + text, flush=True)
    return None


def _bracket_and_bisect_to_flux(simFile, nameBase, fpos, target_flux, init_MunitUsed,
                                ipole, inclination, Rhigh, freq_Hz, fov, npixel,
                                counterjet, rmax_geo, electronModel, sigma_transition, sigma_cut,
                                max_bracket_steps=12, max_bisect_steps=20, tol=0.05):
    """find MunitUsed giving ~target_flux Jy by geometric bracketing + bisection"""
    f0 = _run_ipole_once(simFile, nameBase, init_MunitUsed, fpos, ipole, inclination, Rhigh,
                         freq_Hz, fov, npixel, counterjet, rmax_geo,
                         electronModel, sigma_transition, sigma_cut)
    if f0 is None:
        return None, None
    if abs(f0 - target_flux) <= tol:
        return init_MunitUsed, f0

    lo_M, hi_M = None, None
    lo_F, hi_F = None, None
    M, F = init_MunitUsed, f0
    if F < target_flux:
        for _ in range(max_bracket_steps):
            M_up = M * 10.0
            F_up = _run_ipole_once(simFile, nameBase, M_up, fpos, ipole, inclination, Rhigh,
                                   freq_Hz, fov, npixel, counterjet, rmax_geo,
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
            F_dn = _run_ipole_once(simFile, nameBase, M_dn, fpos, ipole, inclination, Rhigh,
                                   freq_Hz, fov, npixel, counterjet, rmax_geo,
                                   electronModel, sigma_transition, sigma_cut)
            if F_dn is None:
                return None, None
            if F_dn <= target_flux:
                lo_M, lo_F, hi_M, hi_F = M_dn, F_dn, M, F
                break
            M, F = M_dn, F_dn

    if lo_M is None or hi_M is None:
        print("bisect: FAIL could not bracket target flux.", flush=True)
        return None, None

    for _ in range(max_bisect_steps):
        mid_M = math.sqrt(lo_M * hi_M)
        mid_F = _run_ipole_once(simFile, nameBase, mid_M, fpos, ipole, inclination, Rhigh,
                                freq_Hz, fov, npixel, counterjet, rmax_geo,
                                electronModel, sigma_transition, sigma_cut)
        if mid_F is None:
            return None, None
        if abs(mid_F - target_flux) <= tol:
            print(f"bisect: success M*={mid_M:.3e}, F*={mid_F:.4g}", flush=True)
            return mid_M, mid_F
        if mid_F < target_flux:
            lo_M, lo_F = mid_M, mid_F
        else:
            hi_M, hi_F = mid_M, mid_F

    if abs(lo_F - target_flux) < abs(hi_F - target_flux):
        return lo_M, lo_F
    else:
        return hi_M, hi_F


# main positron tuning function
def makePositronImages(
        simFile, Munit, MunitOffset_guess, MunitSlope_guess,
        positronRatios=np.linspace(0,1,2), nameBase=None,
        inclination=163.0, Rhigh=20, freq_Hz=228e9,
        ipole='/work/vmo703/aricarte/run_ipole.sh', fov=160.0,
        npixel=320, counterjet=0, rmax_geo=50, electronModel=2,
        sigma_transition=0.5, sigma_cut=2.0, row=0
):
    print(f"ROWSTART: row={row} simFile={simFile} Rhigh={Rhigh}", flush=True)

    # bracket & bisect for anchor solve
    target_flux = 0.5
    denom0, denom1 = (1 + 2*0.0), (1 + 2*1.0)
    init0 = MunitOffset_guess + MunitSlope_guess*Munit/denom0
    init1 = MunitOffset_guess + MunitSlope_guess*Munit/denom1
    print(f"bisect: anchors fpos0={init0:.3e}, fpos1={init1:.3e}", flush=True)

    M0_star, F0_star = _bracket_and_bisect_to_flux(
        simFile, nameBase, 0.0, target_flux, init0, ipole, inclination,
        Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
        electronModel, sigma_transition, sigma_cut)
    M1_star, F1_star = _bracket_and_bisect_to_flux(
        simFile, nameBase, 1.0, target_flux, init1, ipole, inclination,
        Rhigh, freq_Hz, fov, npixel, counterjet, rmax_geo,
        electronModel, sigma_transition, sigma_cut)

    if M0_star and M1_star:
        denom_diff = (1/denom0 - 1/denom1)
        B = (M0_star - M1_star) / (Munit * denom_diff)
        A = M0_star - B * Munit/denom0
        print(f"bisect: solved A={A:.3e}, B={B:.3e}", flush=True)

        pd.DataFrame({
            'Row': [row], 'MunitOffset_guess': [A], 'MunitSlope_guess': [B]
        }).to_csv(
            'converged_parameters.csv', mode='a',
            header=not os.path.exists('converged_parameters.csv'),
            index=False
        )
        print(f"saved (anchor-matched) parameters for row {row}", flush=True)
        return

    # fallback iterative tuner
    errors = [1.0]
    tolerance, max_iterations, iterations = 0.05, 100, 0

    while max(np.abs(np.array(errors))) > tolerance and iterations < max_iterations:
        fluxes = []
        for fpos in positronRatios:
            print(f"\ndebug: running row {row} with fpos: {fpos}", flush=True)
            MunitUsed = MunitOffset_guess + MunitSlope_guess*Munit / (1 + 2*fpos)
            result = subprocess.run(
                [
                    'python', '-c',
                    f"from ipole_many_models import runIPOLE; "
                    f"runIPOLE('{simFile}', "
                    f"'{nameBase.replace('.h5', '_{0:1.3f}.h5'.format(fpos))}', "
                    f"{MunitUsed}, ipoleExecutable='{ipole}', thetacam={inclination}, "
                    f"Rhigh={Rhigh}, freq_Hz={freq_Hz}, fov={fov}, npixel={npixel}, "
                    f"counterjet={counterjet}, rmax_geo={rmax_geo}, positronRatio={fpos}, "
                    "beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, "
                    f"electronModel={electronModel}, sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"
                ],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True
            )
            text = result.stdout + '\n' + result.stderr
            flux = None
            m = re.search(r'I,Q,U,V \[Jy\]:\s+([\d.\-eE]+)', text)
            if m:
                flux = float(m.group(1))
            if flux is None:
                m = re.search(r'Ftot:\s+([\d.\-eE]+)\s+Jy', text)
                if m:
                    flux = float(m.group(1))
            if flux:
                print(f"debug: fpos={fpos}, flux={flux}", flush=True)
                fluxes.append(flux)
            else:
                print(f"warn: could not parse flux for fpos={fpos}", flush=True)

        errors = [(f - 0.5) for f in fluxes]
        print(f"errors: {errors}", flush=True)
        if not errors:
            break
        if max(np.abs(np.array(errors))) < tolerance:
            break
        if len(errors) >= 2:
            if errors[1] > errors[0]:
                new_MunitSlope_guess = MunitSlope_guess * 1.1 * abs(1 - errors[0])
                new_MunitOffset_guess = MunitOffset_guess * 0.9
            else:
                new_MunitSlope_guess = MunitSlope_guess * 0.9
                new_MunitOffset_guess = MunitOffset_guess * 1.1 * abs(1 - errors[1])
            MunitOffset_guess, MunitSlope_guess = new_MunitOffset_guess, new_MunitSlope_guess
        iterations += 1

    if iterations >= max_iterations:
        print(f"warn: row {row} hit iteration cap ({max_iterations}) without converging", flush=True)

    pd.DataFrame({
        'Row': [row], 'MunitOffset_guess': [MunitOffset_guess], 'MunitSlope_guess': [MunitSlope_guess]
    }).to_csv(
        'converged_parameters.csv', mode='a',
        header=not os.path.exists('converged_parameters.csv'), index=False
    )
    print(f"saved converged parameters after {iterations} iterations for row {row}", flush=True)


# main entry point
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ipole positron optimization')
    parser.add_argument('--start_row', '-start_row', type=int, default=None)
    args = parser.parse_args()

    start_row = int(os.getenv('SLURM_ARRAY_TASK_ID', args.start_row or 0))
    print(f"processing row {start_row}", flush=True)

    MunitOffset_guess, MunitSlope_guess = 5e24, 3.0
    params = pd.read_csv('/work/vmo703/data/test_example.csv', engine='python', encoding='utf-8-sig').rename(columns=str.strip)
    print("loaded columns:", list(params.columns), flush=True)

    model_dict = {'RBETA':[2,2.0],'RBETAwJET':[2,0.5],'CRITBETA':[4,2.0],'CRITBETAwJET':[4,0.5]}
    if start_row >= len(params):
        print(f"warn: row {start_row} exceeds csv length {len(params)}...exiting", flush=True)
        exit()

    row = params.iloc[start_row]
    timestep, evolution, model, spin, Rhigh = row['Timestep'], row['MAD/SANE'], row['Model'], row['Spin'], row['Rhigh']
    electronModel, sigma_cut = model_dict[model]
    positronRatios = np.linspace(0, row['positron frac'], int(row['positron frac']) + 1)
    print(f"positronRatios: {positronRatios}", flush=True)

    sim_prefix = 'Ma' if evolution == 'MAD' else 'Sa'
    Munit_val = 7.48745e+24 if evolution == 'MAD' else 1.82523e+27
    simFile = f'/work/vmo703/grmhd_dump_samples/{sim_prefix}{spin:+}_{timestep}.h5'
    nameBase = f'/work/vmo703/ipole_outputs/output_{sim_prefix}{spin:+}_{timestep}_model_{model}_Rhigh_{Rhigh}.h5'

    makePositronImages(simFile=simFile, Munit=Munit_val, MunitOffset_guess=MunitOffset_guess,
                       MunitSlope_guess=MunitSlope_guess, positronRatios=positronRatios,
                       nameBase=nameBase, Rhigh=Rhigh, electronModel=electronModel,
                       sigma_cut=sigma_cut, sigma_transition=2.0, row=start_row,
                       ipole='/work/vmo703/aricarte/run_ipole.sh')
