import os, re, argparse, subprocess, numpy as np, pandas as pd

from ipole_many_models import runIPOLE


def makePositronImages(
        simFile,
        Munit,
        MunitOffset_guess,
        MunitSlope_guess,
        positronRatios=np.linspace(0,1,2),
        nameBase=None,
        inclination=163.0,
        Rhigh=20,
        freq_Hz=228e9,
        ipole='/work/vmo703/aricarte/run_ipole.sh',
        fov=160.0,
        npixel=320,
        counterjet=0,
        rmax_geo=50,
        electronModel=2,
        sigma_transition=0.5,
        sigma_cut=2.0,

        row=0
):
    errors = [1.0]
    tolerance = 0.05
    iterations = 0
    while max(np.abs(np.array(errors))) > tolerance: # we want the original errors array to keep sign, but to consider max of abs error
        fluxes = []
        # for a single dump file, make images for a variety of positron ratios. #rmax_geo=20 previous default
        for fpos in positronRatios:
            print(f"\n [debug] running row {row} with fpos: {fpos}")
            MunitUsed = MunitOffset_guess + MunitSlope_guess*Munit / (1. + 2.*fpos)

            # the positron ratio is put in the filename. MunitUsed = Munit / (1. + 2.*positronRatio); (FluxMax-FluxMin)/(listMax-listMin)=fpos
            # run ipole as subprocess to capture output
            print(f"running runIPOLE({simFile}, {MunitUsed}, ipoleExecutable={ipole}, thetacam={1.0*inclination+0.0*140}, Rhigh={Rhigh}, freq_Hz={freq_Hz}, fov={fov}, npixel={npixel}, counterjet={counterjet}, rmax_geo={rmax_geo}, positronRatio={fpos}, beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, electronModel={electronModel}, sigma_transition={sigma_transition}, sigma_cut={sigma_cut})")
            result = subprocess.run(
                [
                    'python',
                    '-c',
                    f"from ipole_many_models import runIPOLE; "
                    f"runIPOLE('{simFile}', "
                    f"'{nameBase.replace('.h5', '_{0:1.3f}.h5'.format(fpos))}', "
                    f"{MunitUsed}, ipoleExecutable='{ipole}', thetacam={1.0*inclination+0.0*140}, "
                    f"Rhigh={Rhigh}, freq_Hz={freq_Hz}, fov={fov}, npixel={npixel}, "
                    f"counterjet={counterjet}, rmax_geo={rmax_geo}, positronRatio={fpos}, "
                    f"beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, "
                    f"electronModel={electronModel}, sigma_transition={sigma_transition}, "
                    f"sigma_cut={sigma_cut})"
                ],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                universal_newlines=True
            )

            """result = subprocess.run(
                ['python', '-c', f"from ipole_many_models import runIPOLE; runIPOLE('{simFile}', '{nameBase.replace('.h5', '_{0:1.3f}.h5'.format(fpos))}', {MunitUsed}, ipoleExecutable='{ipole}', thetacam={1.0*inclination+0.0*140}, Rhigh={Rhigh}, freq_Hz={freq_Hz}, fov={fov}, npixel={npixel}, counterjet={counterjet}, rmax_geo={rmax_geo}, positronRatio={fpos}, beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, electronModel={electronModel}, sigma_transition={sigma_transition}, sigma_cut={sigma_cut})"],
                 capture_output=True,
                 text=True
            )"""

            # combine stdout + stderr for parsing
            output_text = result.stdout + '\n' + result.stderr
            # parse output for flux value
            flux = None
            
            # format 1: old ipole printout
            match = re.search(r'I,Q,U,V \[Jy\]:\s+([\d\.\-eE]+)', output_text)
            if match:
                flux = float(match.group(1))

            # format 2: ARC printout with "Ftot: ... Jy"
            if flux is None:
                match = re.search(r'Ftot:\s+([\d\.\-eE]+)\s+Jy', output_text)
                if match:
                    flux = float(match.group(1))

            if flux is not None:
                print(f'[debug] fpos: {fpos}, flux: {flux}')
                fluxes.append(flux)
            else:
                print(f'[warn] could not parse flux for fpos={fpos}')
                print("--- raw ipole output ---")
                print(output_text)

        errors = [(flux-0.5) for flux in fluxes]
        print(f"errors: {errors}")

        if not errors:
            # nothing parsed from ipole output
            print("warn: no flux values parsed -> errors list is empty")
            print("debug: ipole stdout:\n", result.stdout)
            print("debug: ipole stderr:\n", result.stderr)
            # bail out instead of crashing
            break

        if max(np.abs(np.array(errors))) < tolerance: # otherwise there will be new guesses
            break

        if len(errors) >= 2:
            if errors[1] > errors[0]:
                print("decreasing MunitOffset guess")
                new_MunitSlope_guess = MunitSlope_guess * 1.1 * abs(1 - errors[0])
                new_MunitOffset_guess = MunitOffset_guess * 0.9 
            elif errors[0] > errors[1]:
                print("decreasing MunitSlope guess")
                new_MunitSlope_guess = MunitSlope_guess * 0.9
                new_MunitOffset_guess = MunitOffset_guess * 1.1 * abs(1 - errors[1])
            # update guesses
            print(f"MunitOffset_guess: {MunitOffset_guess}->{new_MunitOffset_guess}, "
                  f"MunitSlope_guess: {MunitSlope_guess}->{new_MunitSlope_guess}")
            MunitOffset_guess = new_MunitOffset_guess
            MunitSlope_guess = new_MunitSlope_guess
        else:
            print("warn: not enough errors to update guesses (need at least 2)")
        '''else:
            new_MunitOffset_guess = MunitOffset_guess - 2e24*(0.3*errors[0] + 0.7*errors[1]) # weighted averages where slope mostly affects nonpositron. Won't exactly work for fpos > 1
            new_MunitSlope_guess = MunitSlope_guess - 1e0*(0.7*errors[0] + 0.3*errors[1])'''

        iterations += 1
        
    result_data = {'Row': [row],
                   'MunitOffset_guess': [MunitOffset_guess],
                   'MunitSlope_guess': [MunitSlope_guess]}
    result_df = pd.DataFrame(result_data)
    result_df.to_csv('converged_parameters.csv',
                     mode='a',
                     header = not os.path.exists('converged_parameters.csv'),
                     index=False)
    print(f"saved converged parameters after {iterations} iterations: MunitOffset_guess={MunitOffset_guess}, MunitSlope_guess={MunitSlope_guess}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run ipole positron optimization')
    parser.add_argument('--start_row', '-start_row', type=int, default=None,
                        help='Row index (default: use SLURM_ARRAY_TASK_ID if available)')
    args = parser.parse_args()

    if args.start_row is None:
        slurm_task_id = os.getenv('SLURM_ARRAY_TASK_ID')
        if slurm_task_id is not None:
            start_row = int(slurm_task_id)
        else:
            start_row = 0
    else:
        start_row = args.start_row

    print(f"processing row {start_row}")

    MunitOffset_guess = 5e24
    MunitSlope_guess = 3.0
    params = pd.read_csv(
        '/work/vmo703/data/example.csv',
        engine='python'
    )

    model_dict = {
        'RBETA': [2, 2.0],
        'RBETAwJET': [2, 0.5],
        'CRITBETA': [4, 2.0],
        'CRITBETAwJET': [4, 0.5]
    }


    if start_row >= len(params):
        print(f"warning: row {start_row} exceeds csv length {len(params)}...exiting")
        exit()

    row = params.iloc[start_row]
    timestep = row['Timestep']
    evolution = row['MAD/SANE']
    model = row['Model']
    spin = row['Spin']
    Rlow = 1
    Rhigh = row['Rhigh']
    inclination = 163
    sigma_transition = 2.0
    electronModel, sigma_cut = model_dict[model]
    positronRatios = np.linspace(0, row['positron frac'], row['positron frac']+1)
    print(f"positronRatios: {positronRatios}")

    if evolution == 'MAD':
        simFile = f'/work/vmo703/grmhd_dump_samples/Ma{spin:+}_{timestep}.h5'
        nameBase = f'/work/vmo703/ipole_outputs/output_Ma{spin:+}_{timestep}_model_{model}_Rhigh_{Rhigh}.h5' # if it's CritBeta, user should just ignore Rhigh
        makePositronImages(
            simFile=simFile,
            Munit=7.48745e+24,
            MunitOffset_guess=MunitOffset_guess,
            MunitSlope_guess=MunitSlope_guess,
            positronRatios=positronRatios,
            nameBase=nameBase,
            Rhigh=Rhigh,
            electronModel=electronModel,
            sigma_cut=sigma_cut,
            sigma_transition=sigma_transition,
            row=start_row,
            ipole='/work/vmo703/aricarte/run_ipole.sh'
        )
    elif evolution == 'SANE':
        simFile = f'/work/vmo703/grmhd_dump_samples/Sa{spin:+}_{timestep}.h5'
        nameBase = f'/work/vmo703/ipole_outputs/output_Sa{spin:+}_{timestep}_model_{model}_Rhigh_{Rhigh}.h5'
        makePositronImages(
            simFile=simFile, 
            Munit=1.82523e+27, 
            MunitOffset_guess=MunitOffset_guess, 
            MunitSlope_guess=MunitSlope_guess, 
            positronRatios=positronRatios, 
            nameBase=nameBase, 
            Rhigh=Rhigh, 
            electronModel=electronModel, 
            sigma_cut=sigma_cut, 
            sigma_transition=sigma_transition, 
            row=start_row,
            ipole='/work/vmo703/aricarte/run_ipole.sh'
        )

'''def makePositronImages(simFile, Munit, positronRatios=np.linspace(0,1,2), nameBase='/home1/09963/shaum/ipole_output/positron_test.h5', inclination=163.0, Rhigh=20, freq_Hz=228e9, \
    ipole='/home1/09963/shaum/ipole/ipole',fov=160.0, npixel=320, counterjet=0, rmax_geo=50):
#/work/epc964/ipole_HarvardCannonHPC_copy/aricarte/ipole+e-/ipole
#'/n/home13/anantua/ipole/aricarte/ipole+e-/ipole'

    #For a single dump file, make images for a variety of positron ratios. #rmax_geo=20 previous default
    for positronRatio in positronRatios:
        MunitOffset = 1.333e+25 #MADpt94:5e24,5e24,5e24,3.5e24 #MAD-pt5:1.3e25,8.e24,1.5e25,8.5e24 #SANE-pt5:2.5e28,8.e27,2e28,7.5e27 #Sim:R-Beta,RBetaBetaE,BetaC,BetaCBetaE  #1.2e+28 r50 1e28SaCrBetHyb  
        MunitSlope = 3.0e+0  #MADpt94:.19e1,.14e1,1.5e0,1.5e0 #MAD-pt5:.5e1,4.5e0,3.5e0,4.e0 #SANE-pt5:2.5e1,1.9e1,2.5e1,2.e1 #Sim:R-Beta,RBetaBetaE,BetaC,BetaCBetaE  #1.0e1  r50 1Pt5e1SaCrBetHyb
        MunitUsed = MunitOffset+ MunitSlope*Munit / (1. + 2.*positronRatio)#Optically Thin Approximation; 0.15positronRatio=ne+/(ne+ + ne-)? 
        #The positron ratio is put in the filename. MunitUsed = Munit / (1. + 2.*positronRatio); (FluxMax-FluxMin)/(listMax-listMin)=fpos
        runIPOLE(simFile, nameBase.replace('.h5', '_{0:1.3f}.h5'.format(positronRatio)), MunitUsed, ipoleExecutable=ipole, thetacam=1.0*inclination+0.0*140, Rhigh=Rhigh, \
                 freq_Hz=freq_Hz, fov=fov, npixel=npixel, counterjet=counterjet, rmax_geo=rmax_geo, positronRatio=positronRatio,beta_crit_coefficient=.5,constant_beta_e0=0.01,beta_crit=1,electronModel=2,sigma_transition=0.5,sigma_cut=2.0)#Can modify Munit by using MunitUsed; electron temp model by electronModel

if __name__ == '__main__':
    simFile ='/home1/09963/shaum/Ma+0.94_5600.h5'
             # '/work/epc964/PositronModelingGRMHD/sims/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/MAD/a+0.94/288x128x128_KHARMA_2041/torus.out0.05000.h5'
             #'/work/epc964/GWongBoxSims/Ma-0.5_torus.out0.03500.h5' #05400 #05420 #05440 #05460
             #'/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/MAD/a-0.5/288x128x128_KHARMA_2041/torus.out0.06000.h5'
             #'/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/MAD/a+0.94/384x192x192_IHARM/dump_00002000.h5'
    Munit = 7.48745e+24
    Rhigh = 20
    nameBase = '/home1/09963/shaum/ipole_output/output_Ma+0.94_5600.h5'
             # '/n/home13/anantua/ipole/aricarte/e+Output/positron_test_MADa-Pt5dump06000.h5' 
             # '/n/home13/anantua/ipole/aricarte/e+Output/positron_test_MAD.h5'

    positronRatios = np.linspace(0,1,2)#np.linspace(0,1,51)
    makePositronImages(simFile, Munit, positronRatios, Rhigh=Rhigh, nameBase=nameBase)#,freq_Hz=86e9,fov=320.0)

#    simFile = '/work/epc964/GWongBoxSims/Sa-0.5_4200.h5'
              # '/work/epc964/PositronModelingGRMHD/sims/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/SANE/a-0.5/288x128x128_KHARMA/torus.out0.05000.h5'
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/SANE/a-0.5/288x128x128_KHARMA/torus.out0.05000.h5'
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/SANE/a-0.5/288x128x128_KHARMA/torus.out0.06000.h5'
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/SANE/a-0.5/288x128x128_IHARM/dump_00001500.h5'
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/SANE/a-0.5/288x128x128_IHARM/dump_00001500.h5'
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/SANE/a+0.94/288x128x128_IHARM/dump_00001500.h5' 
              # '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/SANE/a+0.94/288x128x128_IHARM/dump_00002000.h5'
#    Munit = 1.82523e+27
#    Rhigh = 20
#    nameBase = '/work/epc964/ipole_HarvardCannonHPC_copy/aricarte/e+Output/positron_test_SANE_a-Pt5RBeta_dump4200.h5'
#'/work/epc964/ipole_HarvardCannonHPC_copy/aricarte/e+Output/positron_test_SANEa-Pt5RBetaBetaEPt01dump5000.h5'
             #'/n/home13/anantua/ipole/aricarte/e+Output/positron_test_SANEa-Pt5dump5000.h5'
             #'/n/home13/anantua/ipole/aricarte/e+Output/positron_test_SANEa-Pt5dump6000.h5'      
             #'/n/home13/anantua/ipole/aricarte/e+Output/positron_test_SANEtorus_SigmaTransition1_SigmaCut2.h5'
             #'../ipole_output/SANE_a+0.94_R20.0_positrons/positron_test.h5'

#    positronRatios = np.linspace(0,1,2)#np.linspace(0,1,51)
#    makePositronImages(simFile, Munit, positronRatios, Rhigh=Rhigh, nameBase=nameBase)#,freq_Hz=86e9,fov=320.0)
'''
