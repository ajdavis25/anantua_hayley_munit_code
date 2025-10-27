#!/usr/bin/env python3
import os, re, subprocess, pandas as pd

rows_to_redo = {5}  # e.g. the SANE -0.5, Rhigh=40 row
fpos_values = [0.0]  # regenerate the non-positron image(s); add 1.0 if desired
ipole = "/work/vmo703/aricarte/run_ipole.sh"
freq = 228e9
fov = 160.0
npix = 320
counterjet = 0
rmax_geo = 50
incl = 163.0
sigma_transition = 2.0

# use the same electron model mapping the converge code uses
MODEL_DICT = {'RBETA': (2, 2.0), 'RBETAwJET': (2, 0.5), 'CRITBETA': (4, 2.0), 'CRITBETAwJET': (4, 0.5)}

# load the audit (for sim paths) and the converged AB
audit = pd.read_csv("/work/vmo703/data/audit.csv").rename(columns=str.strip)
conv  = pd.read_csv("/work/vmo703/data/audit_subset.csv")

sel = conv[conv["Row"].isin(rows_to_redo)].merge(audit.reset_index().rename(columns={"index":"Row"}), on="Row", how="left")

def munit_for_evolution(evo):
    # match the main script
    return 7.48745e+24 if evo == "MAD" else 1.82523e+27

def run_ipole(simFile, out_h5, Mused, Rhigh, eModel, sCut, fpos):
    cmd = [
        "python","-c",
        (
            "from ipole_many_models import runIPOLE; "
            f"runIPOLE('{simFile}', '{out_h5}', {Mused}, "
            f"ipoleExecutable='{ipole}', thetacam={incl}, Rhigh={Rhigh}, "
            f"freq_Hz={freq}, fov={fov}, npixel={npix}, counterjet={counterjet}, "
            f"rmax_geo={rmax_geo}, positronRatio={fpos}, "
            "beta_crit_coefficient=0.5, constant_beta_e0=0.01, beta_crit=1, "
            f"electronModel={eModel}, sigma_transition={sigma_transition}, sigma_cut={sCut})"
        )
    ]
    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    text = res.stdout + "\n" + res.stderr
    m = re.search(r'I,Q,U,V \[Jy\]:\s+([\d.\-eE]+)', text) or re.search(r'Ftot:\s+([\d.\-eE]+)\s+Jy', text)
    flux = float(m.group(1)) if m else None
    print(f"[redo] fpos={fpos:.1f} Mused={Mused:.3e} -> Flux={flux} Jy")
    return flux

for _, r in sel.iterrows():
    row_id = int(r["Row"])
    evo = r["MAD/SANE"]
    spin = float(r["Spin"])
    t = int(r["Timestep"])
    model = r["Model"]
    Rhigh = int(r["Rhigh"])
    A = float(r["MunitOffset_guess"])
    B = float(r["MunitSlope_guess"])
    eModel, sCut = MODEL_DICT[model]
    Munit = munit_for_evolution(evo)
    sim_prefix = "Ma" if evo == "MAD" else "Sa"
    simFile = f"/work/vmo703/grmhd_dump_samples/{sim_prefix}{spin:+}_{t}.h5"

    # write to a clearly new filename to avoid mixing with old frames
    base = f"/work/vmo703/ipole_outputs/output_{sim_prefix}{spin:+}_{t}_model_{model}_Rhigh_{Rhigh}_REDO"
    os.makedirs(os.path.dirname(base), exist_ok=True)

    print(f"\n[redo] Row={row_id}  {simFile}  Rhigh={Rhigh}  A={A:.3e}  B={B:.3e}")
    for fpos in fpos_values:
        denom = (1.0 + 2.0*fpos)
        Mused = A + B*Munit/denom
        out_h5 = f"{base}_{fpos:0.3f}.h5"
        flux = run_ipole(simFile, out_h5, Mused, Rhigh, eModel, sCut, fpos)
