#!/usr/bin/env python3
import sys, os, importlib.util as ilu, subprocess, shlex

REPO_ROOT = "/work/vmo703/aricarte-copy/aricarte"
POS_DIR   = os.path.join(REPO_ROOT, "PositronIPOLEScripts")

# Make both the repo and the scripts dir importable
for p in (REPO_ROOT, POS_DIR):
    if p not in sys.path:
        sys.path.insert(0, p)

# Preload ipole_many_models as a module named exactly "ipole_many_models"
ipmm_path = os.path.join(POS_DIR, "ipole_many_models.py")
spec = ilu.spec_from_file_location("ipole_many_models", ipmm_path)
ipmm = ilu.module_from_spec(spec); spec.loader.exec_module(ipmm)
sys.modules["ipole_many_models"] = ipmm  # so unlocked script's "from ipole_many_models import ..." works

# Now import their unlocked driver
unlocked_path = os.path.join(POS_DIR, "ipole_positron_test_unlocked.py")
spec2 = ilu.spec_from_file_location("ipole_unlocked", unlocked_path)
unlocked = ilu.module_from_spec(spec2); spec2.loader.exec_module(unlocked)

# ---- Config for the Fig. 9 sanity check ----
sim   = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"  # t ≈ 25,000 M
Munit = 1.83e27                                           # base guess; use your table if you have one
out_base = "/work/vmo703/ipole_outputs/fig9_repro.h5"
ipole = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"

# Run for positronRatio 0 and 1 with their flags (rmax_geo, sigma_* are in the unlocked call)
unlocked.makePositronImages(sim, Munit,
                            positronRatios=[0,1],
                            nameBase=out_base, inclination=163.0, Rhigh=20,
                            ipole=ipole, fov=160.0, npixel=320, counterjet=0, rmax_geo=50)

# Report resulting total fluxes
def show_flux(fpath):
    try:
        import h5py
        import numpy as np
        with h5py.File(fpath, "r") as f:
            print(os.path.basename(fpath), "Ftot_unpol[Jy] =", f["/Ftot_unpol"][()])
    except Exception:
        cmd = f"h5dump -d /Ftot_unpol {shlex.quote(fpath)}"
        print("→", cmd); subprocess.run(cmd, shell=True, check=False)

for pr in (0.0, 1.0):
    show_flux(out_base.replace(".h5", f"_{pr:.3f}.h5"))
