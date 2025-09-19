# /work/vmo703/scratch/fig9_numbers_only.py
import os, subprocess, h5py, numpy as np
from pathlib import Path

# ---- EDIT THESE ----
IP        = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
POS0_PAR  = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos0.par"
POS1_PAR  = "/work/vmo703/ipole_outputs/fig9_rbeta_Sa-0.5_5000_pos1.par"
RUN_IPOLE = False  # set True to re-run ipole each time
# --------------------

TARGETS = {
    "pos0": {"F_unpol": 0.53, "P_over_I": 3.93e-3, "V_over_I": -5.60e-3},
    "pos1": {"F_unpol": 0.52, "P_over_I": 1.62e-3, "V_over_I": -2.81e-3},
}

def parse_par(par_path):
    want = {
        "outfile": None, "dump": None, "freqcgs": None,
        "thetacam": None, "phicam": None, "nx": None, "ny": None,
        "fovx_dsource": None, "fovy_dsource": None, "rcam": None, "counterjet": None,
        "MBH": None, "dsource": None, "emission_type": None, "electronModel": None,
        "trat_small": None, "trat_large": None, "sigma_cut": None, "sigma_transition": None,
        "M_unit": None, "positronRatio": None, "beta_crit": None, "beta_crit_coefficient": None,
        "qu_conv": None,
    }
    with open(par_path) as f:
        for line in f:
            line=line.strip()
            if not line or line.startswith("#"): continue
            k,*vals = line.split()
            if k in want and vals:
                v = vals[0] if k in ["outfile","dump","electronModel"] else vals[0]
                try: v = float(v)
                except: pass
                want[k] = v
    return want

def run_ipole(par_path):
    cmd=[IP,"-par",par_path]
    print("\n=== running:"," ".join(cmd))
    p=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True)
    out,_=p.communicate()
    print(out)

def h5get(G, path, default=np.nan):
    """Safely fetch nested dataset/attr by path like 'header/camera/dx'."""
    cur=G
    for p in path.split("/"):
        if p=="":
            continue
        if isinstance(cur, h5py.Group) and p in cur:
            cur=cur[p]
        else:
            return default
    try:
        return cur[()] if hasattr(cur,"shape") or hasattr(cur,"dtype") else cur
    except Exception:
        return default

def probe_h5(h5_path):
    with h5py.File(h5_path,"r") as H:
        # camera & units actually used (all guarded)
        dx   = float(h5get(H,"/header/camera/dx"))
        dy   = float(h5get(H,"/header/camera/dy"))
        nx   = int(h5get(H,"/header/camera/nx", 0))
        ny   = int(h5get(H,"/header/camera/ny", 0))
        thetacam = float(h5get(H,"/header/camera/thetacam"))
        phicam   = float(h5get(H,"/header/camera/phicam"))
        fovx_d   = float(h5get(H,"/header/camera/fovx_dsource"))
        fovy_d   = float(h5get(H,"/header/camera/fovy_dsource"))
        rcam     = float(h5get(H,"/header/camera/rcam"))
        dsource  = float(h5get(H,"/header/dsource"))
        freq     = float(h5get(H,"/header/freqcgs"))
        scale    = float(h5get(H,"/header/scale"))
        Lunit    = float(h5get(H,"/header/units/L_unit"))
        Munit    = float(h5get(H,"/header/units/M_unit"))
        poshdr   = h5get(H,"/header/electrons/positronRatio", np.nan)

        Ftot       = float(h5get(H,"/Ftot", np.nan))
        Ftot_unpol = float(h5get(H,"/Ftot_unpol", np.nan))

        pol  = np.copy(H["pol"]).transpose((1,0,2)) * scale
        upol = np.copy(H["unpol"]).T * scale
        I,Q,U,V = pol[:,:,0], pol[:,:,1], pol[:,:,2], pol[:,:,3]

        pix = dx*dy*(Lunit/dsource*2.06265e11)**2/(nx*ny) if nx*ny>0 else np.nan

        Iint = I.sum()*pix
        Qint = Q.sum()*pix
        Uint = U.sum()*pix
        Vint = V.sum()*pix
        Pint = (Qint**2+Uint**2)**0.5
        PoverI = Pint/Iint if Iint!=0 else np.nan
        VoverI = Vint/Iint if Iint!=0 else np.nan

        return {
            "nx":nx,"ny":ny,"freq":freq,"thetacam":thetacam,"phicam":phicam,
            "fovx_dsource":fovx_d,"fovy_dsource":fovy_d,"rcam":rcam,"dsource":dsource,
            "L_unit":Lunit,"M_unit":Munit,"scale":scale,"positronRatio_header":poshdr,
            "Ftot":Ftot,"Ftot_unpol":Ftot_unpol,
            "sum_unpol_scale":float(upol.sum()),
            "Iint":Iint,"Qint":Qint,"Uint":Uint,"Vint":Vint,"Pint":Pint,
            "P_over_I":PoverI,"V_over_I":VoverI,
        }

def print_block(title, par, res, targets):
    print("\n" + "-"*28)
    print(title)
    print("PAR electronModel     :", par.get("electronModel"))
    print("PAR positronRatio     :", par.get("positronRatio"))
    print("PAR trat_small/large  :", par.get("trat_small"), par.get("trat_large"))
    print("PAR sigma_cut/trans   :", par.get("sigma_cut"), par.get("sigma_transition"))
    print("PAR beta_crit(_coef)  :", par.get("beta_crit"), par.get("beta_crit_coefficient"))
    print("PAR M_unit            :", "{:.3E}".format(par.get("M_unit")) if par.get("M_unit") else par.get("M_unit"))

    print("\n--- Header actually used ---")
    for k in ["nx","ny","freq","thetacam","phicam","fovx_dsource","fovy_dsource","rcam","dsource","L_unit","M_unit","scale","positronRatio_header"]:
        v=res[k]
        if isinstance(v,float):
            if k in ["freq","L_unit","M_unit"]:
                print(f"{k:>20}: {v:.6E}")
            else:
                print(f"{k:>20}: {v:.6g}")
        else:
            print(f"{k:>20}: {v}")

    print("\n--- Flux checks (Jy) ---")
    print("Ftot (HDF5)            :", f"{res['Ftot']:.6f}")
    print("Ftot_unpol (HDF5)      :", f"{res['Ftot_unpol']:.6f}")
    print("sum(unpol)*scale       :", f"{res['sum_unpol_scale']:.6f}", " (≈ Ftot_unpol)")

    print("\n--- Integrated Stokes (Jy) ---")
    print("I =", f"{res['Iint']:.6f}",
          "Q =", f"{res['Qint']:.3e}",
          "U =", f"{res['Uint']:.3e}",
          "V =", f"{res['Vint']:.3e}")
    print("Fractions  P/I =", f"{res['P_over_I']:.3e}",
          " V/I =", f"{res['V_over_I']:.3e}")

    print("\n--- Compare to Fig. 9 targets ---")
    print(f"Δ F_unpol = {res['Ftot_unpol']-targets['F_unpol']:+.3e}  (obs {res['Ftot_unpol']:.3f} vs {targets['F_unpol']:.3f})")
    print(f"Δ P/I     = {res['P_over_I']-targets['P_over_I']:+.3e}  (obs {res['P_over_I']:.3e} vs {targets['P_over_I']:.3e})")
    print(f"Δ V/I     = {res['V_over_I']-targets['V_over_I']:+.3e}  (obs {res['V_over_I']:.3e} vs {targets['V_over_I']:.3e})")

def main(par_path, tag):
    par=parse_par(par_path)
    if RUN_IPOLE:
        run_ipole(par_path)
    outf=str(par.get("outfile"))
    if not outf or not Path(outf).exists():
        raise SystemExit(f"Output file missing: {outf}")
    res=probe_h5(outf)
    target_key="pos0" if "pos0" in outf else "pos1"
    print_block(tag, par, res, TARGETS[target_key])

if __name__=="__main__":
    main(POS0_PAR,"POS0")
    main(POS1_PAR,"POS1")
