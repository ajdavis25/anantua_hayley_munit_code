#!/usr/bin/env python3
import csv, subprocess, textwrap, h5py, numpy as np, pathlib

CSV = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/munits_table.csv"
IP  = "/work/vmo703/aricarte-copy/aricarte/ipole+e-/ipole"
OUT = pathlib.Path("/work/vmo703/ipole_outputs")
OUT.mkdir(parents=True, exist_ok=True)

DUMP = "/work/vmo703/grmhd_dump_samples/Sa-0.5_5000.h5"
COMMON_PAR = textwrap.dedent(f"""
dump {DUMP}
freqcgs 2.28e11
thetacam 163
phicam 0
nx 320
ny 320
fovx_dsource 160
fovy_dsource 160
rcam 1000
counterjet 0

MBH 6.2e9
dsource 1.69e7

emission_type 4
electronModel RBETA
trat_small 1
trat_large 20
sigma_cut 2
""").strip()

def get_munit_used(csv_path, timestep, model, spin, pairs, flow="SANE"):
    with open(csv_path, newline="") as f:
        rows = list(csv.reader(f))
    for r in rows[3:]:                    # data start on line 4
        if not r or len(r) < 9:
            continue
        ts, flow_, model_, spin_, pairs_, *_rest, munit_used, _ = r
        try:
            ts_i = int(ts)
        except:
            continue
        if (ts_i == timestep and flow_.strip()==flow and model_.strip()==model
                and spin_.strip()==str(spin) and pairs_.strip()==str(pairs)):
            return r, float(munit_used.replace("E","e"))
    raise RuntimeError("Row not found in CSV")

def write_par(out_h5, munit, pos):
    txt = COMMON_PAR + f"\n\nM_unit {munit:.6E}\npositronRatio {pos}\noutfile {out_h5}\n"
    par = pathlib.Path(out_h5).with_suffix(".par")
    par.write_text(txt)
    return str(par)

def run_ipole(par):
    print(f"\n=== running: {IP} -par {par}")
    p = subprocess.run(
        [IP, "-par", par],
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        universal_newlines=True,   # Py3.6-safe
    )
    print(p.stdout)
    return p.stdout

def probe(h5path, fallback_pos=np.nan):
    with h5py.File(h5path, "r") as H:
        s = float(H["header"]["scale"][()])
        ftot    = float(H["Ftot"][()])           # ipole’s own total
        ftot_un = float(H["Ftot_unpol"][()])

        # images in "code" intensity units -> multiply by 'scale' to get Jy/pixel
        unpol = np.array(H["unpol"]).T * s
        pol   = np.array(H["pol"]).transpose(1,0,2) * s

        Iimg, Qimg, Uimg, Vimg = pol[:,:,0], pol[:,:,1], pol[:,:,2], pol[:,:,3]
        It, Qt, Ut, Vt = Iimg.sum(), Qimg.sum(), Uimg.sum(), Vimg.sum()
        Pinv = np.hypot(Qt, Ut)/It if It else np.nan
        Vinv = Vt/It                if It else np.nan

        nx = int(H["header"]["camera"]["nx"][()])
        ny = int(H["header"]["camera"]["ny"][()])
        freq = float(H["header"]["freqcgs"][()])

        # Robust header reads (some builds omit these)
        try:
            pos  = float(H["header"]["electrons"]["positronRatio"][()])
        except Exception:
            pos = float(fallback_pos)

        try:
            em   = H["header"]["electrons"]["electronModel"][()]
            if hasattr(em, "decode"): em = em.decode()
        except Exception:
            em = "unknown"

    return {
        "nx":nx, "ny":ny, "freq":freq, "pos":pos, "model":em, "scale":s,
        "Ftot":ftot, "Ftot_unpol":ftot_un,
        "sum_unpol":float(unpol.sum()), "sum_I":float(Iimg.sum()),
        "It":It, "Qt":Qt, "Ut":Ut, "Vt":Vt, "Pinv":Pinv, "Vinv":Vinv,
    }

def pretty(label, row, res):
    print("\n==========================")
    print(label)
    print("CSV row:", row)
    print(f"electronModel={res['model']}  positronRatio={res['pos']}")
    print(f"nx={res['nx']} ny={res['ny']}  freq={res['freq']:.2e} Hz  scale={res['scale']:.6g}")
    print(f"Ftot (HDF5)         : {res['Ftot']:.6f}")
    print(f"Ftot_unpol (HDF5)   : {res['Ftot_unpol']:.6f}")
    print(f"sum(unpol)*scale    : {res['sum_unpol']:.6f}")
    print(f"sum(Iimg)*scale     : {res['sum_I']:.6f}")
    print("Stokes sums [Jy]    : I={It:.6f}  Q={Qt:.3e}  U={Ut:.3e}  V={Vt:.3e}".format(**res))
    print("Fractions           : P/I={Pinv:.3e}  V/I={Vinv:.3e}".format(**res))

# --- pos0 ---
row0, MU0 = get_munit_used(CSV, 5000, "RBETA", -0.5, 0, flow="SANE")
out0 = str(OUT / "fig9_rbeta_Sa-0.5_5000_pos0.h5")
par0 = write_par(out0, MU0, 0)
run_ipole(par0)
res0 = probe(out0, fallback_pos=0)

# --- pos1 ---
row1, MU1 = get_munit_used(CSV, 5000, "RBETA", -0.5, 1, flow="SANE")
out1 = str(OUT / "fig9_rbeta_Sa-0.5_5000_pos1.h5")
par1 = write_par(out1, MU1, 1)
run_ipole(par1)
res1 = probe(out1, fallback_pos=1)

pretty("POS0", row0, res0)
pretty("POS1", row1, res1)
