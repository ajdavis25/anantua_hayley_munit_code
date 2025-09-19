# /work/vmo703/scratch/paper_metrics.py
#!/usr/bin/env python3
import sys, h5py, numpy as np

def paper_metrics(path):
    with h5py.File(path, "r") as H:
        scale = float(H["/header/scale"][()])
        Ftot  = float(H["/Ftot"][()])                      # total polarized flux [Jy]
        pol   = np.array(H["/pol"]).transpose(1,0,2) * scale  # Jy/pixel

        It = pol[:,:,0].sum()
        Qt = pol[:,:,1].sum()
        Ut = pol[:,:,2].sum()
        Vt = pol[:,:,3].sum()

    # These are the quantities used in the Fig. 9 annotations
    PoverI = (Qt**2 + Ut**2)**0.5 / Ftot
    VoverI = Vt / Ftot
    return Ftot, PoverI, VoverI

def main(paths):
    for p in paths:
        F, P_I, V_I = paper_metrics(p)
        print("{}\n  Ftot={:.6f} Jy  P/I={:.3e}  V/I={:.3e}\n".format(p, F, P_I, V_I))

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python paper_metrics.py <file1.h5> [file2.h5 ...]")
        sys.exit(1)
    main(sys.argv[1:])
