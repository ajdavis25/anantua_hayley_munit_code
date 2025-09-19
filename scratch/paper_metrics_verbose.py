#!/usr/bin/env python3
# /work/vmo703/scratch/paper_metrics_verbose.py
import sys, h5py, numpy as np, os

ARCSEC_PER_RAD = 206265.0
MUAS_PER_RAD   = ARCSEC_PER_RAD * 1e6

def probe(path):
    print(path)
    with h5py.File(path, "r") as H:
        # ——— header bits we care about
        scale   = float(H["/header/scale"][()])
        nx      = int(H["/header/camera/nx"][()])
        ny      = int(H["/header/camera/ny"][()])
        dx      = float(H["/header/camera/dx"][()])
        dy      = float(H["/header/camera/dy"][()])
        dsource = float(H["/header/dsource"][()])
        L_unit  = float(H["/header/units/L_unit"][()])

        # ——— datasets
        Ftot        = float(H["/Ftot"][()])          # polarized total [Jy]
        Ftot_unpol  = float(H["/Ftot_unpol"][()])    # unpolarized total [Jy]
        pol_raw     = np.array(H["/pol"])            # (nx,ny,4) or (ny,nx,4) – order doesn't matter for sums
        unpol_raw   = np.array(H["/unpol"])          # (nx,ny) or (ny,nx)

    # shapes/dtypes
    print("  nx,ny                :", nx, ny)
    print("  scale                :", scale)
    print("  dx,dy                : {:.6e} {:.6e}".format(dx,dy))
    print("  dsource [cm]         : {:.6e}".format(dsource))
    print("  L_unit  [cm]         : {:.6e}".format(L_unit))

    # ——— sums BEFORE applying scale (raw)
    # sums are invariant to transposes for totals
    sumI_raw = pol_raw[:,:,0].sum()
    sumQ_raw = pol_raw[:,:,1].sum()
    sumU_raw = pol_raw[:,:,2].sum()
    sumV_raw = pol_raw[:,:,3].sum()
    sumU_pol_raw = unpol_raw.sum()

    print("\n-- RAW sums (no scale) --")
    print("  sum I,Q,U,V (raw)    : {:.6e}  {:.6e}  {:.6e}  {:.6e}".format(
        sumI_raw, sumQ_raw, sumU_raw, sumV_raw))
    print("  sum unpol   (raw)    : {:.6e}".format(sumU_pol_raw))

    # ——— apply scale (converts raw -> Jy/pixel)
    sumI_scl = sumI_raw * scale
    sumQ_scl = sumQ_raw * scale
    sumU_scl = sumU_raw * scale
    sumV_scl = sumV_raw * scale
    sumU_pol_scl = sumU_pol_raw * scale

    print("\n-- AFTER scale (Jy/pixel summed over pixels) --")
    print("  sum I,Q,U,V (scaled) : {:.6e}  {:.6e}  {:.6e}  {:.6e}".format(
        sumI_scl, sumQ_scl, sumU_scl, sumV_scl))
    print("  sum unpol (scaled)   : {:.6e}".format(sumU_pol_scl))
    print("  Ftot (HDF5)          : {:.6e}".format(Ftot))
    print("  Ftot_unpol (HDF5)    : {:.6e}".format(Ftot_unpol))
    print("  diff sum(I)-Ftot     : {:+.3e}  (rel {:+.3e})".format(
        sumI_scl - Ftot, 0.0 if Ftot==0 else (sumI_scl-Ftot)/Ftot))
    print("  diff sum(unpol)-Funp : {:+.3e}  (rel {:+.3e})".format(
        sumU_pol_scl - Ftot_unpol, 0.0 if Ftot_unpol==0 else (sumU_pol_scl-Ftot_unpol)/Ftot_unpol))

    # ——— pixel-area factor (this is what often trips people up; for ipole outputs you DO NOT multiply by this)
    pix_area_factor = dx * dy * (L_unit / dsource * MUAS_PER_RAD)**2 / (nx * ny)
    print("\n-- Pixel geometry (diagnostic only; do NOT multiply ipole flux by this) --")
    print("  derived 'pixelSize'  : {:.6e}".format(pix_area_factor))
    wrong_I = sumI_scl * pix_area_factor
    print("  sum(I)*pixelSize     : {:.6e}  (WRONG for ipole .h5; just shown to debug unit confusions)".format(wrong_I))

    # ——— Fig. 9 “paper” metrics computed with *polarized* Ftot
    P_I = (sumQ_scl**2 + sumU_scl**2)**0.5 / Ftot if Ftot != 0 else np.nan
    V_I = sumV_scl / Ftot if Ftot != 0 else np.nan

    print("\n-- Paper-style metrics (using Ftot) --")
    print("  P/I                  : {:.3e}".format(P_I))
    print("  V/I                  : {:.3e}".format(V_I))

    # ——— Alternate: compute fractions dividing by sumI_scl (should be numerically identical to Ftot)
    P_I_alt = (sumQ_scl**2 + sumU_scl**2)**0.5 / sumI_scl if sumI_scl != 0 else np.nan
    V_I_alt = sumV_scl / sumI_scl if sumI_scl != 0 else np.nan
    print("  P/I (vs sumI)        : {:.3e}".format(P_I_alt))
    print("  V/I (vs sumI)        : {:.3e}".format(V_I_alt))

    # ——— Just for reference, if someone (incorrectly) used Ftot_unpol for the denominator
    P_over_Funp = (sumQ_scl**2 + sumU_scl**2)**0.5 / Ftot_unpol if Ftot_unpol != 0 else np.nan
    V_over_Funp = sumV_scl / Ftot_unpol if Ftot_unpol != 0 else np.nan
    print("\n-- Fractions if you (incorrectly) divide by Ftot_unpol --")
    print("  P/F_unpol            : {:.3e}".format(P_over_Funp))
    print("  V/F_unpol            : {:.3e}".format(V_over_Funp))
    print("-"*60 + "\n")

def main(paths):
    for p in paths:
        probe(p)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python paper_metrics_verbose.py <file1.h5> [file2.h5 ...]")
        sys.exit(1)
    main(sys.argv[1:])
