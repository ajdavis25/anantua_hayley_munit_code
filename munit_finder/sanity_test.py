# quick check script
import numpy as np
from munit_oskey import offset_slope_key

# choose test case
timestep = "23000"
sim_key = "SANE+0.94"
Munit = 1.82523e27  # 1.82523e27 || 7.48745e24 depending on MAD/SANE

offsets = offset_slope_key[timestep]["MunitOffset"][sim_key]
slopes  = offset_slope_key[timestep]["MunitSlope"][sim_key]
fpos_values = [0.0, 0.25, 0.5, 1.0]

print(f"Testing {sim_key} at timestep {timestep}")
for i, (A, B) in enumerate(zip(offsets, slopes)):
    print(f"\nSet {i}: offset={A:.3e}, slope={B:.3e}")
    for f in fpos_values:
        Mused = A + B * Munit / (1 + 2*f)
        print(f"  fpos={f:.2f}  ->  Mused={Mused:.3e}")
