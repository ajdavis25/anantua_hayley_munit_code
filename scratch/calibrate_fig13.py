# calibrate_fig13.py
#!/usr/bin/env python3
import h5py, numpy as np, os


# knowns from csv for MAD, a=-0.5, RBETA, t=5000
A = 1.3e25 # MunitOffset
M = 7.49e24 # base Munit
S = 5.00 # MunitSlope
B = S*M


# targets (Jy)
T0, T1 = 0.55, 0.47 #r=0, r=1 positronRatio 0/1


# read measured flux & adjust if needed
p0 = '/work/vmo703/ipole_outputs/positron_test_Ma-0.5_dump5000_0.000.h5'
p1 = '/work/vmo703/ipole_outputs/positron_test_Ma-0.5_dump5000_1.000.h5'
with h5py.File(p0, 'r') as H: F0 = float(H['/Ftot'][()])
with h5py.File(p1, 'r') as H: F1 = float(H['/Ftot'][()])


print(f"measured Ftot: {F0:.6f} F(1)={F1:.6f} | targets: 0.55, 0.47")


# check simple common scaling
C0, C1 = T0/F0, T1/F1
print(f"common-scale factors: C0={C0:.6f} C1={C1:.6f} (rel. diff. {(C0-C1)/((C0+C1)/2):.1%})")


# solve for separate scale on offset (x) and slope (y)
# model: F0 approx. k(A+B), F1 approx k(A+B/3)
# after scaling: F0' = k(xA + yB) = T0 ; F1' = k(xA + yB/3) = T1
k = F0/(A+B)
y = (T0 - T1)/(F0 - F1)
x = ((T0 / k) - y*B)/A


new_offset = x*A
new_slope = y*B/M


print(f"\nsuggested:")
print(f"new offset: {new_offset:.3e}")
print(f"new slope:  {new_slope:.4e}")
