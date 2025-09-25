import h5py
from typing import Tuple, Optional

def read_flux(h5path: str) -> Tuple[float, Optional[float], Optional[float], Optional[float]]:
    with h5py.File(h5path, "r") as H:
        F = float(H["/Ftot"][()])
        Q = float(H["/Qt"][()]) if "/Qt" in H else None
        U = float(H["/Ut"][()]) if "/Ut" in H else None
        V = float(H["/Vt"][()]) if "/Vt" in H else None
    return F, Q, U, V
