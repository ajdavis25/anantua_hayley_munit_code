# src/config.py
#!/usr/bin/env python3
from pathlib import Path
from dataclasses import dataclass


@dataclass(frozen=True)
class SimDefaults:
    Rhigh: int = 20
    sigma_cut: float = 2.0
    thetacam: float = 163.0
    rmax_geo: int = 50
    freq_Hz: float = 228e9
    fov: float = 160.0
    npixel: int = 320
    counterjet: int = 0
    T0: float = 0.5
    T1: float = 0.5


defaults = SimDefaults()


IPATH: str = "/work/vmo703/aricarte/run_ipole.sh"
DUMPS: str = "/work/vmo703/grmhd_dump_samples/"
OUTPUT: str = "/work/vmo703/ipole_outputs/"
CSV_MUNIT: str = "/work/vmo703/data/munits_table.csv"
CSV_SUMMARY: Path = Path("/work/vmo703/data/summary_table.csv")


MODEL_SETTINGS = {
    "RBETA":        {"electronModel": 2, "sigma_transition": 2.0},
    "RBETAWJET":    {"electronModel": 2, "sigma_transition": 1.0},
    "CRITBETA":     {"electronModel": 4, "sigma_transition": 2.0},
    "CRITBETAWJET": {"electronModel": 4, "sigma_transition": 1.0},
}
