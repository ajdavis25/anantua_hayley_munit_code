import csv
from dataclasses import dataclass
from typing import Dict, Tuple
from src.config import CSV_MUNIT

@dataclass
class MunitParams:
    timestep: int
    mad_sane: str
    model: str
    spin: float
    positron_frac: int
    munit_offset: float
    munit_slope: float
    munit: float
    munit_used: float
    spectrum_flag: str

    @property
    def A0(self): return self.munit_offset
    @property
    def S0(self): return self.munit_slope
    @property
    def M_base(self): return self.munit

def load_munit_params(csv_path=CSV_MUNIT) -> Dict[Tuple[int, str, str, float], MunitParams]:
    params_dict = {}
    with open(csv_path, newline='') as csvfile:
        for _ in range(2): next(csvfile)
        reader = csv.DictReader(csvfile)
        reader.fieldnames = [h.strip() for h in reader.fieldnames]
        for row in reader:
            row = {k.strip(): (v.strip() if v else "") for k, v in row.items()}
            if not row["Timestep"]: continue
            try:
                key = (int(row["Timestep"]), row["MAD/SANE"].upper(), row["Model"].upper(), float(row["Spin"].replace("+", "").replace("(+)", "")))
                if int(row["positron frac"]) != 0: continue
                params_dict[key] = MunitParams(
                    timestep=int(row["Timestep"]),
                    mad_sane=row["MAD/SANE"].upper(),
                    model=row["Model"].upper(),
                    spin=float(row["Spin"]),
                    positron_frac=int(row["positron frac"]),
                    munit_offset=float(row["MunitOffset"]),
                    munit_slope=float(row["MunitSlope"]),
                    munit=float(row["Munit"]),
                    munit_used=float(row["MunitUsed"]),
                    spectrum_flag=row["M87 GRMONTY spectrum?"]
                )
            except Exception as e:
                print(f"Skipping row: {e}")
    return params_dict
