import csv
from pathlib import Path


results = {
    4000: {
        "MAD": {
            "RBETA": {
                -0.5: {
                    0.0: {"MunitOffset": 1.2e25, "MunitSlope": 2.5, "Munit_base": 7.5e24, "Munit_used": 3.4e25,
                          "flux_I": 0.51, "P/I": 3.9e-3, "V/I": -5.5e-3},
                    1.0: {"MunitOffset": 1.0e25, "MunitSlope": 2.5, "Munit_base": 7.5e24, "Munit_used": 2.0e25,
                          "flux_I": 0.48, "P/I": 1.6e-3, "V/I": -2.8e-3}
                }
            }
        }
    },
}


csv_path = Path("/work/vmo703/data/summary_table.csv")


# convert nested dict to flat rows for CSV
rows = []
for timestep, mad_dict in results.items():
    for mad_sane, model_dict in mad_dict.items():
        for model_name, spin_dict in model_dict.items():
            for spin, frac_dict in spin_dict.items():
                for positron_frac, data in frac_dict.items():
                    row = {
                        "Timestep": timestep,
                        "MAD/SANE": mad_sane,
                        "Model": model_name,
                        "Spin": spin,
                        "positron ratio": positron_frac,
                        "MunitOffset": data["MunitOffset"],
                        "MunitSlope": data["MunitSlope"],
                        "Munit_base": data["Munit_base"],
                        "Munit_used": data["Munit_used"],
                        "Flux I": data["flux_I"],
                        "P/I": data["P/I"],
                        "V/I": data["V/I"],
                    }
                    rows.append(row)


# write to CSV
fieldnames = list(rows[0].keys())
csv_path.parent.mkdir(parents=True, exist_ok=True)
with open(csv_path, "w", newline="") as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    writer.writerows(rows)


print(f"✅ Written summary to {csv_path}")
