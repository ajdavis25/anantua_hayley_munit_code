#!/usr/bin/env python3
import sys, csv

# usage: python extract_munit.py CSV TIMESTEP MODEL SPIN POS FAMILY
# example: python extract_munit.py munits.csv 5000 RBETA -0.5 0 SANE

if len(sys.argv) != 7:
    print("usage: python extract_munit.py CSV TIMESTEP MODEL SPIN POS FAMILY", file=sys.stderr)
    sys.exit(2)

csv_path, timestep, model, spin, pos, family = sys.argv[1:]
timestep = str(int(float(timestep)))  # e.g. "5000"
model    = model.upper().strip()      # RBETA
spin     = spin.strip()               # "-0.5"
pos      = str(int(float(pos)))       # "0" or "1"
family   = family.upper().strip()     # SANE or MAD

with open(csv_path, "r", newline="") as f:
    # skip the first two non-data lines exactly as you said
    f.readline()
    f.readline()
    rdr = csv.DictReader(f)
    for row in rdr:
        # column names are taken exactly from your third line header
        if row["Timestep"].strip() != timestep:      continue
        if row["Model"].strip().upper() != model:    continue
        if row["Spin"].strip() != spin:              continue
        if row["positron frac"].strip() != pos:      continue
        if row["MAD/SANE"].strip().upper() != family:continue
        mu = row["MunitUsed"].strip()
        print(mu)
        sys.exit(0)

print("NOT_FOUND", file=sys.stderr)
sys.exit(1)
