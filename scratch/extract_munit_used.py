import sys, csv
csv_path, dump, model, spin, pos, family = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4], int(sys.argv[5]), sys.argv[6]
with open(csv_path, newline='') as f:
    r = csv.reader(f)
    rows = list(r)[3:]  # data starts at line 4
    for row in rows:
        if not row or len(row)<9: 
            continue
        t, fam, mod, sp, p = int(row[0]), row[1].strip(), row[2].strip(), row[3].strip(), int(row[4])
        if t==dump and fam==family and mod==model and sp==spin and p==pos:
            off, slope, base = float(row[5]), float(row[6]), float(row[7])
            used = off + slope*base/(1+2*pos)
            print(f"{used:.6E}")
            sys.exit(0)
print("NA")
