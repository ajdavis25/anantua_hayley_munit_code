#!/usr/bin/env python
"""Production pass/fail table for the Phase-4 campaign (all 72 branches).

One row per expected branch (36 CSV rows x pos0/pos1), whether landed or not:
  F230 gate   -- tuner-convention 230 GHz flux within 10% of 0.5 Jy
  LX 4pi gate -- L_X(2-10 keV, sky-averaged) <= 4.4e40 erg/s (2017 core limit)
  LX cone17   -- same limit in the 17-degree observer frame (informational)
  provenance  -- Ns/flag/status as production demands
  OVERALL     -- PASS iff F230 + LX(4pi) + provenance all pass; else FAIL;
                 'pending' until the branch lands.

Outputs: p4_passfail.csv and p4_passfail.md next to this script.
The QA watchdog reruns this after every landing, so the table is live.

Run: /work/vmo703/ipole_venv/bin/python p4_passfail.py
"""

import csv
import os
from datetime import datetime

HERE = os.path.dirname(os.path.abspath(__file__))
QA_CSV = os.path.join(HERE, "p4_qa_run_2026-09-15.csv")
BRANCH_CSV = "/work/vmo703/data/final_grmonty_paper.csv"
OUT_CSV = os.path.join(HERE, "p4_passfail.csv")
OUT_MD = os.path.join(HERE, "p4_passfail.md")
LX_GATE = 4.4e40


def spin_tag(spin):
    v = float(spin)
    return ("+" if v >= 0 else "") + spin.lstrip("+")


def expected_name(row, pos):
    prefix = {"SANE": "S", "MAD": "M"}[row["state"].strip().upper()]
    model = row["model"].strip()
    tail = "_rh%s" % row["Rhigh"].strip()
    if model.upper().startswith("CRITBETA"):
        tail += "_bc%s_f%s" % (row["beta_crit"].strip(), row["f"].strip())
    return "spectrum_%sa%s_%s_%s%s_pos%s.h5" % (
        prefix, spin_tag(row["spin"].strip()), row["dump_index"].strip(),
        model, tail, pos)


def main():
    with open(BRANCH_CSV) as fh:
        branches = list(csv.DictReader(fh))
    qa = {}
    if os.path.exists(QA_CSV):
        with open(QA_CSV) as fh:
            for r in csv.DictReader(fh):
                qa[r["spectrum"]] = r

    rows = []
    for b in branches:
        family = "%s %s" % (b["state"].strip(), b["model"].strip())
        for pos in ("0", "1"):
            name = expected_name(b, pos)
            r = qa.get(name)
            if r is None:
                rows.append(dict(branch=name.replace("spectrum_", "").replace(".h5", ""),
                                 family=family, pos=pos, F230="-", LX_4pi="-",
                                 LX_cone17="-", provenance="-", OVERALL="pending",
                                 LX_over_gate=""))
                continue
            lx_cone = float(r["L_X_cone17"]) if r.get("L_X_cone17") else float("nan")
            cone_gate = "PASS" if lx_cone <= LX_GATE else "FAIL"
            overall = ("PASS" if (r["gate_F230"] == "PASS" and
                                  r["gate_LX_4pi"] == "PASS" and
                                  r["gate_provenance"] == "PASS") else "FAIL")
            rows.append(dict(branch=name.replace("spectrum_", "").replace(".h5", ""),
                             family=family, pos=pos,
                             F230=r["gate_F230"], LX_4pi=r["gate_LX_4pi"],
                             LX_cone17=cone_gate, provenance=r["gate_provenance"],
                             OVERALL=overall,
                             LX_over_gate="%.2f" % (float(r["L_X_4pi"]) / LX_GATE)))

    cols = ["branch", "family", "pos", "F230", "LX_4pi", "LX_cone17",
            "provenance", "LX_over_gate", "OVERALL"]
    with open(OUT_CSV, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)

    landed = [r for r in rows if r["OVERALL"] != "pending"]
    npass = sum(r["OVERALL"] == "PASS" for r in landed)
    nfail = len(landed) - npass
    with open(OUT_MD, "w") as fh:
        fh.write("# Phase-4 production pass/fail table\n\n")
        fh.write("Updated %s UTC by p4_passfail.py (auto-refreshed by the QA "
                 "watchdog on every landing).\n\n" %
                 datetime.utcnow().strftime("%Y-%m-%d %H:%M"))
        fh.write("**%d/72 landed: %d PASS, %d FAIL, %d pending.** "
                 % (len(landed), npass, nfail, len(rows) - len(landed)))
        fh.write("Gates: F230 = 0.5 Jy anchor (tuner convention, 10%% band); "
                 "LX = 2-10 keV vs the 2017 core limit 4.4e40 erg/s at 4pi "
                 "(cone17 = same limit at the 17-deg observer frame, "
                 "informational); LX/gate column is the 4pi ratio.\n\n")
        fh.write("| branch | pos | F230 | LX 4pi | LX cone17 | prov | LX/gate | OVERALL |\n")
        fh.write("|---|---|---|---|---|---|---|---|\n")
        for r in rows:
            fh.write("| %s | %s | %s | %s | %s | %s | %s | %s |\n" % (
                r["branch"], r["pos"], r["F230"], r["LX_4pi"], r["LX_cone17"],
                r["provenance"], r["LX_over_gate"], r["OVERALL"]))
    print("[passfail] %d/72 landed: %d PASS, %d FAIL, %d pending" %
          (len(landed), npass, nfail, len(rows) - len(landed)))
    print("[passfail] %s" % OUT_MD)


if __name__ == "__main__":
    main()
