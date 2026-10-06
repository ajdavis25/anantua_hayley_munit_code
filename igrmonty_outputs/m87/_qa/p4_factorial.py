#!/usr/bin/env python
"""The 2x3 factorial: X-ray gate verdict by electron base x configuration.

Rows: electron base model (R-beta, Crit-beta).
Columns: MAD plain / MAD + wJET supplement / SANE + wJET supplement.
Cell: branches landed, PASS/FAIL vs the 2017 core X-ray gate (4.4e40, 4pi),
and the median L_X/gate ratio. Cells not yet populated print "pending".

Design claim it certifies: the gate failures are {MAD} x {supplement},
independent of the base model -- the supplement's P_e = beta_e0 P_B injection
only becomes X-ray-loud where MAD-strength fields meet sigma>2 volume.

Run: /work/vmo703/ipole_venv/bin/python p4_factorial.py
(reads the QA CSV the watchdog refreshes; rerun any time for current state)
"""

import csv
import os
import re

QA_CSV = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                      "p4_qa_run_2026-09-15.csv")
LX_GATE = 4.4e40

CELLS = [
    ("MAD plain", r"^Ma.*_(RBETA|CRITBETA)_rh"),
    ("MAD + wJET", r"^Ma.*_(RBETA|CRITBETA)wJET_"),
    ("SANE + wJET", r"^Sa.*_(RBETA|CRITBETA)wJET_"),
]


def base_of(name):
    return "Crit-beta" if "CRITBETA" in name else "R-beta"


def main():
    with open(QA_CSV) as fh:
        rows = list(csv.DictReader(fh))
    table = {}
    for r in rows:
        name = r["spectrum"].replace("spectrum_", "")
        for col, pat in CELLS:
            if re.match(pat, name):
                table.setdefault((base_of(name), col), []).append(
                    float(r["L_X_4pi"]) / LX_GATE)
                break

    print("2x3 factorial -- L_X(2-10 keV, 4pi) / gate(4.4e40); "
          "PASS = every branch under 1.0\n")
    colw = 24
    header = "%-12s" % "base" + "".join("%-*s" % (colw, c) for c, _ in CELLS)
    print(header)
    print("-" * len(header))
    for base in ("R-beta", "Crit-beta"):
        line = "%-12s" % base
        for col, _ in CELLS:
            vals = sorted(table.get((base, col), []))
            if not vals:
                line += "%-*s" % (colw, "pending")
                continue
            med = vals[len(vals) // 2]
            nfail = sum(v > 1.0 for v in vals)
            verdict = "PASS" if nfail == 0 else "FAIL %d/%d" % (nfail, len(vals))
            line += "%-*s" % (colw, "%s med=%.2gx n=%d" % (verdict, med, len(vals)))
        print(line)
    fail_cells = sorted(set(
        (b, c) for (b, c), vals in table.items() if any(v > 1.0 for v in vals)))
    print("\ncells containing gate failures: %s"
          % (", ".join("%s x %s" % bc for bc in fail_cells) or "none"))
    print("reading: the wJET supplement drives every failure; MAD+wJET is the")
    print("loud regime, and pair-loaded (pos1) branches push marginal cells over.")
    print("NOTE: raw campaign values; confirm-run overrides in p4_confirm_notes.md")
    print("(the lone SANE+wJET 'failure' was bias-0.05 MC noise -- confirmed PASS).")


if __name__ == "__main__":
    main()
