#!/usr/bin/env python
"""Generate ipole pars for the Phase-4 movie frames (gated ipole).

24 branches (6 production families x 2 spins x pos0/1) x 11 dumps
(4000..6000 step 200) = 264 frames, consistent with the Phase-4 SPECTRA:
  - gated ipole: constant_beta_paper_literal 1 on wJET families,
    thetae_floor 3e-2 on Crit-beta families;
  - MBH 6.5e9, dsource 16.8 Mpc (the grmonty production values -- the older
    image set used 6.2e9/16.9; changed deliberately for spectra consistency);
  - thetacam 163 deg (group imaging convention), 230 GHz, 320x320, fov 160 uas,
    sigma_cut 10 (ipole has no grmonty-style jet Thetae override; sigma>10
    carries ~1.5% of MAD emission);
  - M_unit per dump: ln-linear interpolation across the three tuned anchors
    (4000/5000/6000) read from _qa/p4_qa_run_2026-09-15.csv.

Run: /work/vmo703/ipole_venv/bin/python generate_pars.py
"""

import csv
import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))
QA_CSV = "/work/vmo703/igrmonty_outputs/m87/_qa/p4_qa_run_2026-09-15.csv"
DUMP_DIR = "/work/vmo703/grmhd_dump_samples"
DUMPS = list(range(4000, 6001, 200))
ANCHORS = (4000, 5000, 6000)

# family -> (electronModel, extra par lines, rh label)
FAMILIES = {
    "MAD_RBETA": (2, ["trat_small 1", "trat_large 20",
                      "sigma_transition 1e20"], 20),
    "MAD_RBETAwJET": (2, ["trat_small 1", "trat_large 80",
                          "sigma_transition 2.0", "constant_beta_e0 0.1",
                          "constant_beta_e0_exponent 1.0",
                          "constant_beta_paper_literal 1"], 80),
    "MAD_CRITBETA": (4, ["trat_small 1", "trat_large 20",
                         "beta_crit 0.01", "beta_crit_coefficient 0.5",
                         "sigma_transition 1e20", "thetae_floor 3.0e-2"], 20),
    "MAD_CRITBETAwJET": (4, ["trat_small 1", "trat_large 20",
                             "beta_crit 0.01", "beta_crit_coefficient 0.5",
                             "sigma_transition 2.0", "constant_beta_e0 0.1",
                             "constant_beta_e0_exponent 1.0",
                             "constant_beta_paper_literal 1",
                             "thetae_floor 3.0e-2"], 20),
    "SANE_RBETAwJET": (2, ["trat_small 1", "trat_large 160",
                           "sigma_transition 2.0", "constant_beta_e0 0.1",
                           "constant_beta_e0_exponent 1.0",
                           "constant_beta_paper_literal 1"], 160),
    "SANE_CRITBETAwJET": (4, ["trat_small 1", "trat_large 20",
                              "beta_crit 1", "beta_crit_coefficient 0.5",
                              "sigma_transition 2.0", "constant_beta_e0 0.1",
                              "constant_beta_e0_exponent 1.0",
                              "constant_beta_paper_literal 1",
                              "thetae_floor 3.0e-2"], 20),
}

COMMON = """thetacam 163.0
phicam 0.0
freqcgs 230.0e9
MBH 6.5e9
dsource 16.8e6
nx 320
ny 320
fovx_dsource 160.0
fovy_dsource 160.0
counterjet 0
rmax_geo 50
emission_type 4
sigma_cut 10.0
"""


def qa_munits():
    """(family, spin, pos, dump) -> tuned M_unit from the campaign finals."""
    out = {}
    with open(QA_CSV) as fh:
        for r in csv.DictReader(fh):
            n = r["spectrum"].replace("spectrum_", "")
            # e.g. Ma+0.94_4000_RBETAwJET_rh80_pos0.h5
            state = "MAD" if n.startswith("Ma") else "SANE"
            spin = n[2:n.index("_")]
            rest = n[n.index("_") + 1:]
            dump = int(rest.split("_")[0])
            model = rest.split("_")[1]
            pos = "pos1" if "_pos1" in n else "pos0"
            fam = "%s_%s" % (state, model)
            out[(fam, spin, pos, dump)] = float(r["M_unit"])
    return out


def interp_munit(mu, fam, spin, pos, dump):
    a = {d: mu[(fam, spin, pos, d)] for d in ANCHORS}
    if dump in a:
        return a[dump]
    lo = max(d for d in ANCHORS if d < dump)
    hi = min(d for d in ANCHORS if d > dump)
    t = (dump - lo) / float(hi - lo)
    return math.exp((1 - t) * math.log(a[lo]) + t * math.log(a[hi]))


def main():
    mu = qa_munits()
    os.makedirs(os.path.join(HERE, "pars"), exist_ok=True)
    os.makedirs(os.path.join(HERE, "img"), exist_ok=True)
    branches = []
    for fam, (emodel, extra, rh) in sorted(FAMILIES.items()):
        state = fam.split("_")[0]
        prefix = "Ma" if state == "MAD" else "Sa"
        for spin in ("+0.94", "-0.5"):
            for pos in ("pos0", "pos1"):
                branch = "%s%s_%s_rh%d_%s" % (prefix, spin,
                                              fam.split("_", 1)[1], rh, pos)
                lst = []
                for dump in DUMPS:
                    m = interp_munit(mu, fam, spin, pos, dump)
                    par = os.path.join(HERE, "pars",
                                       "%s_%d.par" % (branch, dump))
                    img = os.path.join(HERE, "img",
                                       "%s_%d.h5" % (branch, dump))
                    with open(par, "w") as fh:
                        fh.write("dump %s/%s%s_%d.h5\n"
                                 % (DUMP_DIR, prefix, spin, dump))
                        fh.write("outfile %s\n" % img)
                        fh.write("M_unit %.8e\n" % m)
                        fh.write("positronRatio %d\n"
                                 % (1 if pos == "pos1" else 0))
                        fh.write("electronModel %d\n" % emodel)
                        for line in extra:
                            fh.write(line + "\n")
                        fh.write(COMMON)
                    lst.append(par)
                lf = os.path.join(HERE, "pars", branch + ".list")
                with open(lf, "w") as fh:
                    fh.write("\n".join(lst) + "\n")
                branches.append(branch)
    with open(os.path.join(HERE, "branches.txt"), "w") as fh:
        fh.write("\n".join(branches) + "\n")
    print("%d branches x %d dumps = %d pars"
          % (len(branches), len(DUMPS), len(branches) * len(DUMPS)))


if __name__ == "__main__":
    main()
