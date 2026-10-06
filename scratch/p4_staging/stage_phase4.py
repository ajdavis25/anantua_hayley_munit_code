#!/usr/bin/env python
"""Stage Phase-4 production: fix the tuner warm-start CSV + add MAD wJET rows.

1. Replaces the stale April M_units in data/final_grmonty_paper.csv
   (columns MunitUsed_pos0/pos1 -- the values auto_munit_bracket.py actually
   warm-starts from) with the July-2026 tuned values read from each archived
   spectrum's /params/M_unit in igrmonty_outputs/m87/pre_bugfix_2026-07/.
2. Records each branch's fitted bias (/params/bias) in new informational
   columns bias_jul_pos0/pos1 (names deliberately do NOT match the tuner's
   ^MunitUsed_pos\\d$ discovery regex).
3. Flags branches whose July campaign never converged (from the committed
   munits_tuning_history.csv) in the notes column.
4. Appends 6 new MAD RBETAwJET rh80 rows (spins +0.94/-0.5 x dumps
   4000/5000/6000) -- the go-forward paper configs, which the edited
   auto_munit_bracket.py will run with constant_beta_paper_literal 1.
   Their warm-start M_unit is prescaled from the completed paper-literal
   A/B spectrum's 230 GHz flux (tuner power-law, p=2); pos1 uses the median
   July MAD pos1/pos0 ratio.

The original CSV is backed up beside this script before writing.
Run: /work/vmo703/ipole_venv/bin/python stage_phase4.py
"""

import csv
import importlib.util
import shutil
from pathlib import Path

import h5py
import numpy as np

CSV_PATH = Path("/work/vmo703/data/final_grmonty_paper.csv")
HISTORY_CSV = Path("/work/vmo703/data/munits_tuning_history.csv")
ARCHIVE = Path("/work/vmo703/igrmonty_outputs/m87/pre_bugfix_2026-07")
DUMP_DIR = Path("/work/vmo703/grmhd_dump_samples")
STAGING = Path("/work/vmo703/scratch/p4_staging")
PAPERLIT_SPEC = Path("/work/vmo703/scratch/pb_fix_ab/out/spec_mad_wjet_paperlit.h5")
PAPERLIT_MUNIT = 1.3488e25  # M_unit the A/B pair ran at (P31 legacy-tuned)
F_TARGET_JY = 0.5
P_FACTOR = 2.0  # tuner default --p-factor
NEW_RHIGH = 80  # MAD wJET production Rhigh (P3 jetcheck lineage)

VCP_SPEC = importlib.util.spec_from_file_location(
    "vcp", "/work/vmo703/igrmonty/tools/viewing_cone_postprocess.py")
vcp = importlib.util.module_from_spec(VCP_SPEC)
VCP_SPEC.loader.exec_module(vcp)


def spin_tag(spin_str):
    v = float(spin_str)
    return ("+" if v >= 0 else "") + spin_str.lstrip("+")


def archived_h5(state, spin, dump_idx, model, rhigh, pos):
    prefix = {"SANE": "S", "MAD": "M"}[state]
    pat = "spectrum_{}a{}_{}_{}_rh{}*_pos{}.h5".format(
        prefix, spin_tag(spin), dump_idx, model, rhigh, pos)
    hits = sorted(ARCHIVE.glob(pat))
    if len(hits) != 1:
        raise RuntimeError("expected 1 archive match for %s, got %d: %s"
                           % (pat, len(hits), [h.name for h in hits]))
    return hits[0]


def read_params(h5path):
    with h5py.File(str(h5path), "r") as f:
        return float(f["/params/M_unit"][()]), float(f["/params/bias"][()])


def july_converged_map():
    """(state, model, spin, dump, pos) -> True if any current-era history row converged.

    The history file mixes schemas: the 2025 header row has 19 columns, while
    rows appended by the current tuner (the July campaign and later) have 27
    (positron_ratio + 7 wJET columns inserted before 'converged'). DictReader
    against the old header misparses them, so parse positionally by width.
    """
    conv = {}
    with HISTORY_CSV.open() as fh:
        for i, row in enumerate(csv.reader(fh)):
            if i == 0 or not row:
                continue
            if len(row) >= 27:          # current-era schema
                state, model, spin, dump_idx, pos = row[2], row[3], row[4], row[5], row[6]
                converged = row[22]
            elif len(row) == 19:        # 2025-era schema (pre-campaign)
                continue
            else:
                continue
            key = (state.strip().upper(), model.strip(),
                   spin.strip(), dump_idx.strip(), pos.strip())
            conv[key] = conv.get(key, False) or converged.strip() in ("1", "True", "true")
    return conv


def paperlit_f230():
    spec = vcp.load_spectrum(str(PAPERLIT_SPEC))
    mask = np.ones(spec["dOmega_sr"].shape[0], dtype=bool)
    modes = vcp.compute_spectra_modes(
        spec["nu_hz"], spec["nuLnu_theta_freq_cgs"], spec["dOmega_sr"],
        mask, vcp.D_MPC if hasattr(vcp, "D_MPC") else 16.8)
    fnu = np.asarray(modes["fnu_4pi_jy"], dtype=float)
    nu = spec["nu_hz"]
    good = fnu > 0
    return float(np.exp(np.interp(np.log(230.0e9), np.log(nu[good]), np.log(fnu[good]))))


def main():
    STAGING.mkdir(parents=True, exist_ok=True)
    backup = STAGING / (CSV_PATH.name + ".bak_apr30")
    if not backup.exists():
        shutil.copy2(str(CSV_PATH), str(backup))
        print("[backup] %s" % backup)

    with CSV_PATH.open() as fh:
        reader = csv.DictReader(fh)
        fieldnames = list(reader.fieldnames)
        rows = list(reader)
    print("[read] %d rows from %s" % (len(rows), CSV_PATH))

    for col in ("bias_jul_pos0", "bias_jul_pos1"):
        if col not in fieldnames:
            fieldnames.append(col)

    conv = july_converged_map()
    mad_ratios = []

    print("\n%-42s %-4s %12s %12s %7s %6s %s" % (
        "branch", "pos", "M_apr", "M_jul", "ratio", "bias", "jul_conv"))
    for row in rows:
        state = row["state"].strip().upper()
        model = row["model"].strip()
        spin = row["spin"].strip()
        dump_idx = row["dump_index"].strip()
        rhigh = row["Rhigh"].strip()
        tag = "%s %s a%s d%s rh%s" % (state, model, spin, dump_idx, rhigh)
        unconv = []
        for pos in ("0", "1"):
            h5p = archived_h5(state, spin, dump_idx, model, rhigh, pos)
            m_jul, bias_jul = read_params(h5p)
            key_old = "MunitUsed_pos%s" % pos
            m_apr = float(row[key_old])
            key = (state, model, spin, dump_idx, pos)
            is_conv = conv.get(key, False)
            if not is_conv:
                unconv.append(pos)
            print("%-42s %-4s %12.4e %12.4e %7.3f %6.3g %s" % (
                tag, pos, m_apr, m_jul, m_jul / m_apr, bias_jul,
                "yes" if is_conv else "NO"))
            row[key_old] = "%.8e" % m_jul
            row["bias_jul_pos%s" % pos] = "%.6g" % bias_jul
        if state == "MAD":
            mad_ratios.append(float(row["MunitUsed_pos1"]) / float(row["MunitUsed_pos0"]))
        note_add = "; munit=jul2026 pre_bugfix h5"
        if unconv:
            note_add += " (jul_unconverged pos%s: best-so-far)" % ",".join(unconv)
        row["notes"] = (row.get("notes") or "").rstrip() + note_add

    # ---- new MAD wJET production rows -------------------------------------
    f230 = paperlit_f230()
    m_start = PAPERLIT_MUNIT * (F_TARGET_JY / f230) ** (1.0 / P_FACTOR)
    pos1_ratio = float(np.median(mad_ratios))
    print("\n[paperlit] F230(4pi) = %.4g Jy at M_unit %.4e" % (f230, PAPERLIT_MUNIT))
    print("[warmstart] new-row MunitUsed_pos0 = %.4e (p=%g rescale to %.2g Jy)"
          % (m_start, P_FACTOR, F_TARGET_JY))
    print("[warmstart] pos1/pos0 median over July MADs = %.4f" % pos1_ratio)

    next_id = max(int(r["row_id"]) for r in rows) + 1
    for spin in ("0.94", "-0.5"):
        for dump_idx in ("4000", "5000", "6000"):
            dump_file = DUMP_DIR / ("Ma%s_%s.h5" % (spin_tag(spin), dump_idx))
            if not dump_file.exists():
                raise RuntimeError("missing GRMHD dump: %s" % dump_file)
            new = {k: "" for k in fieldnames}
            new.update({
                "row_id": str(next_id),
                "dump_index": dump_idx,
                "timestep": str(int(dump_idx) * 5),
                "state": "MAD",
                "model": "RBETAwJET",
                "spin": spin,
                "Rhigh": str(NEW_RHIGH),
                "Munit": "%.8e" % 7.49e24,
                "MunitUsed_pos0": "%.8e" % m_start,
                "MunitUsed_pos1": "%.8e" % (m_start * pos1_ratio),
                "converged": "False",
                "notes": ("staged 2026-09-09: new MAD wJET production branch, "
                          "paper-literal P_B (constant_beta_paper_literal=1 via "
                          "WJET_DEFAULTS); warm start = P31 M_unit rescaled by "
                          "paperlit F230=%.3gJy; pos1 via median MAD ratio" % f230),
            })
            rows.append(new)
            next_id += 1
    print("[append] 6 MAD RBETAwJET rh%d rows (row_id 24-29)" % NEW_RHIGH)

    with CSV_PATH.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print("[write] %s: %d rows x pos0/pos1 = %d production tasks"
          % (CSV_PATH, len(rows), 2 * len(rows)))
    print("[launch] sbatch --array=0-%d%%5 run_auto_munit.slurm" % (2 * len(rows) - 1))


if __name__ == "__main__":
    main()
