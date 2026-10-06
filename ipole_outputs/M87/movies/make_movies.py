#!/usr/bin/env python3
"""Build per-branch movies from M87 ipole PNGs (dump sequence as time axis).

Adapted 2026-10-03 from aricarte PositronIPOLEScripts/folderToMovie.py:
  - regex updated to the current naming:
      output_{Ma|Sa}{spin}_{dump}_model_{MODEL}_Rhigh_{N}_{pos}.png
    (the original pattern predates the _model_/_Rhigh_ segments and skips
    every current file);
  - grouping key now includes spin and Rhigh so different branches never
    merge into one movie;
  - frames ordered by dump index (4000 -> 5000 -> 6000), framerate 2 fps
    (three dumps = the whole time axis for this grid);
  - mpeg4 .mp4 output (this ffmpeg build has no libx264);
  - output tree: movies/{MAD|SANE}/{MODEL}/a{spin}_rh{N}_pos{P}.mp4

Sources: M87/images/ and M87/{Ma,Sa}*/ per-dump dirs; scrap/ and betacrit/
are excluded. Groups with fewer than 2 frames are listed, not rendered.

Run: python3 make_movies.py
"""

import os
import re
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

ROOT = Path("/work/vmo703/ipole_outputs/M87")
OUT = ROOT / "movies"
FFMPEG = "/usr/bin/ffmpeg"  # system build: has libx264 (browser-playable)
FRAMERATE = 2

PAT = re.compile(
    r"output_(?P<sa>Ma|Sa)(?P<spin>[+-][\d.]+)_(?P<dump>\d+)_model_"
    r"(?P<model>RBETAwJET|RBETA|CRITBETAwJET|CRITBETA)_Rhigh_(?P<rh>\d+)_"
    r"(?P<pos>[\d.]+)\.png", re.IGNORECASE)

# second naming convention (per-dump production dirs): Rhigh lives in the
# PATH (e.g. Sa+0.94_5000/RBETA/Rhigh20/pos1/SANE_spin+0.94_t5000_RBETA_pos1.png)
PAT2 = re.compile(
    r"(?P<state>MAD|SANE)_spin(?P<spin>[+-][\d.]+)_t(?P<dump>\d+)_"
    r"(?P<model>RBETAWJET|RBETA|CRITBETAWJET|CRITBETA)_pos(?P<pos>\d)\.png",
    re.IGNORECASE)
RH_IN_PATH = re.compile(r"Rhigh(\d+)", re.IGNORECASE)

# third convention: Sa+0.94_6000_RBETAwJET_Rhigh80_1.00.png
PAT3 = re.compile(
    r"(?P<sa>Ma|Sa)(?P<spin>[+-][\d.]+)_(?P<dump>\d+)_"
    r"(?P<model>RBETAwJET|RBETA|CRITBETAwJET|CRITBETA)_Rhigh(?P<rh>\d+)_"
    r"(?P<pos>[\d.]+)\.png", re.IGNORECASE)


def main():
    pngs = []
    for p in ROOT.rglob("*.png"):
        rel = p.relative_to(ROOT).parts
        if rel[0] in ("scrap", "betacrit", "movies"):
            continue
        pngs.append(p)

    groups = defaultdict(dict)  # key -> {dump: path}
    skipped = 0
    # frames_p4 (gated-ipole, production-consistent) wins any collision with
    # the older image sets for the same (branch, dump)
    for p in sorted(pngs, key=lambda q: ("frames_p4" not in str(q), str(q))):
        m = PAT.search(p.name)
        if m:
            state = "MAD" if m.group("sa") == "Ma" else "SANE"
            pos = "pos1" if float(m.group("pos")) >= 0.5 else "pos0"
            rh = int(m.group("rh"))
        else:
            m = PAT2.search(p.name)
            if m:
                state = m.group("state").upper()
                pos = "pos" + m.group("pos")
                mrh = RH_IN_PATH.search(str(p))
                rh = int(mrh.group(1)) if mrh else 0
            else:
                m = PAT3.search(p.name)
                if not m:
                    skipped += 1
                    continue
                state = "MAD" if m.group("sa") == "Ma" else "SANE"
                pos = "pos1" if float(m.group("pos")) >= 0.5 else "pos0"
                rh = int(m.group("rh"))
        key = (state, m.group("model").upper(), m.group("spin"), rh, pos)
        # first occurrence wins per (key, dump): duplicates across source
        # dirs carry the same content for the same branch+dump
        groups[key].setdefault(int(m.group("dump")), p)

    made, thin = 0, []
    for key in sorted(groups):
        state, model, spin, rh, pos = key
        frames = [groups[key][d] for d in sorted(groups[key])]
        name = "a%s_rh%d_%s" % (spin, rh, pos)
        if len(frames) < 2:
            thin.append(("%s/%s/%s" % (state, model, name), len(frames)))
            continue
        outdir = OUT / state / model
        outdir.mkdir(parents=True, exist_ok=True)
        outfile = outdir / (name + ".mp4")
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as tf:
            for f in frames:
                tf.write("file '%s'\n" % f)
                tf.write("duration %.3f\n" % (1.0 / FRAMERATE))
            tf.write("file '%s'\n" % frames[-1])  # concat quirk: repeat last
            listfile = tf.name
        cmd = [FFMPEG, "-y", "-f", "concat", "-safe", "0", "-i", listfile,
               "-vf", "scale=1280:-2", "-vcodec", "libx264",
               "-pix_fmt", "yuv420p", "-crf", "20", "-movflags", "+faststart",
               "-r", "8", str(outfile)]
        res = subprocess.run(cmd, stdout=subprocess.DEVNULL,
                             stderr=subprocess.PIPE)
        os.unlink(listfile)
        if res.returncode == 0 and outfile.exists():
            made += 1
            print("[ok] %s (%d frames)" % (outfile.relative_to(ROOT), len(frames)))
        else:
            print("[FAIL] %s\n%s" % (outfile, res.stderr.decode()[-400:]),
                  file=sys.stderr)
    print("\n%d movies made; %d unmatched pngs skipped" % (made, skipped))
    if thin:
        print("groups with <2 frames (no movie): %d" % len(thin))
        for g, n in thin[:12]:
            print("  %s (%d frame)" % (g, n))


if __name__ == "__main__":
    main()
