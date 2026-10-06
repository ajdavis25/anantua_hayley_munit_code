#!/usr/bin/env python
"""Render PNG frames from the Phase-4 ipole images in the HOUSE STYLE.

Reuses the user's own two-panel plotter (sgrA/create_images.py:
plotPositronTestFrame -- Stokes I in afmhot with polarization ticks and P/I,
plus a Stokes V panel in seismic, colorbars, 8x4 in) verbatim, with the same
arguments the production driver uses (cpMax=0.1, fractionalCircular=False).
The function definition is exec'd from the file WITHOUT the module-level
driver loops at the bottom (no __main__ guard there).

PNG names follow the movie script's primary convention:
  output_{Ma|Sa}{spin}_{dump}_model_{MODEL}_Rhigh_{rh}_{0.000|1.000}.png

Run after the render array: /work/vmo703/ipole_venv/bin/python plot_frames.py
"""

import os
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
PNG = os.path.join(HERE, "png")
CREATE_IMAGES = "/work/vmo703/sgrA/create_images.py"

# import plotPositronTestFrame without executing the sgrA driver loops
src = open(CREATE_IMAGES).read().split("# handle .h5 files")[0]
ns = {}
exec(compile(src, CREATE_IMAGES, "exec"), ns)
plot_fn = ns["plotPositronTestFrame"]


# our image files end in the dump index (..._pos0_5000.h5), which defeats the
# original last-token fallback in _infer_positron_ratio (it printed the dump
# time as the pair ratio). Override it in the function's namespace with a
# version keyed on the _pos<N>_ segment anywhere in the name.
def _infer_pos_fixed(imageFile):
    m = re.search(r"_pos([01])(?:[_.]|$)", os.path.basename(imageFile))
    return float(m.group(1)) if m else float("nan")


ns["_infer_positron_ratio"] = _infer_pos_fixed

BR = re.compile(r"^(?P<sa>Ma|Sa)(?P<spin>[+-][\d.]+)_(?P<model>\w+?)_rh"
                r"(?P<rh>\d+)_(?P<pos>pos[01])$")


def main():
    os.makedirs(PNG, exist_ok=True)
    with open(os.path.join(HERE, "branches.txt")) as fh:
        branches = [l.strip() for l in fh if l.strip()]
    nmade = nfail = 0
    for branch in branches:
        m = BR.match(branch)
        for dump in range(4000, 6001, 200):
            h5p = os.path.join(HERE, "img", "%s_%d.h5" % (branch, dump))
            if not os.path.exists(h5p):
                continue
            posv = "1.000" if m.group("pos") == "pos1" else "0.000"
            name = ("output_%s%s_%d_model_%s_Rhigh_%s_%s.png"
                    % (m.group("sa"), m.group("spin"), dump,
                       m.group("model"), m.group("rh"), posv))
            try:
                plot_fn(h5p, cpMax=0.1, fractionalCircular=False,
                        output=os.path.join(PNG, name))
                nmade += 1
            except Exception as e:
                print("[fail] %s: %s" % (name, e))
                nfail += 1
            plt.close("all")
        print("[branch] %s done" % branch)
    print("%d PNGs (%d failed) -> %s" % (nmade, nfail, PNG))


if __name__ == "__main__":
    main()
