#!/usr/bin/env python
"""Numpy prediction of the paper-literal (B^2/8pi) wJET Thetae field.

At constant_beta_e0_exponent = 1 the paper-literal form equals the legacy form
with e0 -> e0/(12 pi) exactly, so the validated numpy port predicts the fixed
field without new physics code. Saves full-grid assigned Thetae for both
variants (later compared zone-by-zone against the C binary's debug dumps) and
renders the same slice figure for the paper-literal case.

Run: /work/vmo703/ipole_venv/bin/python predict_fields.py
"""

import importlib.util

import numpy as np

SPEC = importlib.util.spec_from_file_location(
    "tsc", "/work/vmo703/igrmonty_outputs/m87/_qa/thetae_slice_check.py")
tsc = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(tsc)
stb = tsc.stb

OUT = "/work/vmo703/scratch/pb_fix_ab/out"

cfg = stb.PRESETS[0]  # MAD Ma+0.94_4000 RBETAwJET rh80, exponent = 1.0
assert cfg["constant_beta_e0_exponent"] == 1.0

tsc.run_preset(cfg, save_npy=f"{OUT}/assigned_numpy_legacy.npy")

cfg_pl = dict(cfg, tag=cfg["tag"] + "_PAPERLIT",
              constant_beta_e0=cfg["constant_beta_e0"] / (12.0 * np.pi))
tsc.run_preset(cfg_pl, save_npy=f"{OUT}/assigned_numpy_paperlit.npy")
