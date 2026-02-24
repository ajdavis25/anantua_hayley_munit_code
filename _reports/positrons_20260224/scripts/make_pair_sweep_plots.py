#!/usr/bin/env python3

import csv
from pathlib import Path

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPORT_DIR = Path("/work/vmo703/_reports/positrons_20260224")
CSV_PATH = REPORT_DIR / "tables" / "pair_sweep_20260224.csv"
FIG_DIR = REPORT_DIR / "figures"


def read_rows(csv_path: Path):
    rows = []
    with csv_path.open("r", newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            row["positron_ratio"] = float(row["positron_ratio"])
            for key in (
                "nu_used_hz",
                "flux_jy",
                "tau_scatt_weighted",
                "tau_abs_weighted",
                "tau_scatt_mean",
                "tau_abs_mean",
                "flux_scale_vs_pos0",
                "tau_scatt_weighted_scale_vs_pos0",
                "tau_abs_weighted_scale_vs_pos0",
                "tau_scatt_mean_scale_vs_pos0",
                "tau_abs_mean_scale_vs_pos0",
            ):
                row[key] = float(row[key])
            rows.append(row)
    rows.sort(key=lambda r: r["positron_ratio"])
    return rows


def load_isotropic_sed(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        numin = float(f["/params/NUMIN"][()])
        numax = float(f["/params/NUMAX"][()])
        n_ebins = int(f["/params/N_EBINS"][()])
        nu = np.exp(np.linspace(np.log(numin), np.log(numax), n_ebins))

        nu_lnu = f["/output/nuLnu"][:]  # [typebin, ebin, thbin]
        d_omega = f["/output/dOmega"][:]  # [thbin]
        ang_weights = d_omega / (4.0 * np.pi)

        # Angle-integrated spectrum used consistently with model luminosity integration.
        sed_by_type = np.tensordot(nu_lnu, ang_weights, axes=([2], [0]))
        sed_total = sed_by_type.sum(axis=0)
    return nu, sed_total


def plot_sed_overlay(rows):
    fig, ax = plt.subplots(figsize=(8.0, 5.2), dpi=170)

    for row in rows:
        fpos = row["positron_ratio"]
        nu, sed = load_isotropic_sed(Path(row["spectrum"]))
        ax.loglog(nu, sed, lw=2.0, label=fr"$f_+={fpos:g}$")

    nu_ref = rows[0]["nu_used_hz"]
    ax.axvline(nu_ref, color="0.35", ls="--", lw=1.1, label=fr"ref $\nu$={nu_ref:.3e} Hz")
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel(r"angle-averaged $\nu L_\nu$ [$L_\odot$]")
    ax.set_title("Pair sweep SED overlay (same dump, same seed)")
    ax.grid(alpha=0.25, which="both")
    ax.legend(frameon=False)
    fig.tight_layout()
    out = FIG_DIR / "pair_sed_overlay_nuLnu_Lsun.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_relative_scalings(rows):
    fpos = np.array([r["positron_ratio"] for r in rows])
    flux_scale = np.array([r["flux_scale_vs_pos0"] for r in rows])
    tau_s_mean_scale = np.array([r["tau_scatt_mean_scale_vs_pos0"] for r in rows])
    tau_s_w_scale = np.array([r["tau_scatt_weighted_scale_vs_pos0"] for r in rows])

    fig, ax = plt.subplots(figsize=(7.6, 4.8), dpi=170)
    ax.plot(fpos, flux_scale, "o-", lw=2, ms=6, label="flux scale")
    ax.plot(fpos, tau_s_mean_scale, "s-", lw=2, ms=6, label="tau_scatt mean scale")
    ax.plot(fpos, tau_s_w_scale, "^-", lw=2, ms=6, label="tau_scatt weighted scale")
    ax.set_xlabel("positron_ratio")
    ax.set_ylabel("scale relative to positron_ratio=0")
    ax.set_title("Relative pair-scaling checks")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False)
    fig.tight_layout()
    out = FIG_DIR / "pair_scaling_relative.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def plot_absolute_metrics(rows):
    fpos = np.array([r["positron_ratio"] for r in rows])
    flux_jy = np.array([r["flux_jy"] for r in rows])
    tau_s_mean = np.array([r["tau_scatt_mean"] for r in rows])
    tau_s_w = np.array([r["tau_scatt_weighted"] for r in rows])

    fig, ax1 = plt.subplots(figsize=(7.8, 4.9), dpi=170)
    l1 = ax1.plot(fpos, flux_jy, "o-", lw=2, ms=6, color="#1f77b4", label="flux [Jy]")
    ax1.set_xlabel("positron_ratio")
    ax1.set_ylabel("flux at ref nu [Jy]", color="#1f77b4")
    ax1.tick_params(axis="y", labelcolor="#1f77b4")
    ax1.grid(alpha=0.25)

    ax2 = ax1.twinx()
    l2 = ax2.plot(fpos, tau_s_mean, "s--", lw=1.8, ms=5, color="#2ca02c", label="tau_scatt mean")
    l3 = ax2.plot(fpos, tau_s_w, "^--", lw=1.8, ms=5, color="#d62728", label="tau_scatt weighted")
    ax2.set_ylabel("scattering optical depth", color="0.2")

    lines = l1 + l2 + l3
    labels = [ln.get_label() for ln in lines]
    ax1.legend(lines, labels, loc="upper left", frameon=False)

    ax1.set_title("Absolute flux and tau trends vs pair fraction")
    fig.tight_layout()
    out = FIG_DIR / "pair_scaling_absolute.png"
    fig.savefig(out)
    plt.close(fig)
    return out


def main():
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    rows = read_rows(CSV_PATH)
    out1 = plot_sed_overlay(rows)
    out2 = plot_relative_scalings(rows)
    out3 = plot_absolute_metrics(rows)
    print(out1)
    print(out2)
    print(out3)


if __name__ == "__main__":
    main()
