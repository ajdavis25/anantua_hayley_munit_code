#!/usr/bin/env python3
import re, h5py
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.axes_grid1 import make_axes_locatable


def _infer_positron_ratio(imageFile: str):
    """try to infer positron ratio from filename or header"""
    stem = Path(imageFile).stem
    m = re.search(r'(?:^|_)pos(?P<val>[\d.]+)', stem)
    if m:
        return float(m.group("val"))
    try:
        with h5py.File(imageFile, "r") as H:
            for k in [
                "/header/positronRatio",
                "/header/positron_ratio",
                "/header/positron",
                "/header/electrons/positronRatio",
            ]:
                if k in H:
                    return float(H[k][()])
    except Exception:
        pass
    return float("nan")


def colorbar(mappable):
    """attach a colorbar to an axis"""
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def plotPositronTestFrame(imageFile, cpMax=1e-2, fractionalCircular=True):
    """reproduce creating_images.ipynb image logic"""

    with h5py.File(imageFile, "r") as hfp:
        dx = hfp["header"]["camera"]["dx"][()]
        dy = hfp["header"]["camera"]["dy"][()]
        dsource = hfp["header"]["dsource"][()]
        lunit = hfp["header"]["units"]["L_unit"][()]
        scale = hfp["header"]["scale"][()]

        fov_muas = dx / dsource * lunit * 2.06265e11
        evpa_0 = "W"
        if "evpa_0" in hfp["header"]:
            evpa_0 = hfp["header"]["evpa_0"][()]

        unpol = np.copy(hfp["unpol"]).transpose((1, 0)) * scale
        imagep = np.copy(hfp["pol"]).transpose((1, 0, 2)) * scale
        pixelSize = (
            dx
            * dy
            * (lunit / dsource * 2.06265e11) ** 2
            / (imagep.shape[0] * imagep.shape[1])
        )

        I = imagep[:, :, 0] / pixelSize
        Q = imagep[:, :, 1] / pixelSize
        U = imagep[:, :, 2] / pixelSize
        V = imagep[:, :, 3] / pixelSize
        fov_muas = 1 * fov_muas

    extent = [-fov_muas / 2, fov_muas / 2, -fov_muas / 2, fov_muas / 2]
    fig, axarr = plt.subplots(1, 2, figsize=(8, 4))
    ax1, ax2 = axarr

    # stokes I intensity map
    im1 = ax1.imshow(
        I,
        cmap="afmhot",
        vmin=0.0,
        vmax=np.max(I),
        origin="lower",
        extent=extent,
    )
    colorbar(im1)
    ax1.set_title("Stokes I [Jy $\\mu$as$^{-2}$]")

    # circular polarization map
    if fractionalCircular:
        cpfrac = 100.0 * V / I
        im2 = ax2.imshow(
            cpfrac,
            cmap="seismic",
            vmin=-cpMax,
            vmax=cpMax,
            origin="lower",
            extent=extent,
        )
        colorbar(im2)
        ax2.set_title("CP [%]")
    else:
        im2 = ax2.imshow(
            V,
            cmap="seismic",
            vmin=-cpMax / 100.0 * np.max(I),
            vmax=cpMax / 100.0 * np.max(I),
            origin="lower",
            extent=extent,
        )
        colorbar(im2)
        ax2.set_title("CP [Jy $\\mu$as$^{-2}$]")

    # compute evpa and polarization vectors
    evpa = (180.0 / np.pi) * 0.5 * np.arctan2(U, Q)
    if evpa_0 == "W":
        evpa += 90.0
        evpa[evpa > 90.0] -= 180.0

    npix = I.shape[0]
    xs = np.linspace(-fov_muas / 2, fov_muas / 2, npix)
    Xs, Ys = np.meshgrid(xs, xs)
    lpscal = np.max(np.sqrt(Q**2 + U**2))
    vxp = np.sqrt(Q**2 + U**2) * np.sin(evpa * np.pi / 180.0) / lpscal
    vyp = -np.sqrt(Q**2 + U**2) * np.cos(evpa * np.pi / 180.0) / lpscal
    skip = int(npix / 32)
    ax1.quiver(
        Xs[::skip, ::skip],
        Ys[::skip, ::skip],
        vxp[::skip, ::skip],
        vyp[::skip, ::skip],
        headwidth=1,
        headlength=1,
        width=0.005,
        color="#00ff00",
        units="width",
        scale=4,
        pivot="mid",
    )

    # axes formatting
    for axis in [ax1, ax2]:
        axis.set_aspect("equal")
        axis.set_xlim([-40, 40])
        axis.set_ylim([-40, 40])
        axis.set_xlabel("$\\mu$as")
        axis.set_xticks(np.linspace(-40, 40, 5))
        axis.set_yticks(np.linspace(-40, 40, 5))
    ax1.set_ylabel("$\\mu$as")

    # text labels
    positronRatio = _infer_positron_ratio(imageFile)
    bbox = {"boxstyle": "round", "facecolor": "wheat", "alpha": 0.8}
    ax2.text(
        0.05,
        0.05,
        r"$n_{pairs}/(n_-)_0=$" + f"{positronRatio:1.2f}",
        ha="left",
        va="bottom",
        transform=ax2.transAxes,
        fontsize=12,
    )
    ax2.text(
        0.05,
        0.95,
        f"V/I = {np.sum(V)/np.sum(I):1.2e}",
        ha="left",
        va="top",
        transform=ax2.transAxes,
        fontsize=12,
    )
    ax1.text(
        0.05,
        0.95,
        f"I={np.sum(I)*pixelSize:3.2e} Jy",
        ha="left",
        va="top",
        transform=ax1.transAxes,
        fontsize=12,
        color="white",
    )
    ax1.text(
        0.05,
        0.05,
        f"P/I={np.sqrt(np.sum(Q)**2 + np.sum(U)**2)/np.sum(I):1.2e}",
        ha="left",
        va="bottom",
        transform=ax1.transAxes,
        fontsize=12,
        color="white",
    )

    fig.tight_layout()
    outpath = imageFile.replace(".h5", ".png")
    fig.savefig(outpath, dpi=150)
    plt.close(fig)
    print(f"saved {outpath}")


def batch_make_images(h5_dir):
    """run over all .h5 files in a directory"""
    h5_dir = Path(h5_dir)
    for h5file in sorted(h5_dir.glob("*.h5")):
        if "output_" in h5file.name or "MAD_" in h5file.name or "SANE_" in h5file.name:
            plotPositronTestFrame(str(h5file), cpMax=0.1, fractionalCircular=False)


if __name__ == "__main__":
    batch_make_images("/work/vmo703/ipole_outputs")
