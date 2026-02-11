#!/usr/bin/env python3
import argparse
import re, h5py, matplotlib, numpy as np, matplotlib.pyplot as plt
from pathlib import Path
matplotlib.use('Agg')


root_dir = Path("/work/vmo703/ipole_outputs/M87")


def infer_positron_ratio(path):
    # try filename like ..._pos0_ or ..._pos1_
    m = re.search(r'_pos(\d)(?:_|\.|$)', path)
    if m:
        return float(m.group(1))
    # fallback: read from HDF5 header if present
    try:
        with h5py.File(path, "r") as f:
            for k in ["/header/positronRatio", "/header/positron_ratio", "/header/positron"]:
                if k in f:
                    return float(f[k][()])
    except Exception:
        pass
    return None


def _infer_positron_ratio(imageFile: str):
    stem = Path(imageFile).stem  # no .h5
    # case 1: our pattern ..._pos<val>
    m = re.search(r'(?:^|_)pos(?P<val>[\d.]+)$', stem)
    if m:
        return float(m.group('val'))
    # case 2: last token is purely numeric ..._<val>
    last = stem.split('_')[-1]
    try:
        return float(last)
    except ValueError:
        pass
    # case 3: try to read from file (if the build writes it)
    try:
        with h5py.File(imageFile, "r") as H:
            return float(H["header"]["electrons"]["positronRatio"][()])
    except Exception:
        return float("nan")  # unknown but don't crash


def colorbar(mappable):
    """ the way matplotlib colorbar should have been implemented """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def plotPositronTestFrame(
    imageFile,
    intensityMax=3e-3,
    cpMax=1e-2,
    output=None,
    EVPA_CONV="EofN",
    fractionalCircular=True,
    interpolation="nearest",
    dpi=300,
):

    #Open the IPOLE output and extract relevant data.
    
    with h5py.File(imageFile, 'r') as hfp:
        dx = hfp['header']['camera']['dx'][()]
        dy = hfp['header']['camera']['dy'][()]
        dsource = hfp['header']['dsource'][()]
        lunit = hfp['header']['units']['L_unit'][()]
        fov_muas = dx / dsource * lunit * 2.06265e11
        scale = hfp['header']['scale'][()]
        evpa_0 = 'W'
        if 'evpa_0' in hfp['header']:
            evpa_0 = hfp['header']['evpa_0'][()]
        unpol = np.copy(hfp['unpol']).transpose((1,0)) * scale
        imagep = np.copy(hfp['pol']).transpose((1,0,2)) * scale
        pixelSize = dx * dy * (lunit / dsource * 2.06265e11)**2 / (imagep.shape[0] * imagep.shape[1])
        I = imagep[:,:,0] / pixelSize
        Q = imagep[:,:,1] / pixelSize
        U = imagep[:,:,2] / pixelSize
        V = imagep[:,:,3] / pixelSize
        fov_muas = 1*fov_muas # Find 86GHz FoV, 230GHz FoV in ipole_positron_Test_unlocked.py; Increasing fov_muas here zooms in                                         
    extent = [ -1*fov_muas/2, 1*fov_muas/2, -1*fov_muas/2, 1*fov_muas/2 ]

    #Initialize a plot with two panels.
    fig, axarr = plt.subplots(1, 2, figsize=(8,4))
    ax1 = axarr[0]
    ax2 = axarr[1]

    #Total intensity and linear polarization ticks.
    #Is = np.sum(I)
    im1 = ax1.imshow(
        I,
        cmap='afmhot',
        vmin=0.,
        vmax=np.max(I),
        origin='lower',
        extent=extent,
        interpolation=interpolation
    ) #vmax=intensityMax
    # intensityMax; vmax=np.max(I)
    colorbar(im1)
    #print(Is)
    #print('V/I = {0:1.2e}'.format(np.sum(V)/np.sum(I)))
    
    #Circular polarization fraction
    if fractionalCircular:
        cpfrac = 100.*V/I
        im2 = ax2.imshow(
            cpfrac,
            cmap='seismic',
            vmin=-cpMax,
            vmax=cpMax,
            origin='lower',
            extent=extent,
            interpolation=interpolation
        )
        colorbar(im2)
        ax2.set_title("CP [%]")
    else:
        im2 = ax2.imshow(
            V,
            cmap='seismic',
            vmin=-cpMax/100.0*intensityMax,
            vmax=cpMax/100.0*intensityMax,
            origin='lower',
            extent=extent,
            interpolation=interpolation
        )
        colorbar(im2)
        ax2.set_title("CP [Jy $\mu$as$^{-2}$]")

    #evpa
    evpa = (180./3.14159)*0.5*np.arctan2(U,Q)
    if evpa_0 == "W":
        evpa += 90.
        evpa[evpa > 90.] -= 180.
    if EVPA_CONV == "NofW":
        evpa += 90.
        evpa[evpa > 90.] -= 180.

    #quiver on intensity
    npix = I.shape[0]
    xs = np.linspace(-fov_muas/2,fov_muas/2,npix)
    Xs,Ys = np.meshgrid(xs,xs)
    lpscal = np.max(np.sqrt(Q*Q+U*U))
    vxp = np.sqrt(Q*Q+U*U)*np.sin(evpa*3.14159/180.)/lpscal
    vyp = -np.sqrt(Q*Q+U*U)*np.cos(evpa*3.14159/180.)/lpscal
    skip = int(npix/32)
    ax1.quiver(Xs[::skip,::skip],Ys[::skip,::skip],vxp[::skip,::skip],vyp[::skip,::skip],
      headwidth=1, headlength=1,
      width=0.005,
      color='#00ff00',
      units='width',
      scale=4,
      pivot='mid')

    #Formatting
    ax1.set_title("Stokes I [Jy $\mu$as$^{-2}$]")
    ax1.set_ylabel('$\mu$as')
    ax1.set_xlabel('$\mu$as')
    ax2.set_xlabel('$\mu$as')
    for axis in [ax1,ax2]:
        axis.set_aspect('equal')
        axis.set_xlim([-1*40,1*40])
        axis.set_ylim([-1*40,1*40])
        axis.set_xticks(np.linspace(-1*40,1*40,1*5))
        axis.set_yticks(np.linspace(-1*40,1*40,1*5))

    #Label.  Note that the positron fraction was included as part of the file name.  See some basic calculations here.
    positronRatio = _infer_positron_ratio(imageFile)
    bbox = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.8}
    ax2.text(0.05, 0.05, r'$n_{pairs}/(n_-)_0=$' + '{0:1.2f}'.format(positronRatio), ha='left', va='bottom', transform=ax2.transAxes, fontsize=12)
    ax2.text(0.05, 0.95, 'V/I = {0:1.2e}'.format(np.sum(V)/np.sum(I)), ha='left', va='top', transform=ax2.transAxes, fontsize=12)
    ax1.text(0.05, 0.95, 'I={0:3.2e} Jy'.format(np.sum(I) * pixelSize), ha='left', va='top', transform=ax1.transAxes, fontsize=12, color='white')
    #ax1.text(0.05, 0.95, 'I={0:3.2f} Jy'.format(np.sum(I) * pixelSize), ha='left', va='top', transform=ax1.transAxes, fontsize=12, color='white')
    ax1.text(0.05, 0.05, 'P/I={0:1.2e}'.format(np.sqrt(np.sum(Q)**2 + np.sum(U)**2)/np.sum(I)), ha='left', va='bottom', transform=ax1.transAxes, fontsize=12, color='white')
    print(positronRatio, '{0:1.2e}'.format(np.sum(V)/np.sum(I)))
    fig.tight_layout()
    fig.savefig(output or imageFile.replace(".h5",".png"), dpi=dpi)
    plt.close(fig)


def iter_h5_files(root: Path, include_root: bool = True, include_subdirs: bool = True):
    """Yield .h5 files in root_dir and/or recursively in each root_dir/* directory."""
    seen = set()

    # Files directly in root_dir (e.g. ~/M87/*.h5).
    if include_root:
        for h5_file in sorted(root.glob("*.h5")):
            if not h5_file.is_file():
                continue
            resolved = h5_file.resolve()
            if resolved in seen:
                continue
            seen.add(resolved)
            yield h5_file

    # Files under each child directory in root_dir, recursively.
    if include_subdirs:
        for child in sorted(root.glob("*")):
            if not child.is_dir() or child.name == "images":
                continue
            for h5_file in sorted(child.rglob("*.h5")):
                if not h5_file.is_file():
                    continue
                resolved = h5_file.resolve()
                if resolved in seen:
                    continue
                seen.add(resolved)
                yield h5_file


def parse_args():
    parser = argparse.ArgumentParser(
        description="Create polarization images from IPOLE .h5 outputs."
    )
    parser.add_argument(
        "--root-dir",
        type=Path,
        default=root_dir,
        help=f"Top-level directory to scan (default: {root_dir}).",
    )
    parser.add_argument(
        "--scope",
        choices=["all", "root", "subdirs"],
        default="all",
        help="Which .h5 files to process: all (default), root-level only, or subdirs only.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Regenerate PNGs even if output images already exist.",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="PNG output DPI (default: 300).",
    )
    parser.add_argument(
        "--interpolation",
        default="nearest",
        choices=["nearest", "none", "antialiased", "bilinear", "bicubic", "hanning", "hamming", "hermite", "kaiser", "quadric", "catrom", "gaussian", "bessel", "mitchell", "sinc", "lanczos"],
        help="Matplotlib imshow interpolation mode (default: nearest).",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    root = args.root_dir.expanduser().resolve()

    include_root = args.scope in ("all", "root")
    include_subdirs = args.scope in ("all", "subdirs")

    processed = 0
    skipped_existing = 0
    failed = 0

    for h5_file in iter_h5_files(root, include_root=include_root, include_subdirs=include_subdirs):
        image_output_dir = h5_file.parent / "images"
        image_output_dir.mkdir(exist_ok=True)
        output_image_path = image_output_dir / (h5_file.stem + ".png")

        if output_image_path.exists() and not args.overwrite:
            skipped_existing += 1
            continue

        try:
            print(f"processing {h5_file}...")
            plotPositronTestFrame(
                str(h5_file),
                cpMax=0.1,
                fractionalCircular=False,
                output=str(output_image_path),
                interpolation=args.interpolation,
                dpi=args.dpi
            )
            processed += 1
        except Exception as e:
            failed += 1
            print(f"failed to process {h5_file.name}: {e}")

    print(
        f"done: processed={processed}, skipped_existing={skipped_existing}, "
        f"failed={failed}, scope={args.scope}, root_dir={root}"
    )


if __name__ == "__main__":
    main()
