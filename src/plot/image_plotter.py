import re, h5py, matplotlib, matplotlib.pyplot as plt, numpy as np
from pathlib import Path
matplotlib.use('Agg')


def _infer_positron_ratio(imageFile: str):
    stem = Path(imageFile).stem
    m = re.search(r'(?:^|_)pos(?P<val>[\d.]+)$', stem)
    if m:
        return float(m.group('val'))
    last = stem.split('_')[-1]
    try:
        return float(last)
    except ValueError:
        pass
    try:
        with h5py.File(imageFile, "r") as H:
            return float(H["header"]["electrons"]["positronRatio"][()])
    except Exception:
        return float("nan")


def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)


def plotPositronTestFrame(imageFile, intensityMax=3e-3, cpMax=1e-2, output=None, EVPA_CONV="EofN", fractionalCircular=True):
    with h5py.File(imageFile, 'r') as hfp:
        dx = hfp['header']['camera']['dx'][()]
        dy = hfp['header']['camera']['dy'][()]
        dsource = hfp['header']['dsource'][()]
        lunit = hfp['header']['units']['L_unit'][()]
        fov_muas = dx / dsource * lunit * 2.06265e11
        scale = hfp['header']['scale'][()]
        evpa_0 = hfp['header'].get('evpa_0', 'W')
        imagep = np.copy(hfp['pol']).transpose((1,0,2)) * scale
        pixelSize = dx * dy * (lunit / dsource * 2.06265e11)**2 / (imagep.shape[0] * imagep.shape[1])
        I, Q, U, V = [imagep[:,:,i] / pixelSize for i in range(4)]

    extent = [-fov_muas/2, fov_muas/2, -fov_muas/2, fov_muas/2]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8,4))

    im1 = ax1.imshow(I, cmap='afmhot', vmin=0., vmax=np.max(I), origin='lower', extent=extent)
    colorbar(im1)

    if fractionalCircular:
        cpfrac = 100. * V / I
        im2 = ax2.imshow(cpfrac, cmap='seismic', vmin=-cpMax, vmax=cpMax, origin='lower', extent=extent)
        ax2.set_title("CP [%]")
    else:
        im2 = ax2.imshow(V, cmap='seismic', vmin=-cpMax/100.0*intensityMax, vmax=cpMax/100.0*intensityMax, origin='lower', extent=extent)
        ax2.set_title("CP [Jy $\mu$as$^{-2}$]")
    colorbar(im2)

    evpa = (180./np.pi)*0.5*np.arctan2(U, Q)
