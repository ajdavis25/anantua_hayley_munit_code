import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from folderToMovie import *

def colorbar(mappable):
    """ the way matplotlib colorbar should have been implemented """
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    return fig.colorbar(mappable, cax=cax)

def plotPositronTestFrame(imageFile, intensityMax=3e-3, cpMax=20, output=None, EVPA_CONV="EofN", fractionalCircular=True):

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
    extent = [ -fov_muas/2, fov_muas/2, -fov_muas/2, fov_muas/2 ]

    #Initialize a plot with two panels.
    fig, axarr = plt.subplots(1, 2, figsize=(8,4))
    ax1 = axarr[0]
    ax2 = axarr[1]

    #Total intensity and linear polarization ticks.
    im1 = ax1.imshow(I, cmap='afmhot', vmin=0., vmax=intensityMax, origin='lower', extent=extent)
    colorbar(im1)

    #Circular polarization fraction
    if fractionalCircular:
        cpfrac = 100.*V/I
        im2 = ax2.imshow(cpfrac, cmap='seismic', vmin=-cpMax, vmax=cpMax, origin='lower', extent=extent)
        colorbar(im2)
        ax2.set_title("CP [%]")
    else:
        im2 = ax2.imshow(V, cmap='seismic', vmin=-cpMax/100.0*intensityMax, vmax=cpMax/100.0*intensityMax, origin='lower', extent=extent)
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
        axis.set_xlim([-40,40])
        axis.set_ylim([-40,40])
        axis.set_xticks(np.linspace(-40,40,5))
        axis.set_yticks(np.linspace(-40,40,5))

    #Label.  Note that the positron fraction was included as part of the file name.  See some basic calculations here.
    positronRatio = float(imageFile.split('_')[-1][:-3])
    bbox = {'boxstyle': 'round', 'facecolor': 'wheat', 'alpha': 0.8}
    ax2.text(0.05, 0.05, r'$n_{pairs}/(n_-)_0=$' + '{0:1.2f}'.format(positronRatio), ha='left', va='bottom', transform=ax2.transAxes, fontsize=12)
    ax2.text(0.05, 0.95, 'V/I = {0:1.2e}'.format(np.sum(V)/np.sum(I)), ha='left', va='top', transform=ax2.transAxes, fontsize=12)
    ax1.text(0.05, 0.95, 'I={0:3.2f} Jy'.format(np.sum(I) * pixelSize), ha='left', va='top', transform=ax1.transAxes, fontsize=12, color='white')
    ax1.text(0.05, 0.05, 'P/I={0:1.2e}'.format(np.sqrt(np.sum(Q)**2 + np.sum(U)**2)/np.sum(I)), ha='left', va='bottom', transform=ax1.transAxes, fontsize=12, color='white')

    fig.tight_layout()
    if output is None:
        fig.show()
    else:
        fig.savefig(output)
        plt.close(fig)

def makePositronTestMovie(imageFolder, movieName, temporaryFolder='../images/frames/', fractionalCircular=True, cpMax=20.0, intensityMax=3e-3):

    #Sanitize input by adding a slash to the end if necessary.
    if imageFolder[-1] != '/':
        imageFolder += '/'
    if temporaryFolder[-1] != '/':
       temporaryFolder += '/'
        
    #WARNING:  Destroys all frames in the temporary folder.
    os.system('rm '+temporaryFolder+'frame*png')

    #In the image folder, there's a series of images which contain the positron fraction in the file name.  Interpret them.
    listOfImages = os.listdir(imageFolder)
    positronRatios = [float(imageFile.split('_')[-1][:-3]) for imageFile in listOfImages]
    ordering = np.argsort(positronRatios)
    listOfImages = [listOfImages[order] for order in ordering]

    #Make a plot for each frame.
    for frame in range(len(listOfImages)):
        print ( 'Creating frame {0} of {1}.'.format(frame+1, len(listOfImages)) )#Added parentheses after download from Angelo
        plotPositronTestFrame(imageFolder + listOfImages[frame], output=temporaryFolder + 'frame{0:04}.png'.format(frame), \
        fractionalCircular=fractionalCircular, cpMax=cpMax, intensityMax=intensityMax)

    #Convert the folder into a movie.
    folderToMovie(temporaryFolder, movieName)

if __name__ == '__main__':
    #These are example commands that I used...
#positronTest_MADaPt94_RHigh20BetaEPt01MOffset5Pt0e24MSlope1Pt5e0SigmaTransition0Pt5SigmaCut2KHARMAcpMaxPt1.mp4
#    makePositronTestMovie('../ipole_output/MAD_a+0.94_R20.0_positrons/', '../movies/positronTest_MAD.mp40',)/n/home13/anantua/ipole/aricarte/e+Output
#    makePositronTestMovie('/n/home13/anantua/ipole/aricarte/e+Output/', 'movies/positronTest_MAD_BetaE0Pt01SigmaTransition1SigmaCut2.mp4', fractionalCircular=False, cpMax=4)
    makePositronTestMovie('/work/epc964/ipole_HarvardCannonHPC_copy/aricarte/e+Output/', 'movies/positronTest_MADa+Pt94_fPt5BetaC1MOffset7Pt5e24MSlope1Pt5e0SigmaTransition2SigmaCut2KHARMAcpMaxPt1T5000rMax50.mp4', fractionalCircular=False, cpMax=1e-1)#cpMax=4     

#positronTest_SANEa-Pt5_RHigh20MOffset8Pt0e27MSlope2Pt0e1SigmaTransitionPt5SigmaCut2KHARMAcpMax1.mp4
#    makePositronTestMovie('../ipole_output/SANE_a+0.94_R20.0_positrons/', '../movies/positronTest_SANE.mp4', fractionalCircular=False, cpMax=20, intensityMax=1e-3)
#    makePositronTestMovie('ipole_output/SANE_a+0.94_R20.0_positrons/', 'movies/positronTest_SANE.mp4', fractionalCircular=False, cpMax=20, intensityMax=1e-3)
#    makePositronTestMovie('/n/home13/anantua/ipole/aricarte/e+Output/', 'movies/positronTest_SANEa-Pt5_fPt5BetaC1BetaEPt01MOffset7Pt5e27MSlope2Pt0e1SigmaTransition0Pt5SigmaCut2KHARMAcpMax1rMax50.mp4',fractionalCircular=False, cpMax=1e-0, intensityMax=1e-3)#cpMax=20,intensityMax=1e-3
