import h5py
import numpy as np

def computeTauRhoV(hdf5_filename):

    #Open file.  All we need is the polarization.
    with h5py.File(hdf5_filename, 'r') as imageFile:
        imageCube = imageFile['pol'][()]

    #Stokes I is the 1st index, and the Faraday depth is the 5th for no good reason.
    stokesI = imageCube[:,:,0]
    faradayRotationDepth = imageCube[:,:,4]

    #Compute the intensity-weighted average Faraday depth.  
    #Note that the Faraday depth is intensity-weighted among pixels, but not within pixels.
    tauRhoV = np.sum(faradayRotationDepth * stokesI) / np.sum(stokesI)

    return tauRhoV
