import numpy as np
import h5py
import os
from scipy.optimize import brentq

def computeFluxToGoalRatio(logMunit, dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux):

    M_unit = 10.0**logMunit

    #First, run IPOLE
    with open('./fitMunit.par', 'w') as myfile:
        myfile.write('#Created automatically by fitMunit. \n \n')
        myfile.write('dump '+dumpfile+'\n')
        myfile.write('outfile ../temporary_output/fitMunit.h5\n')
        myfile.write('thetacam '+str(thetacam)+'\n')
        myfile.write('phicam '+str(phicam)+'\n')
        myfile.write('freqcgs '+str(freqcgs)+'\n')
        myfile.write('MBH {0:3.3e}'.format(MBH)+'\n')
        myfile.write('dsource {0:3.3e}'.format(DBH)+'\n')
        myfile.write('M_unit '+str(M_unit)+'\n')
        myfile.write('trat_large '+str(Rhigh)+'\n')
        myfile.write('nx 160\n')
        myfile.write('ny 160\n')
    os.system("../ipole_fewerapprox/ipole -par ./fitMunit.par -unpol")

    #Then read the output and extract the flux.
    with h5py.File('../temporary_output/fitMunit.h5', 'r') as hfp:
        I = hfp['unpol'][:,:]
        scale = hfp['header']['scale'][()]

    print("logMunit={0:3.4f} implies I={1:3.4f}".format(logMunit, np.sum(I) * scale))
    return np.log10(np.sum(I) * scale / goalFlux)

def fitMunit(dumpfile, differenceTolerance=0.01, output='../temporary_output/fitMunit.h5', \
    thetacam=163, phicam=0, freqcgs=230.0e9, MBH=6.5e9, guesses=[24,26], absoluteBounds=[10,40], \
    Rhigh=3, goalFlux=0.5, DBH=16.9e6):

    #Initialize
    window = [guesses[0], guesses[1]]
    logFluxRatios = [computeFluxToGoalRatio(window[0], dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux), \
            computeFluxToGoalRatio(window[1], dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux)]

    #Check how the initial values are doing.
    diff = np.abs(window[1] - window[0])
    while True:
        if (logFluxRatios[0] < 0) & (logFluxRatios[1] < 0):
            #If both fluxes are too low, there are two possibilities...
            if logFluxRatios[1] > logFluxRatios[0]:
                #Just needs more accretion.  Go up.
                if window[1] >= absoluteBounds[1]:
                    raise ValueError("Can't raise the Munit any higher.")
                else:
                    window = [window[1], np.min([window[1]+diff,absoluteBounds[1]])]
                    logFluxRatios = [logFluxRatios[1], computeFluxToGoalRatio(window[1], dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux)]
            elif logFluxRatios[1] < logFluxRatios[0]:
                #Too optically thick.  Go down.
                if window[0] <= absoluteBounds[0]:
                    raise ValueError("Can't decrease the Munit any lower.")
                else:
                    window = [np.max([absoluteBounds[0],window[0]-diff]), window[0]]
                    logFluxRatios = [computeFluxToGoalRatio(window[0], dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux), logFluxRatios[0]]
        elif (logFluxRatios[0] > 0) & (logFluxRatios[1] > 0):
            #If both fluxes are too high, always go down.
            if window[0] <= absoluteBounds[0]:
                raise ValueError("Can't decrease the Munit any lower.")
            else:
                window = [np.max([absoluteBounds[0],window[0]-diff]), window[0]]
                logFluxRatios = [computeFluxToGoalRatio(window[0], dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux), logFluxRatios[0]]
        else:
            #Ready to enter secant method
            break
    
    while True:
        #Do the bisection method.  Secant was going crazy.
        middle = 0.5 * (window[0] + window[1])
        logFluxRatioMiddle = computeFluxToGoalRatio(middle, dumpfile, thetacam, phicam, freqcgs, MBH, DBH, Rhigh, goalFlux)
        if np.abs(logFluxRatios[1]) < differenceTolerance:
            #Success!
            break
        else:
            if logFluxRatioMiddle > 0:
                window = [window[0], middle]
                logFluxRatios = [logFluxRatios[0], logFluxRatioMiddle]
            elif logFluxRatioMiddle < 0:
                window = [middle, window[1]]
                logFluxRatios = [logFluxRatioMiddle, logFluxRatios[1]]

    return 10.0**middle

if __name__ == '__main__':
    #Fitting parameters
    goalFlux = 0.5
    differenceTolerance = 0.005
    guesses = [27,30]

    #IPOLE parameters
    freqcgs = 230e9
    MBH = 6.5e9
    phicam = 0

    '''
    #MAD a+0.94
    dumpfile = '/bd3/eht/GRMHD/MAD/a+0.94/384x192x192_IHARM/dumps/dump_00002000.h5'
    Rhigh = 160
    thetacam = 163
    '''

    '''
    #MAD a-0.5
    dumpfile = '/bd3/eht/GRMHD/MAD/a-0.5/384x192x192_IHARM/dumps/dump_00002000.h5'
    Rhigh = 20
    thetacam = 17
    '''

    '''
    #SANE a+0.94
    dumpfile = '/bd3/eht/GRMHD/SANE/a+0.94/288x128x128_IHARM/dumps/dump_00002000.h5'
    Rhigh = 160
    thetacam = 163
    '''

    #SANE a-0.94
    dumpfile = '/bd3/eht/GRMHD/SANE/a-0.94/288x128x128_IHARM/dumps/dump_00002000.h5'
    Rhigh = 10
    thetacam = 17

    print fitMunit(dumpfile, Rhigh=Rhigh, phicam=phicam, thetacam=thetacam, MBH=MBH, guesses=guesses, \
    freqcgs=freqcgs, goalFlux=goalFlux, differenceTolerance=differenceTolerance)
