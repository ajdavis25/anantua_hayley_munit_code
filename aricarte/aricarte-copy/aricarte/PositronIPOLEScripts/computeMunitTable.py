"""
Compute Munit for a list of models.
"""

from fitMunit import *
import numpy as np

def computeMunitTable(dumpList, RhighList, thetaList, goalFlux=0.5, differenceTolerance=0.005, \
        output='Munit_table.txt', MBH=6.5e9, DBH=16.9e6, freqcgs=230e9, guessBracket=(22,26), phicam=0):

    with open(output, 'w') as myfile:
        header = "GRMHD Dump\tM_unit\tRhigh\ttheta\n"
        myfile.write(header)
        for i in range(len(dumpList)):
            value = fitMunit(dumpList[i], Rhigh=RhighList[i], phicam=phicam, thetacam=thetaList[i], MBH=MBH, DBH=DBH, guesses=guessBracket, \
            freqcgs=freqcgs, goalFlux=goalFlux, differenceTolerance=differenceTolerance)
            myfile.write(dumpList[i] + '\t' + '{0:3.3e}\t{1:3.2f}\t{2:3.2f}\n'.format(value, RhighList[i], thetaList[i]))

if __name__ == '__main__':
    '''
    dumpList = ['/bd3/eht/GRMHD/SANE/a-0.94/288x128x128_IHARM/dumps/dump_00002000.h5', '/bd3/eht/GRMHD/SANE/a-0.94/288x128x128_IHARM/dumps/dump_00002000.h5', \
            '/bd3/eht/GRMHD/SANE/a+0.94/288x128x128_IHARM/dumps/dump_00002000.h5', '/bd3/eht/GRMHD/MAD/a-0.5/384x192x192_IHARM/dumps/dump_00002000.h5', \
            '/bd3/eht/GRMHD/MAD/a-0.5/384x192x192_IHARM/dumps/dump_00002000.h5', '/bd3/eht/GRMHD/MAD/a+0.94/384x192x192_IHARM/dumps/dump_00002000.h5', \
            '/bd3/eht/GRMHD/MAD/a+0.94/384x192x192_IHARM/dumps/dump_00002000.h5']
    '''
    dumpList = ['/n/home11/aricarte/projects/grmhd_library/SANE/a-0.94/288x128x128_IHARM/dump_00002000.h5', '/n/home11/aricarte/projects/grmhd_library/SANE/a-0.94/288x128x128_IHARM/dump_00002000.h5', \
            '/n/home11/aricarte/projects/grmhd_library/SANE/a+0.94/288x128x128_IHARM/dump_00002000.h5', '/n/home11/aricarte/projects/grmhd_library/MAD/a-0.5/384x192x192_IHARM/dump_00002000.h5', \
            '/n/home11/aricarte/projects/grmhd_library/MAD/a-0.5/384x192x192_IHARM/dump_00002000.h5', '/n/home11/aricarte/projects/grmhd_library/MAD/a+0.94/384x192x192_IHARM/dump_00002000.h5', \
            '/n/home11/aricarte/projects/grmhd_library/MAD/a+0.94/384x192x192_IHARM/dump_00002000.h5']
    RhighList = [10, 80, 160, 20, 160, 20, 160]
    #thetaList = [17, 17, 163, 17, 17, 163, 163]
    thetaList = [60, 60, 60, 60, 60, 60, 60]

    #Sgr A* mass and distance from Boehle et al. (2016)
    output = './Munit_table_sgra.txt'
    computeMunitTable(dumpList, RhighList, thetaList, output=output, MBH=4.02e6, DBH=7.86e3, goalFlux=3.4)
