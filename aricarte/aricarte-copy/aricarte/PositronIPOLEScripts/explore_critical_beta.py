from fit_Munits import *
import os

def makeCriticalBetaImages(simFile, nameBase='../ipole_output/test/test.h5', inclination=50.0, freq_Hz=230e9,\
    ipole='../ipole_critical_beta/ipole', fov=200.0, npixel=400, counterjet=0, rmax_geo=20, beta_crit_coefficient_list=[0.1,0.5,0.9], \
    beta_crit_list=[1.0/3.0,1.0,3.0], keepImageLocation='../ipole_output/optimizeMunit/', fluxGoal=2.0):

    if not os.path.isdir(keepImageLocation):
        os.system('mkdir ' + keepImageLocation)

    for beta_crit in beta_crit_list:
        for beta_crit_coefficient in beta_crit_coefficient_list:
            optimizeMunit([simFile], inclination, fluxGoal, ipoleExecutable=ipole, \
            freq_Hz=freq_Hz, fov=fov, npixel=npixel, counterjet=counterjet, rmax_geo=rmax_geo, beta_crit_coefficient=beta_crit_coefficient, \
            beta_crit=beta_crit, keepImage=True, keepImageLocation=keepImageLocation)

if __name__ == '__main__':
    library = '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/'
    simFiles = [library + name for name in ['MAD/a+0.94/384x192x192_IHARM/dump_00002000.h5', 'SANE/a+0.94/288x128x128_IHARM/dump_00002000.h5']]
    for simFile in simFiles:
        makeCriticalBetaImages(simFile, fov=200.0, npixel=400, inclination=50.0, keepImageLocation='../ipole_output/optimizeMunit_inc50/')
        makeCriticalBetaImages(simFile, fov=200.0, npixel=400, inclination=10.0, keepImageLocation='../ipole_output/optimizeMunit_inc10/')
        makeCriticalBetaImages(simFile, fov=200.0, npixel=400, inclination=90.0, keepImageLocation='../ipole_output/optimizeMunit_inc90/')
