from __future__ import print_function
import os
import time
import numpy as np
from ipole_many_models import runIPOLE
import h5py

#This is a dictionary that returns the dump range for a given folder.
library = '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library/'
folderToDumpRange = {}

#MADs
folderToDumpRange[library + 'MAD/a+0.94/384x192x192_IHARM'] = [1000,2000]
folderToDumpRange[library + 'MAD/a+0.5/384x192x192_IHARM'] = [1000,2000]
folderToDumpRange[library + 'MAD/a0/384x192x192_IHARM'] = [1000,1999]
folderToDumpRange[library + 'MAD/a-0.5/384x192x192_IHARM'] = [1000,1800]
folderToDumpRange[library + 'MAD/a-0.94/384x192x192_IHARM'] = [1400,2000]

#SANEs
folderToDumpRange[library + 'SANE/a+0.94/288x128x128_IHARM'] = [600,1200]
folderToDumpRange[library + 'SANE/a+0.5/288x128x128_IHARM'] = [600,1100]
folderToDumpRange[library + 'SANE/a0/288x128x128_IHARM'] = [1000,2000]
folderToDumpRange[library + 'SANE/a-0.5/288x128x128_IHARM'] = [1000,1600]
folderToDumpRange[library + 'SANE/a-0.94/288x128x128_IHARM'] = [1200,1800]

def computeTotalFlux(imageFile):

	"""
	Given an image, return the total flux density in Jy.
	"""

	with h5py.File(imageFile, 'r') as myfile:
		unscaledFlux = myfile['unpol'][()]
		scale = myfile['header']['scale'][()]

	#There shouldn't be any nans, but I don't want a single rogue pixel to completely crash everything.
	return np.nansum(unscaledFlux) * scale

def computeAverageTotalFlux(dumpFiles, Munit, thetacam, temporaryImageFileName='../ipole_output/optimizeMunit.h5', freq_Hz=230e9, MBH=4.14e6, npixel=400, \
	ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, phicam=0.0, fov=200.0, rmax_geo=50, parameterFileName='./runIPOLE.par', \
	counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, kappa_kappa=3.5, emission_type=4, unpol=True, Rhigh=20.0):

	fluxes_to_average = []
	for dumpFile in dumpFiles:
		#Make an image
		runIPOLE(dumpFile, temporaryImageFileName, Munit, freq_Hz=freq_Hz, MBH=MBH, npixel=npixel, ipoleExecutable=ipoleExecutable, \
		dsource=dsource, thetacam=thetacam, phicam=phicam, fov=fov, rmax_geo=rmax_geo, parameterFileName=parameterFileName, counterjet=counterjet, \
		target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, unpol=unpol, kappa_slope=kappa_kappa, \
		emission_type=emission_type, Rhigh=Rhigh)

		#Compute the total flux of that image.
		fluxes_to_average.append(computeTotalFlux(temporaryImageFileName))

	return np.mean(fluxes_to_average)

def optimizeMunit(dumpFiles, thetacam, fluxGoal, fractionalTolerance=0.05, Munit_guess=1e17, temporaryImageFileName='../ipole_output/optimizeMunit.h5', \
	freq_Hz=230e9, MBH=4.14e6, npixel=400, ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, \
	phicam=0.0, fov=200.0, rmax_geo=50, parameterFileName='./runIPOLE.par', counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, \
	kappa_kappa=3.5, emission_type=4, Rhigh=20.0):

	"""
	Use the secant method to obtain the optimal Munit for a list of files.
	"""

	#First, compute two fluxes for two Munits, separated arbitrarily around Munit_guess
	Munit_trial_list = [Munit_guess / np.sqrt(10), Munit_guess * np.sqrt(10)]
	flux_trial_list = []
	for Munit in Munit_trial_list:
		flux_trial_list.append(computeAverageTotalFlux(dumpFiles, Munit, thetacam, temporaryImageFileName=temporaryImageFileName, freq_Hz=freq_Hz, MBH=MBH, \
		npixel=npixel, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, fov=fov, rmax_geo=rmax_geo, parameterFileName=parameterFileName, \
		counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, kappa_kappa=kappa_kappa, emission_type=emission_type, \
		Rhigh=Rhigh))

	#Now, we start applying the secant method.
	iteration = 2
	M0 = Munit_trial_list[0] 
	M1 = Munit_trial_list[1]
	F0 = flux_trial_list[0]
	F1 = flux_trial_list[1]
	while True:
		iteration += 1

		#Compute slope and produce a new values.
		alpha = np.log(F1/F0) / np.log(M1/M0)
		A0 = F0 / np.power(M0, alpha)
		Munit = np.exp(np.log(fluxGoal/A0)/alpha)
		flux = computeAverageTotalFlux(dumpFiles, Munit, thetacam, temporaryImageFileName=temporaryImageFileName, freq_Hz=freq_Hz, MBH=MBH, \
		npixel=npixel, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, fov=fov, rmax_geo=rmax_geo, parameterFileName=parameterFileName, \
		counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, kappa_kappa=kappa_kappa, emission_type=emission_type, \
		Rhigh=Rhigh)
		Munit_trial_list.append(Munit)
		flux_trial_list.append(flux)
		if flux < fluxGoal:
			M0 = Munit
			F0 = flux
		else:
			M1 = Munit
			F1 = flux
		if np.abs(flux - fluxGoal) / fluxGoal <= fractionalTolerance:
			break
	return Munit_trial_list, flux_trial_list

def addToMunitTable(folder, thetacam, Munit_table_name, fluxGoal=2.0, fractionalTolerance=0.05, Munit_guess=1e17, temporaryImageFolderName='../ipole_output/optimizeMunit', \
	dumpSamples=100, dumpRange=[1000,2000], freq_Hz=230e9, MBH=4.14e6, npixel=400, ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, \
	phicam=0.0, fov=200.0, rmax_geo=50, parameterFileName=None, counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, keepImage=False, \
	keepImageLocation='/bd4/eht/Ricarte_SgrA/230GHz/', overwrite=False, kappa_kappa=3.5, emission_type=4, Rhigh=20.0, formatting='Rhigh'):

	"""
	Compute an Munit for a given inclination and folder of dumpfiles.  Add to the table.
	"""

	if not os.path.isdir(temporaryImageFolderName):
		os.system('mkdir '+temporaryImageFolderName)

	#Get a list of files.
	if folder[-1] != '/':
		folder += '/'
	availableDumpFiles = np.array([file for file in os.listdir(folder) if 'dump' in file])
	availableDumpNumbers = np.array([int(file.split('_')[1].split('.')[0]) for file in availableDumpFiles])
	desiredDumpNumbers = np.linspace(dumpRange[0], dumpRange[1], np.min((dumpSamples,dumpRange[1]-dumpRange[0]+1))).astype(int)
	isToBeUsed = np.in1d(availableDumpNumbers, desiredDumpNumbers)
	usedDumpFiles = np.array([folder + file for file in availableDumpFiles[isToBeUsed]])
	usedDumpNumbers = availableDumpNumbers[isToBeUsed]
	order = np.argsort(usedDumpNumbers)
	usedDumpFiles = usedDumpFiles[order]
	usedDumpNumbers = usedDumpNumbers[order]

	#Make sure we haven't done this one already.
	with open(Munit_table_name, 'r') as myfile:
		existingTable = np.loadtxt(myfile, dtype=str)
	if existingTable.shape[0] > 0:
		existingTable = np.atleast_2d(existingTable)
		existingFilenames = existingTable[:,0]
		if formatting == 'Rhigh':
			existingRhigh = existingTable[:,2].astype(float)
			existingInclinations = existingTable[:,3].astype(float)
			comboExists = np.any((existingFilenames == usedDumpFiles[-1]) & (existingRhigh == Rhigh) & (existingInclinations == thetacam))
		elif formatting == 'criticalBeta':
			existingBetaCritCoefficients = existingTable[:,2].astype(float)
			existingBetaCrits = existingTable[:,3].astype(float)
			existingInclinations = existingTable[:,4].astype(float)
			comboExists = np.any((existingFilenames == usedDumpFiles[-1]) & (existingBetaCritCoefficients == beta_crit_coefficient) & (existingBetaCrits == beta_crit) & (existingInclinations == thetacam))

		if comboExists:
			if overwrite:
				print("This combination exists in the table, but we're doing it anyway.")
			else:
				print("This combination exists in the table, so we're skipping it.")
				return

	#Naming scheme
	if temporaryImageFolderName[-1] != '/':
		temporaryImageFolderName += '/'
	outname = 'image_'
	splitFolderName = folder.split('/')
	if 'MAD' in splitFolderName:
		outname += 'M'
		startingIndex = splitFolderName.index('MAD')
	elif 'SANE' in splitFolderName:
		outname += 'S'
		startingIndex = splitFolderName.index('SANE')
	#Spin
	outname += splitFolderName[startingIndex+1]
	#Dump
	outname += '_{0:04d}'.format(usedDumpNumbers[-1])
	#Inclination
	outname += '_{0:d}'.format(int(thetacam))
	#phi?
	outname += '_{0:d}'.format(int(phicam))
	#Frequency
	outname += '_{0:1.3e}'.format(freq_Hz)
	#Mass
	outname += '_{0:1.3e}'.format(MBH)
	#Munit
	outname += '_MunitTBD'
	#Distance
	outname += '_{0:1.3e}'.format(dsource)
	if formatting == 'Rhigh':
		#Rhigh
		outname += '_{0:3.1f}'.format(Rhigh)
	elif formatting == 'criticalBeta':
		#Critical beta coefficient
		outname += '_{0:1.3f}'.format(beta_crit_coefficient)
		#Critical beta
		outname += '_{0:1.3f}'.format(beta_crit)
	#nx
	outname += '_{0:d}'.format(npixel)
	#ny
	outname += '_{0:d}'.format(npixel)
	#And finally, the suffix.
	outname += '.h5'
	temporaryImageFileName = temporaryImageFolderName + outname

	#The actual computation, where Munit is fit.
	Munit_trials, flux_trials = optimizeMunit(usedDumpFiles, thetacam, fluxGoal, fractionalTolerance=fractionalTolerance, Munit_guess=Munit_guess, \
	temporaryImageFileName=temporaryImageFileName, freq_Hz=freq_Hz, MBH=MBH, npixel=npixel, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, \
	fov=fov, rmax_geo=rmax_geo, parameterFileName=parameterFileName, counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, \
	beta_crit=beta_crit, kappa_kappa=kappa_kappa, emission_type=emission_type, Rhigh=Rhigh)
	print("Obtained an Munit of {0:1.3e} after {1:d} iterations.".format(Munit_trials[-1], len(Munit_trials)))

	if keepImage:
		if keepImageLocation[-1] != '/':
			keepImageLocation += '/'
		if not os.path.isfile(outname):
			os.system('mv '+temporaryImageFileName+' '+keepImageLocation + outname.replace('MunitTBD', '{0:1.5e}'.format(Munit_trials[-1])))

	with open(Munit_table_name, 'a') as tableFile:
		entry = '\n' + usedDumpFiles[-1] + '\t' + str(Munit_trials[-1])
		if formatting == 'Rhigh':
			entry += '\t' + '{0:3.2f}'.format(Rhigh)
		elif formatting == 'criticalBeta':
			entry += '\t' + '{0:3.2f}'.format(beta_crit_coefficient) + '\t' + '{0:3.2f}'.format(beta_crit)
		entry += '\t' + '{0:3.2f}'.format(thetacam)
		print(entry)
		tableFile.write(entry)

def initializeMunitTable(filename, formatting='Rhigh'):

	if os.path.isfile(filename):
		print("Will not initialize.  This file already exists.")
		return

	with open(filename, 'w') as myfile:
		if formatting == 'Rhigh':
			header = '#dumpPath' + '\t' + 'Munit' + '\t' + 'Rhigh' + '\t' + 'thetacam'
		elif formatting == 'criticalBeta':
			header = '#dumpPath' + '\t' + 'Munit' + '\t' + 'beta_crit_coefficient' + '\t' + 'beta_crit' + '\t' + 'thetacam'
		myfile.write(header)

def createEntireMunitTable(outputName, inclinationList=[50], dumpSamples=100, npixel=400, fov=200, Munit_guess=1e17, extraTextAfterKey='/', foldersToDo=None, \
	MBH=4.14e6, dsource=8.127e3, kappa_kappa=3.5, emission_type=4, ipoleExecutable='../ipole_critical_beta/ipole', flipRetrograde=True, fluxGoal=2.0, \
	formatting='Rhigh', Rhigh_list=[1.,10.,20.,40.,80.,160.]):

	initialTime = time.time()

	initializeMunitTable(outputName, formatting=formatting)
	if foldersToDo is None:
		folders = folderToDumpRange.keys()
	else:
		folders = foldersToDo

	if formatting == 'Rhigh':
		for folder in folders:
			dumpRange = folderToDumpRange[folder]
			for inclination in inclinationList:
				thetacam = inclination
				if flipRetrograde:
					if 'a-' in folder:
						thetacam = 180.0-inclination
				for Rhigh in Rhigh_list:
					#extraTextAfterKey should be set to '/dumps/' on Illinois, but not on Harvard.
					addToMunitTable(folder + extraTextAfterKey, thetacam, outputName, dumpSamples=dumpSamples, npixel=npixel, Munit_guess=Munit_guess, fov=fov, \
					dumpRange=dumpRange, MBH=MBH, dsource=dsource, kappa_kappa=kappa_kappa, emission_type=emission_type, ipoleExecutable=ipoleExecutable, \
					fluxGoal=fluxGoal, Rhigh=Rhigh, formatting=formatting)
	elif formatting == 'criticalBeta':
		for folder in folders:
			dumpRange = folderToDumpRange[folder]
			for inclination in inclinationList:
				thetacam = inclination
				if flipRetrograde:
					if 'a-' in folder:
						thetacam = 180.0-inclination
				#extraTextAfterKey should be set to '/dumps/' on Illinois, but not on Harvard.
				addToMunitTable(folder + extraTextAfterKey, thetacam, outputName, dumpSamples=dumpSamples, npixel=npixel, Munit_guess=Munit_guess, fov=fov, \
				dumpRange=dumpRange, MBH=MBH, dsource=dsource, kappa_kappa=kappa_kappa, emission_type=emission_type, ipoleExecutable=ipoleExecutable, \
				fluxGoal=fluxGoal, formatting=formatting)

	finalTime = time.time()
	print("Finished in {0:1.3e} hours.".format((finalTime-initialTime)/3600))

if __name__ == '__main__':
	import sys
	keyIndexToDo = int(sys.argv[1])
	foldersToDo = [folderToDumpRange.keys()[keyIndexToDo]]

	#Critical beta, Sgr A*
	'''
	Munit_table_name = 'Munit_table_critical_beta.txt'
	inclinationList = [10,90]
	createEntireMunitTable(Munit_table_name, inclinationList=inclinationList, dumpSamples=100, npixel=400, foldersToDo=foldersToDo)
	'''

	#Kappa, M87*
	'''
	Munit_table_name = './Munit_tables/Munit_table_kappa7.0_M87.txt'
	inclinationList = [163.]
	createEntireMunitTable(Munit_table_name, inclinationList=inclinationList, dumpSamples=11, npixel=320, fov=160, foldersToDo=foldersToDo, fluxGoal=0.5, \
	MBH=6.2e9, dsource=16.9e6, kappa_kappa=7.0, emission_type=2, ipoleExecutable='../ipole_versions/ipole_kappa/ipole', formatting='Rhigh', Rhigh_list=[1.,10.,40.,160.], \
	Munit_guess=1e25)
	'''

	#Kappa, Sgr A*
	kappaList = [3.5, 5.0, 7.0]
	for kappa in kappaList:
		tableName = './Munit_tables/Munit_table_kappa' + str(kappa) + '_SgrA.txt'
		inclinationList = [10.0, 50.0, 90.0]
		createEntireMunitTable(tableName, inclinationList=inclinationList, dumpSamples=11, npixel=400, fov=200, foldersToDo=foldersToDo, fluxGoal=2.0, \
		MBH=4.14e6, dsource=8.127e3, kappa_kappa=kappa, emission_type=2, ipoleExecutable='../ipole_versions/ipole_kappa/ipole', formatting='Rhigh', Rhigh_list=[1.,10.,40.,160.], \
		Munit_guess=1e17)
