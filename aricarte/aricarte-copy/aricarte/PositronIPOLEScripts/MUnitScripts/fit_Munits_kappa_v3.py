from __future__ import print_function
import os
import time
import numpy as np
import h5py
import subprocess
import time

#This is a dictionary that returns the dump range for a given folder.
library = '/n/holystore01/LABS/bhi/Lab/narayan_lab/grmhd_library_v3/'
folderToDumpRange = {}

grmhd_names = [library + name for name in \
'MAD/a+0.94/288x128x128_KHARMA_2041', \
'MAD/a+0.5/288x128x128_KHARMA_2041', \
'MAD/a0/288x128x128_KHARMA_2041', \
'MAD/a-0.5/288x128x128_KHARMA_2041', \
'MAD/a-0.94/288x128x128_KHARMA_2041', \
'SANE/a+0.94/288x128x128_KHARMA', \
'SANE/a+0.5/288x128x128_KHARMA', \
'SANE/a0/288x128x128_KHARMA', \
'SANE/a-0.5/288x128x128_KHARMA', \
'SANE/a-0.94/288x128x128_KHARMA']

#Same dump range for each model.
for name in grmhd_names:
	folderToDumpRange[name] = [5000,6000]

#Image specifications for each source.
sources = ['M87', 'SgrA']
npixel = {'M87': 240, 'SgrA': 300}        #Cut in half from the usual specifications due to only needing total flux.
fov = {'M87': 240, 'SgrA': 300}           #Expanded from the usual specifications due to kappa.
MBH = {'M87': 6.2e9, 'SgrA': 4.14e6}
dsource = {'M87': 16.9e6, 'SgrA': 8.127e3}
inclinationList = {'M87': [163], 'SgrA': [10,50,90]}
flipRetrograde = {'M87': True, 'SgrA': False}
fluxGoal = {'M87': 0.5, 'SgrA': 2.4}

#Kappa options.
kappaList = [3.5, 5.0, 7.0, 0.0]
emissionType_list = [2, 2, 2, 4]

def computeIPOLEFlux(inputSimulationFile, Munit=1e20, Rhigh=20, freq_Hz=230e9, MBH=6.2e9, npixel_max=129, ipoleExecutable='../ipole_versions/ipole_dev/ipole', \
    thetacam=163.0, phicam=0.0, fov=160.0, rmax_geo=50, counterjet=0, target_nturns=-1, positronRatio=0, \
    beta_crit_coefficient=0.5, beta_crit=1, dsource=16.9e6, unpol=False, emission_type=4, kappa_slope=3.5, adaptiveRefinement=True, npixel_min=17, printCommand=True):

	#IPOLE can take all of the free parameters on the command line.
	command = ipoleExecutable

	#Always do this unpolarized with quenched output.
	command += ' -unpol'
	command += ' -quench'

	#Other parameters.
	command += ' --dump='+inputSimulationFile
	command += ' --thetacam='+str(thetacam)
	command += ' --phicam='+str(phicam)
	command += ' --freqcgs='+str(freq_Hz)
	command += ' --MBH={0:2.2e}'.format(MBH)
	command += ' --dsource={0:2.2e}'.format(dsource)
	command += ' --M_unit='+str(Munit)
	command += ' --trat_large='+str(Rhigh)
	command += ' --nx='+str(npixel_max)
	command += ' --ny='+str(npixel_max)
	command += ' --fov='+str(fov)
	command += ' --counterjet='+str(counterjet)
	command += ' --rmax_geo='+str(rmax_geo)
	command += ' --target_nturns='+str(target_nturns)
	command += ' --positronRatio='+str(positronRatio)
	command += ' --beta_crit='+str(beta_crit)
	command += ' --beta_crit_coefficient='+str(beta_crit_coefficient)
	command += ' --emission_type='+str(emission_type)
	command += ' --kappa='+str(kappa_slope)
	if adaptiveRefinement:
		command += ' --nx_min='+str(npixel_min)
		command += ' --ny_min='+str(npixel_min)
		
	if printCommand:
		print(command)
	proc = subprocess.Popen(command.split(' '), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	
	#Reading the flux from the terminal output.  Assuming things have been printed in a very specific way.
	output = [ z for y in [ str(x)[2:-1].split("\\n") for x in proc.communicate() ] for z in y ]
	Ftot_line = [l for l in output if 'unpol xfer' in l][0]
	st = Ftot_line.split()
	Ftot_unpol = float(st[-2+st.index('unpol')][1:])
	print(" F = {0:1.3e} Jy for Munit = {1:1.3e}.".format(Ftot_unpol, Munit))
	return Ftot_unpol

def computeAverageTotalFlux(dumpFiles, Munit, thetacam, freq_Hz=230e9, MBH=4.14e6, npixel_max=129, \
	npixel_min=17, ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, phicam=0.0, fov=200.0, rmax_geo=50, \
	counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, model_kappa=3.5, emission_type=4, unpol=True, Rhigh=20.0, adaptiveRefinement=True):

	fluxes_to_average = []
	for dumpFile in dumpFiles:
		#Make an image
		totalFlux = computeIPOLEFlux(dumpFile, Munit=Munit, freq_Hz=freq_Hz, MBH=MBH, npixel_max=npixel_max, ipoleExecutable=ipoleExecutable, \
		dsource=dsource, thetacam=thetacam, phicam=phicam, fov=fov, rmax_geo=rmax_geo, counterjet=counterjet, \
		target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, unpol=unpol, kappa_slope=model_kappa, \
		emission_type=emission_type, Rhigh=Rhigh, npixel_min=npixel_min, adaptiveRefinement=adaptiveRefinement)

		#Compute the total flux of that image.
		fluxes_to_average.append(totalFlux)
	averageFlux = np.mean(fluxes_to_average)
	print("  Obtained an average flux of {0:1.3e} Jy for an Munit of {1:1.3e}.".format(averageFlux, Munit))

	return np.mean(fluxes_to_average)

def optimizeMunit(dumpFiles, thetacam, fluxGoal, fractionalTolerance=0.02, Munit_guess=1e19, \
	freq_Hz=230e9, MBH=4.14e6, npixel_max=129, npixel_min=17, ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, \
	phicam=0.0, fov=200.0, rmax_geo=50, counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, \
	model_kappa=3.5, emission_type=4, Rhigh=20.0, adaptiveRefinement=True, smartLimits=[1e10,1e40]):

	"""
	Use the secant method to obtain the optimal Munit for a list of files.
	"""

	iteration = 0
	solutionFound = False
	Munit_trial_list = []
	flux_trial_list = []

	#First, compute two fluxes for two Munits, separated arbitrarily around Munit_guess
	Munit_first_trials = [Munit_guess, Munit_guess*2]
	for Munit in Munit_first_trials:
		if not solutionFound:
			iteration += 1
			flux = computeAverageTotalFlux(dumpFiles, Munit, thetacam, freq_Hz=freq_Hz, MBH=MBH, \
			npixel_max=npixel_max, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, fov=fov, rmax_geo=rmax_geo, \
			counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, model_kappa=model_kappa, emission_type=emission_type, \
			Rhigh=Rhigh, npixel_min=npixel_min, adaptiveRefinement=adaptiveRefinement)
			Munit_trial_list.append(Munit)
			flux_trial_list.append(flux)
			if (flux >= fluxGoal / (1.0 + fractionalTolerance)) & (flux <= fluxGoal * (1.0 + fractionalTolerance)):
				solutionFound = True

	#Now, we start applying the secant method in log space.
	if not solutionFound:
		M0 = Munit_trial_list[0] 
		M1 = Munit_trial_list[1]
		F0 = flux_trial_list[0]
		F1 = flux_trial_list[1]
	while not solutionFound:
		iteration += 1

		#Compute the logarithmic slope.
		alpha = np.log(F1/F0) / np.log(M1/M0)
		if alpha > 0:
			#Propose a new Munit by guessing that there's a straight line in log-space: alpha*logM - logF = constant.
			A0 = F0 / np.power(M0, alpha)
			Munit = np.exp(np.log(fluxGoal/A0)/alpha)
		elif alpha <= 0:
			#This can happen if Munit is so high that you're decreasing the flux by making it more and more optically thick.
			#Try to scale it back by the same factor that you've overshot.
			if M0 < M1:
				Munit = M0 / (fluxGoal/F0)
			else:
				Munit = M1 / (fluxGoal/F1)

		#Sometimes this sampler will push you off to crazy values.  Prevent that.
		Munit = np.maximum(smartLimits[0], np.minimum(smartLimits[1], Munit))

		if Munit == M1:
			#The sampler has gotten stuck.
			raise ValueError("The sampler has gotten stuck at an Munit of {0:1.3e}.".format(Munit))

		flux = computeAverageTotalFlux(dumpFiles, Munit, thetacam, freq_Hz=freq_Hz, MBH=MBH, \
		npixel_max=npixel_max, npixel_min=npixel_min, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, fov=fov, rmax_geo=rmax_geo, \
		counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, model_kappa=model_kappa, emission_type=emission_type, \
		Rhigh=Rhigh, adaptiveRefinement=adaptiveRefinement)
		Munit_trial_list.append(Munit)
		flux_trial_list.append(flux)
		if flux < fluxGoal:
			M0 = Munit
			F0 = flux
		else:
			M1 = Munit
			F1 = flux
		if (flux >= fluxGoal / (1.0 + fractionalTolerance)) & (flux <= fluxGoal * (1.0 + fractionalTolerance)):
			solutionFound = True

	return Munit_trial_list, flux_trial_list

def addToMunitTable(folder, thetacam, Munit_table_name, fluxGoal=2.0, fractionalTolerance=0.02, Munit_guess=1e19, smartLimits=[1e10,1e40], \
	dumpSamples=100, dumpRange=[1000,2000], freq_Hz=230e9, MBH=4.14e6, npixel_max=129, npixel_min=17, ipoleExecutable='../ipole_critical_beta/ipole', dsource=8.127e3, \
	phicam=0.0, fov=200.0, rmax_geo=50, counterjet=0, target_nturns=-1, beta_crit_coefficient=0.5, beta_crit=1, keepImage=False, \
	overwrite=False, model_kappa=3.5, emission_type=4, Rhigh=20.0, formatting='Rhigh', adaptiveRefinement=True, improveGuess=True, \
	improvedGuessFractionalTolerance=0.5):

	"""
	Compute an Munit for a given inclination and folder of dumpfiles.  Add to the table.
	"""

	#Get a list of files.
	if folder[-1] != '/':
		folder += '/'

	#Compatible with two different naming schemes.
	availableDumpFiles = np.array([file for file in os.listdir(folder) if (file[:4] == 'dump') | (file[:5] == 'torus')])
	availableDumpNumbers = np.zeros(len(availableDumpFiles), dtype=int)
	for name_index in range(len(availableDumpFiles)):
		if 'dump' in availableDumpFiles[name_index]:
			availableDumpNumbers[name_index] = int(availableDumpFiles[name_index].split('_')[1].split('.')[0])
		elif 'torus' in file:
			availableDumpNumbers[name_index] = int(availableDumpFiles[name_index].split('.')[2])
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

	#To save time, make sure that the Munit_guess is close for a single snapshot.  Taking the one in the middle.
	if improveGuess:
		print("Refining the guessed Munit by fitting a single snapshot with a fractional tolerance of {0:1.2f}.".format(improvedGuessFractionalTolerance))
		singleDumpFile = [usedDumpFiles[int(np.floor(len(usedDumpFiles)/2))]]
		Munit_trials_for_improved_guess, _ = optimizeMunit(singleDumpFile, thetacam, fluxGoal, fractionalTolerance=improvedGuessFractionalTolerance, Munit_guess=Munit_guess, smartLimits=smartLimits, \
		freq_Hz=freq_Hz, MBH=MBH, npixel_max=npixel_max, npixel_min=npixel_min, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, \
		fov=fov, rmax_geo=rmax_geo, counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, \
		beta_crit=beta_crit, model_kappa=model_kappa, emission_type=emission_type, Rhigh=Rhigh, adaptiveRefinement=adaptiveRefinement)
		Munit_guess = Munit_trials_for_improved_guess[-1]

	#The actual computation, where Munit is fit.
	Munit_trials, flux_trials = optimizeMunit(usedDumpFiles, thetacam, fluxGoal, fractionalTolerance=fractionalTolerance, Munit_guess=Munit_guess, smartLimits=smartLimits, \
	freq_Hz=freq_Hz, MBH=MBH, npixel_max=npixel_max, npixel_min=npixel_min, ipoleExecutable=ipoleExecutable, dsource=dsource, phicam=phicam, \
	fov=fov, rmax_geo=rmax_geo, counterjet=counterjet, target_nturns=target_nturns, beta_crit_coefficient=beta_crit_coefficient, \
	beta_crit=beta_crit, model_kappa=model_kappa, emission_type=emission_type, Rhigh=Rhigh, adaptiveRefinement=adaptiveRefinement)
	print("Obtained an Munit of {0:1.3e} after {1:d} iterations.".format(Munit_trials[-1], len(Munit_trials)))

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

def createEntireMunitTable(outputName, inclinationList=[50], dumpSamples=100, npixel_max=129, npixel_min=17, fov=200, extraTextAfterKey='/', foldersToDo=None, \
	MBH=4.14e6, dsource=8.127e3, model_kappa=3.5, emission_type=4, ipoleExecutable='../ipole_critical_beta/ipole', flipRetrograde=True, fluxGoal=2.0, \
	formatting='Rhigh', Rhigh_list=[1.,10.,20.,40.,80.,160.], beta_crit_coefficient=0.5, beta_crit=1, adaptiveRefinement=True, freq_Hz=230e9, Munit_guess=None, smartLimits=[1e10,1e40]):

	initialTime = time.time()

	initializeMunitTable(outputName, formatting=formatting)
	if foldersToDo is None:
		folders = folderToDumpRange.keys()
	else:
		folders = foldersToDo

	if formatting == 'Rhigh':
		for folder in folders:
			dumpRange = folderToDumpRange[folder]
			if Munit_guess is None:
				if 'MAD' in folder.split('/'):
					Munit_guess = 1e17
				elif 'SANE' in folder.split('/'):
					Munit_guess = 1e20
			for inclination in inclinationList:
				thetacam = inclination
				if flipRetrograde:
					if 'a-' in folder:
						thetacam = 180.0-inclination
				for Rhigh in Rhigh_list:
					#extraTextAfterKey should be set to '/dumps/' on Illinois, but not on Harvard.
					addToMunitTable(folder + extraTextAfterKey, thetacam, outputName, dumpSamples=dumpSamples, npixel_max=npixel_max, npixel_min=npixel_min, Munit_guess=Munit_guess, fov=fov, \
					dumpRange=dumpRange, MBH=MBH, dsource=dsource, model_kappa=model_kappa, emission_type=emission_type, ipoleExecutable=ipoleExecutable, \
					fluxGoal=fluxGoal, Rhigh=Rhigh, formatting=formatting, adaptiveRefinement=adaptiveRefinement, freq_Hz=freq_Hz, smartLimits=smartLimits)
	elif formatting == 'criticalBeta':
		for folder in folders:
			dumpRange = folderToDumpRange[folder]
			if Munit_guess is None:
				if 'MAD' in folder.split('/'):
					Munit_guess = 1e17
				elif 'SANE' in folder.split('/'):
					Munit_guess = 1e20
			for inclination in inclinationList:
				thetacam = inclination
				if flipRetrograde:
					if 'a-' in folder:
						thetacam = 180.0-inclination
				#extraTextAfterKey should be set to '/dumps/' on Illinois, but not on Harvard.
				addToMunitTable(folder + extraTextAfterKey, thetacam, outputName, dumpSamples=dumpSamples, npixel_max=npixel_max, npixel_min=npixel_min, Munit_guess=Munit_guess, fov=fov, \
				dumpRange=dumpRange, MBH=MBH, dsource=dsource, model_kappa=model_kappa, emission_type=emission_type, ipoleExecutable=ipoleExecutable, \
				fluxGoal=fluxGoal, formatting=formatting, beta_crit_coefficient=beta_crit_coefficient, beta_crit=beta_crit, adaptiveRefinement=adaptiveRefinement, freq_Hz=freq_Hz, smartLimits=smartLimits)

	finalTime = time.time()
	print("Finished in {0:1.3e} hours.".format((finalTime-initialTime)/3600))

if __name__ == '__main__':
	import sys

	#The first command line argument tells you which GRMHD simulation to do.
	modelIndexToDo = int(sys.argv[1])
	foldersToDo = [folderToDumpRange.keys()[modelIndexToDo]]

	#The second command line argument tells you which value of Kappa to do.
	kappaIndexToDo = int(sys.argv[2])
	model_kappa = kappaList[kappaIndexToDo]
	emissionType = emissionType_list[kappaIndexToDo]

	#The third command line argument tells you which source to do.
	sourceIndexToDo = int(sys.argv[3])
	source = sources[sourceIndexToDo]

	#This is universal.
	Rhigh_list = [1.0,10.0,40.0,160.0]
	dumpSamples = 100

	#Assign a name based on kappa and the source.
	if emissionType == 4:
		emissionNamePiece = 'thermal'
	elif emissionType == 2:
		emissionNamePiece = 'kappa{0:1.1f}'.format(model_kappa)
	Munit_table_name = './Munit_tables/Munit_table_' + emissionNamePiece + '_' + source + '_v3.txt'
	print("Processing " + foldersToDo[0] + " for " + Munit_table_name + '.')

	#I have some default Munit guesses from an old kappa = 5 table.  Also better to guess lower, since you want to be in the optically thin regime.
	if source == 'M87':
		if 'MAD' in foldersToDo[0]:
			Munit_guess = 2e25 / 10
		elif 'SANE' in foldersToDo[0]:
			Munit_guess = 4e28 / 10
	elif source == 'SgrA':
		if 'MAD' in foldersToDo[0]:
			Munit_guess = 4e17 / 10
		elif 'SANE' in foldersToDo[0]:
			Munit_guess = 6e20 / 10

	#The actual work.
	createEntireMunitTable(Munit_table_name, inclinationList=inclinationList[source], dumpSamples=dumpSamples, npixel_max=npixel[source], fov=fov[source], foldersToDo=foldersToDo, npixel_min=17, \
	formatting='Rhigh', ipoleExecutable='../ipole_versions/ipole_dev_narayan/ipole', adaptiveRefinement=False, flipRetrograde=flipRetrograde[source], fluxGoal=fluxGoal[source], Rhigh_list=Rhigh_list, \
	MBH=MBH[source], dsource=dsource[source], model_kappa=model_kappa, emission_type=emissionType, Munit_guess=Munit_guess)
