#################################################################################
### Script to analyze the weight structure in a network (weight distributions ###
###              and mean weight within and between subpopulations)           ###
#################################################################################

### Copyright 2020-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

### example call from shell: python3 analyzeWeights.py "Weight Distributions and Mean Weight Matrix" "OVERLAP10 no AC, no ABC"

import valueDistributions as vd
import adjacencyFunctions as adj
from overlapParadigms import *
from utilityFunctions import *
import sys

##############################################################################################
### initialize
N_pop = 2500 # number of neurons in the considered population
core_size = 600 # number of excitatory neurons in one cell assembly
MWM = False # specifies whether to create abstract mean weight matrix
MCW = False # specifies whether to create file with mean core weights
WD = False # specifies whether to create weight distribution plots

if len(sys.argv) < 3: # if there are less than 2 commandline arguments
	print("No argument provided! Running the default routine. Paradigm may be wrong.")
	WD = True
	MWM = True
	paradigm = "NOOVERLAP"
else:
	# read strings from argument 2, telling what analyses to perform
	if "Weight Distributions" in str(sys.argv[1]):
		WD = True
	if "Mean Weight Matrix" in str(sys.argv[1]):
		MWM = True
	if "Mean Core Weights" in str(sys.argv[1]):
		MCW = True

	# read argument 3, telling what paradigm to consider
	paradigm = str(sys.argv[2])

print("Paradigm:", paradigm)

##############################################################################################
### get cell assembly definitions
try:
	coreA, coreB, coreC = coreDefinitions(paradigm, core_size)
except:
	raise

##############################################################################################
### look for network output files in this directory
rawpaths = Path(".")

for x in sorted(rawpaths.iterdir()):

	full_path = str(x)
	hpath = os.path.split(full_path)[0] # take head
	tpath = os.path.split(full_path)[1] # take tail

	if not x.is_dir():

		if hasTimestamp(tpath) and "_net_" in tpath and not "_av_" in tpath and not ".png" in tpath:
			timestamp = tpath.split("_net_")[0]
			time_for_readout = tpath.split("_net_")[1].split(".txt")[0]

			if WD:
				print("Plotting weight distributions from dataset", timestamp, "with time", time_for_readout)
				N_pop_row = int(round(np.sqrt(N_pop)))
				vd.plotWeightDistributions3CAs(".", timestamp, "", N_pop_row, time_for_readout, coreA, coreB, coreC)

			if MWM:
				print("Creating abstract mean weight matrix from dataset", timestamp, "with time", time_for_readout)
				adj.meanWeightMatrix(timestamp, time_for_readout, coreA, coreB, coreC, N_pop, pr = True)

			if MCW:
				print("Computing mean core weights from dataset", timestamp, "with time", time_for_readout)
				adj.meanCoreWeights(timestamp, time_for_readout, coreA, coreB, coreC, N_pop, pr = True)
