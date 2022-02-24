#############################################################################
### Recursively moves through directories and extracts parameters and the ###
###                mean weights across neuronal subpopulations            ###
#############################################################################

### Copyright 2021-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import os
import traceback
from utilityFunctions import *
import adjacencyFunctions as adj
from pathlib import Path
from shutil import copyfile

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays
Nl = 40 # the number of excitatory neurons in one line of a quadratic grid

# extractRecursion
# Recursively looks for data directories and extracts parameters and mean weights and weight standard deviations from
# neuronal subpopulations
# directory: the directory to look in
# fout: file handle to output file
def extractRecursion(directory, fout):

	data_found = False # specifies if any data has been found
	rawpaths = Path(directory)
	MI = []

	print("Contents of directory " + directory)
	print([str(x) for x in rawpaths.iterdir()])
	rawpaths = Path(directory)

	for x in sorted(rawpaths.iterdir()): # loop through elements in this directory

		if x.is_dir():

			full_path = str(x)
			path_tail = os.path.split(full_path)[1] # extract the folder name

			if hasTimestamp(path_tail): # simulation data found

				data_found = True
				print("========================")

				if "_TRIPLET" in path_tail: # presumed simulation that contains the learning process and 10s-recall
					[timestamp, desc] = path_tail.split("_TRIPLET", 1)
					weight_matrix_file = os.path.join(full_path, "network_plots", timestamp + "_net_11.0.txt")

				else: # presumed simulation that contains 8h-recall
					[timestamp, desc] = path_tail.split(" ", 1) # extract timestamp and description, to link with a previous data folder for learning
					weight_matrix_file = os.path.join(full_path, "network_plots", timestamp + "_net_28810.1.txt")

				print(timestamp + " in " + directory)
				print("------------------------")
				params = readParams(full_path + os.sep + timestamp + "_PARAMS.txt")
				core = np.arange(params[12]) # define the cell assembly

				# write the parameter values to the output file
				fout.write(timestamp + "\t\t")

				for i in range(len(params)):
					fout.write(str(params[i]) + "\t\t")

				# compute the mean weights and standard deviations and write them to the output file
				all = np.arange(Nl**2)
				noncore = all[np.logical_not(np.in1d(all, core))] # neurons that are not in the cell assembly core

				adj.loadWeightMatrix(weight_matrix_file)

				print("Core -> core ('CA'):")
				hm = adj.meanEarlyPhaseWeight(core, pr=False)
				hsd = adj.sdEarlyPhaseWeight(core, pr=False)
				zm = adj.meanLatePhaseWeight(core, pr=False)
				zsd = adj.sdLatePhaseWeight(core, pr=False)
				fout.write(ftos(hm) + "\t\t" + ftos(hsd) + "\t\t" + ftos(zm) + "\t\t" + ftos(zsd) + "\t\t")
				print(" Mean total weight: " + ftos(hm + zm) + " +- " + ftos(np.sqrt(hsd**2 + zsd**2)))

				print("Core -> non-core ('outgoing'):")
				hm = adj.meanEarlyPhaseWeight(core, noncore, pr=False)
				hsd = adj.sdEarlyPhaseWeight(core, noncore, pr=False)
				zm = adj.meanLatePhaseWeight(core, noncore, pr=False)
				zsd = adj.sdLatePhaseWeight(core, noncore, pr=False)
				fout.write(ftos(hm) + "\t\t" + ftos(hsd) + "\t\t" + ftos(zm) + "\t\t" + ftos(zsd) + "\t\t")
				print(" Mean total weight: " + ftos(hm + zm) + " +- " + ftos(np.sqrt(hsd**2 + zsd**2)))

				print("Non-core -> core ('incoming'):")
				hm = adj.meanEarlyPhaseWeight(noncore, core, pr=False)
				hsd = adj.sdEarlyPhaseWeight(noncore, core, pr=False)
				zm = adj.meanLatePhaseWeight(noncore, core, pr=False)
				zsd = adj.sdLatePhaseWeight(noncore, core, pr=False)
				fout.write(ftos(hm) + "\t\t" + ftos(hsd) + "\t\t" + ftos(zm) + "\t\t" + ftos(zsd) + "\t\t")
				print(" Mean total weight: " + ftos(hm + zm) + " +- " + ftos(np.sqrt(hsd**2 + zsd**2)))

				print("Non-core -> non-core ('non-CA'):")
				hm = adj.meanEarlyPhaseWeight(noncore, pr=False)
				hsd = adj.sdEarlyPhaseWeight(noncore, pr=False)
				zm = adj.meanLatePhaseWeight(noncore, pr=False)
				zsd = adj.sdLatePhaseWeight(noncore, pr=False)
				fout.write(ftos(hm) + "\t\t" + ftos(hsd) + "\t\t" + ftos(zm) + "\t\t" + ftos(zsd) + "\n")
				print(" Mean total weight: " + ftos(hm + zm) + " +- " + ftos(np.sqrt(hsd**2 + zsd**2)))

			else: # recurse into the next directory
				ret = extractRecursion(directory + os.sep + path_tail, fout)
				data_found = data_found or ret

	return data_found

# main

try:
	fout = open("Params_MeanWeight.txt", "a")

except IOError:
	print('Error opening "Params_MeanWeight.txt"')
	exit()

if extractRecursion('.', fout):
	print("========================")

fout.close()
