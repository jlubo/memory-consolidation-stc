#############################################################################
### Recursively moves through directories and extracts parameters and the ###
###         mean protein amount across core and non-core neurons          ###
#############################################################################

### Copyright 2021-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import os
from utilityFunctions import *
from pathlib import Path
from shutil import copyfile

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays
Nl = 40 # the number of excitatory neurons in one line of a quadratic grid
t_readout = "14.5" # "3621" # 7221 # time at which to read out the protein concentrations

# extractRecursion
# Recursively looks for data directories and extracts parameters and protein concentrations at a specified time
# (concomitantly also extracts mean weight within the core and within the non-core population at the same time)
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
					data_file = os.path.join(full_path, timestamp + "_mean_weight.txt")
				else:
					print("Error: unknown stimulation protocol.")
					exit()

				print(timestamp + " in " + directory)
				print("------------------------")
				params = readParams(full_path + os.sep + timestamp + "_PARAMS.txt")

				# read from data file the protein concentrations at time t_readout
				fin = open(data_file)
				all_data = fin.read()
				fin.close()

				all_data = all_data.split("\n")

				extracted_data = None
				for el in all_data:
					if el.find(t_readout + "\t") == 0: # if the row starts with the readout time
						if extracted_data is None:
							extracted_data = el
						else:
							print("Warning: more than one matching data row found in " + data_file)

				if extracted_data is None:
					print("Error: no matching data row found in " + data_file)
					exit()

				# write the parameter values to the output file
				fout.write(timestamp + "\t\t")

				for i in range(len(params)):
					fout.write(str(params[i]) + "\t\t")

				# write extracted result to the output file
				fout.write(str(extracted_data) + "\n")

			else: # recurse into the next directory
				ret = extractRecursion(directory + os.sep + path_tail, fout)
				data_found = data_found or ret

	return data_found

# main

try:
	fout = open("Params_ProteinAmount.txt", "a")

except IOError:
	print('Error opening "Params_ProteinAmount.txt"')
	exit()

if extractRecursion('.', fout):
	print("========================")

fout.close()
