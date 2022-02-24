#############################################################################
### Recursively moves through directories and extracts parameters and the ###
###               mean firing rates of neuronal subpopulations            ###
#############################################################################

### Copyright 2021-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
from pathlib import Path
from subprocess import call
from shutil import copyfile
import os
import traceback
from utilityFunctions import readParams

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays
Nl = 50 # the number of excitatory neurons in one line of a quadratic grid
Ne = 2500
Ni = 625
N = 3125 # the number of neurons in the whole network

# extractRecursion
# Recursively looks for data directories and extracts parameters and the Q measure from them
# directory: the directory to look in
# fout: file handle to output file
def extractRecursion(directory, fout):

	data_found = False # specifies if any data has been found
	rawpaths = Path(directory)
	MI = []

	print("Contents of directory " + directory)
	print([str(x) for x in rawpaths.iterdir()])
	rawpaths = Path(directory)

	for x in sorted(rawpaths.iterdir()):

		dest_file = ""

		if x.is_dir():

			full_path = str(x)
			path_tail = os.path.split(full_path)[1] # extract the folder name

			if hasTimestamp(path):

				data_found = True

				if "_TRIPLET" in path_tail:
					[timestamp, prev_desc] = path_tail.split("_TRIPLET", 1)
				else:
					[timestamp, desc] = path_tail.split(" ", 1)

				print("========================")
				print(timestamp + " in " + directory)
				print("------------------------")
				params = readParams(full_path + os.sep + timestamp + "_PARAMS.txt")

				core = np.arange(params[12]) # define the cell assembly

				# write the parameter values to the output file
				fout.write(timestamp + "\t\t")

				for i in range(len(params)):
					fout.write(str(params[i]) + "\t\t")

				# read out *_spike_raster.txt:
				frpath = full_path + os.sep + timestamp + "_spike_raster.txt"
				try:
					f = open(frpath)

				except IOError:
					print('Error opening "' + frpath + '"')
					exit()

				# read activities from file and determine mean activities for different regions
				tws = np.column_stack((np.arange(7,14.5,0.5)-0.25, np.arange(7,14.5,0.5)+0.25))
				CA_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the cell assembly
				c_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the exc. control population
				inh_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the exc. population
				tot_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the whole network

				for line in f:
					segs = line.split('\t\t')

					if (segs[0] != ""):
						t = float(segs[0])

						for j in range(len(tws)):
							if t >= tws[j][0] and t < tws[j][1]: # time window for firing rate
								#print(segs[1])
								n = int(segs[1])

								tot_spikes[j] += 1

								if n < Na:
									CA_spikes[j] += 1
								elif n < Ne:
									c_spikes[j] += 1
								elif n < N:
									inh_spikes[j] += 1
								else:
									print('Error reading from "' + frpath + '"')

								break # exit the for-loop
				f.close()

				for j in range(len(tws)):
					#fout.write(str(CA_spikes[j] / Na / 0.5) + "\t\t" + str(c_spikes[j] / (Ne-Na) / 0.5) + "\t\t" + str(inh_spikes[j] / Ni / 0.5) + \
					#	   "\t\t" + str(tot_spikes[j] / N / 0.5) + "\n")
					fout.write(str(CA_spikes[j] / Na / 0.5) + "\t\t" + str(c_spikes[j] / (Ne-Na) / 0.5) + "\t\t")

				fout.write("\n")

			else:
				ret = extractRecursion(directory + os.sep + path, fout)
				data_found = data_found or ret

	return data_found

try:
	fout = open("firing_rates.txt", "w")

except IOError:
	print('Error opening "firing_rates.txt"')
	exit()

if extractRecursion('.', fout):
	print("========================")

fout.close()
