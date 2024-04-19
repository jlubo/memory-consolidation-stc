#############################################################################
### Recursively moves through directories and extracts parameters and the ###
###               mean firing rates of neuronal subpopulations            ###
#############################################################################

### Copyright 2021-2023 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
from pathlib import Path
from subprocess import call
from shutil import copyfile
import os
from utilityFunctions import hasTimestamp, getTimestamp, readParams
import json

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays
Ne = 1600 # the number of excitatory neurons 
Ni = 400 # the number of inhibitory neurons 
N = Ne + Ni # the number of neurons in the whole network

# extractRecursion
# Recursively looks for data directories and extracts parameters and the Q measure from them
# sup_directory: the directory to look in
# fout: file handle to output file
# col_sep [optional]: characters separating columns in the data file
# other_simulator [optional]: specifies if alternative simulator (such as Arbor) is being used
def extractRecursion(sup_directory, fout, col_sep = '\t\t', other_simulator = False):

	data_found = False # specifies if any data has been found
	rawpaths = Path(sup_directory)
	MI = []

	print("Contents of directory " + sup_directory)
	print([str(x) for x in rawpaths.iterdir()])
	rawpaths = Path(sup_directory)

	for x in sorted(rawpaths.iterdir()):

		dest_file = ""

		if x.is_dir():

			full_path = str(x)
			path_tail = os.path.split(full_path)[1] # extract the folder name

			if hasTimestamp(path_tail):

				data_found = True

				timestamp = getTimestamp(path_tail)

				print("========================")
				print(timestamp + " in " + sup_directory)
				print("------------------------")

				if other_simulator:
					# load parameter configuration from JSON file
					config = json.load(open(os.path.join(full_path, timestamp + "_config.json"), "r"))
					Na = config['populations']['N_CA'] # core size

					# write the timestamp to the output file
					fout.write(timestamp + col_sep)

				else:
					params = readParams(full_path + os.sep + timestamp + "_PARAMS.txt")
					Na = params[12] # core size

					# write the timestamp to the output file
					fout.write(timestamp + col_sep)

					# write the parameter values to the output file
					for i in range(len(params)):
						fout.write(str(params[i]) + col_sep)
									
				# define the cell assembly
				core = np.arange(Na)

				# read out spike data file
				if other_simulator:
					data_path = os.path.join(full_path, timestamp + "_spikes.txt")
				else:
					data_path = os.path.join(full_path, timestamp + "_spike_raster.txt")
				try:
					f = open(data_path)
				except IOError:
					print('Error opening "' + data_path + '"')
					exit()

				# read activities from file and determine mean activities for different regions
				tws = np.column_stack((np.arange(7,14.5,0.5)-0.25, np.arange(7,14.5,0.5)+0.25))
				CA_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the cell assembly
				c_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the exc. control population
				inh_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the exc. population
				tot_spikes = np.zeros(len(tws), dtype=int) # number of spikes in the whole network

				for line in f:
					segs = line.split(col_sep)

					if (segs[0] != ""):
						if other_simulator:
							t = float(segs[0]) / 1000 # convert ms to s
						else:
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
									print('Error reading from "' + full_path + '"')

								break # exit the for-loop
				f.close()

				for j in range(len(tws)):
					#fout.write(str(CA_spikes[j] / Na / 0.5) + col_sep + str(c_spikes[j] / (Ne-Na) / 0.5) + col_sep + str(inh_spikes[j] / Ni / 0.5) + \
					#	   col_sep + str(tot_spikes[j] / N / 0.5) + "\n")
					fout.write(str(CA_spikes[j] / Na / 0.5) + col_sep + str(c_spikes[j] / (Ne-Na) / 0.5) + col_sep)

				fout.write("\n")

			else:
				ret = extractRecursion(os.path.join(sup_directory, path_tail), fout, col_sep, other_simulator)
				data_found = data_found or ret

	return data_found

try:
	fout = open("firing_rates.txt", "w")

except IOError:
	print('Error opening "firing_rates.txt"')
	exit()

other_simulator = False
if other_simulator:
	col_sep = ' '
else:
	col_sep = '\t\t'
if extractRecursion('.', fout, col_sep, other_simulator):
	print("========================")

fout.close()
