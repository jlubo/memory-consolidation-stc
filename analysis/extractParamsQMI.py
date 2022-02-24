#############################################################################
### Recursively moves through directories and extracts parameters and the ###
###                 measures Q and MI from simulation data                ###
#############################################################################

### Copyright 2018-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import os
import traceback
from utilityFunctions import *
from calculateQ import *
from calculateMIa import *
from pathlib import Path
from shutil import copyfile

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays

Nl_exc = 40 # the number of excitatory neurons in one line of a quadratic grid
ref_time = "11.0" # readout time for the reference firing rate or weight distribution (typically during learning)

# extractRecursion
# Recursively looks for data directories and extracts parameters and the Q and MI measures from them
# For 8h-recall simulation data that itself does not contain the data from during learing, the related data
# from during learning has to be available, with the same description and an earlier timestamp!
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

		dest_file = ""

		if x.is_dir():

			full_path = str(x)
			path_tail = os.path.split(full_path)[1] # extract the folder name

			if hasTimestamp(path_tail): # simulation data found

				data_found = True
				print("========================")

				if "_TRIPLET" in path_tail: # presumed simulation that contains the learning process and 10s-recall
					[timestamp, prev_desc] = path_tail.split("_TRIPLET", 1)
					prev_full_path = full_path # to copy file containing the reference data to another data folder
					prev_path = path_tail
					prev_timestamp = timestamp
					prev_desc = prev_desc.split(" ", 1)[1] # extract description, to later link with a data folder for 8h-recall
				else: # presumed simulation that contains 8h-recall
					[timestamp, desc] = path_tail.split(" ", 1) # extract timestamp and description, to link with a previous data folder for learning
					if desc == prev_desc: # temporarily copy reference data file if the state of this simulation has been loaded from previous one
						dest_file = full_path + os.sep + "network_plots" + os.sep + timestamp + "_net_" + ref_time + ".txt"
						copyfile(prev_full_path + os.sep + "network_plots" + os.sep + prev_timestamp + "_net_" + ref_time + ".txt", \
						         dest_file)
					else:
						print("No previous data folder found, errors may occur...")

				print(timestamp + " in " + directory)
				print("------------------------")
				params = readParams(full_path + os.sep + timestamp + "_PARAMS.txt")

				core = np.arange(params[12]) # define the cell assembly

				print(params[9])
				readout_time = str(float(params[9].split("at")[1]) + 0.1) # define the readout time for recall
				print(readout_time)

				# write the parameter values to the output file
				fout.write(timestamp + "\t\t")

				for i in range(len(params)):
					fout.write(str(params[i]) + "\t\t")

				fout.write(readout_time + "\t\t")

				(Q, Q_err, v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err) = (np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)

				# compute the Q value
				try:
					Q, Q_err, v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err = calculateQ(full_path + os.sep + "network_plots" + os.sep, timestamp, core, Nl_exc, readout_time, params[16])
				except OSError as e:
					print(traceback.format_exc())
				except ValueError as e:
					print(traceback.format_exc())
				else:
					print("v_as = " + str(v_as) + " +- " + str(v_as_err))
					print("v_ans = " + str(v_ans) + " +- " + str(v_ans_err))
					print("v_ctrl = " + str(v_ctrl) + " +- " + str(v_ctrl_err))
					print("Q = " + str(Q) + " +- " + str(Q_err))

				# compute the MI and selfMI value
				try:
					if "STIP" in params[8]: # stipulated CA
						MI, selfMI = calculateMIa(full_path + os.sep + "network_plots" + os.sep, timestamp, Nl_exc, readout_time, "", core)
					else: # learned CA
						MI, selfMI = calculateMIa(full_path + os.sep + "network_plots" + os.sep, timestamp, Nl_exc, readout_time, ref_time)
				except OSError as e:
					print(traceback.format_exc())
				except ValueError as e:
					print(traceback.format_exc())
				else:
					print("Norm. MIa = " + str(MI / selfMI)) # MI normalized by the self-information of the reference distribution

				# write Q and MI value to the output file
				fout.write(str(v_as) + "\t\t" + str(v_as_err) + "\t\t" + str(v_ans) + "\t\t" + str(v_ans_err) + "\t\t" + str(v_ctrl) + "\t\t" + str(v_ctrl_err) + "\t\t")
				fout.write(str(Q) + "\t\t" + str(Q_err))

				fout.write("\t\t" + str(MI) + "\t\t" + str(selfMI))
				fout.write("\n")

				if dest_file != "":
					os.remove(dest_file) # remove temporary reference data file

			else: # recurse into the next directory
				ret = extractRecursion(directory + os.sep + path_tail, fout)
				data_found = data_found or ret

	return data_found

# main

try:
	fout = open("Params_Q_MI.txt", "a")

except IOError:
	print('Error opening "Params_Q_MI.txt"')
	exit()

if extractRecursion('.', fout):
	print("========================")

fout.close()
