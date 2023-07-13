######################################################################################
### Recursively moves through directories to extract parameters, firing rates, and ###
###          the recall quality measures Q and MI from spike raster data           ###
######################################################################################

### Copyright 2018-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import os
import traceback
from utilityFunctions import hasTimestamp, readParams
from calculateQ import getSubpopulations, calculateQ
from calculateMIa import calculateMIa
from pathlib import Path
from shutil import copyfile
import json

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays
Ne = 1600 # the number of excitatory neurons 
Ni = 400 # the number of inhibitory neurons 
N = Ne + Ni # the number of neurons in the whole network
time_window_size = 0.5 # duration of the time window for firing rate computation (in seconds)

# extractRecursion
# Recursively looks for directories with spike raster data and extracts parameters, firing rates, and the recall performance measures Q and MI
# sup_directory: the directory to look in
# fout: file handle to output file
# col_sep [optional]: characters separating columns in the data file
def extractRecursion(sup_directory, fout, col_sep = '\t\t'):

	data_found = False # specifies if any data has been found
	rawpaths = Path(sup_directory)
	MI = []

	print("Contents of directory " + sup_directory)
	print([str(x) for x in rawpaths.iterdir()])
	rawpaths = Path(sup_directory)

	for x in sorted(rawpaths.iterdir()): # loop through elements in this directory

		dest_file = ""

		full_path = str(x)
		path_tail = os.path.split(full_path)[1] # extract the filename

		if hasTimestamp(path_tail) \
		   and ("_spikes.txt" in path_tail or "_spike_raster.txt" in path_tail):

			data_found = True
			if "_spikes.txt" in path_tail:
				other_simulator = True
				col_sep = ' '
			else:
				other_simulator = False
				col_sep = '\t\t'

			[timestamp, _] = path_tail.split("_spike", 1)

			print("========================")
			print(timestamp + " in " + sup_directory)
			print("------------------------")
			
			# setting the parameters
			if other_simulator:
				# load parameter configuration from JSON file
				config = json.load(open(os.path.join(sup_directory, timestamp + "_config.json"), "r"))
				Na = config['populations']['N_CA']
				p_r = config['populations']['p_r']

				# define the readout time for learning
				readout_time_0 = float(config['simulation']['learn_protocol']['time_start']) + 1.0
				print("readout_time_0 =",  readout_time_0)

				# define the readout time for recall
				readout_time_1 = float(config['simulation']['recall_protocol']['time_start']) + 0.1
				print("readout_time_1 =",  readout_time_1)
				
				# get values of important parameters
				w_ei = config['populations']['w_ei']
				w_ie = config['populations']['w_ie']
				w_ii = config['populations']['w_ii']

			else:
				# load parameter configuration from *_PARAMS.txt file
				params = readParams(os.path.join(sup_directory, timestamp + "_PARAMS.txt"))
				Na = int(params[12])
				p_r = float(params[16])

				# define the readout time for learning/reference
				if len(params[8].split("at")) > 2:
					readout_time_0 = float(params[8].split("at")[1]) + 1.0 # 1.0 sec. after onset of learning stimulus
				else: 
					readout_time_0 = 11.0 # default value if not specified
				print("readout_time_0 =",  readout_time_0)

				# define the readout time for recall
				readout_time_1 = float(params[9].split("at")[1]) + 0.1 # 0.1 sec. after onset of recall stimulus
				print("readout_time_1 =",  readout_time_1)

				# get values of important parameters
				w_ei = params[0]
				w_ie = params[1]
				w_ii = params[2]

			# define the cell assembly (the subpopulation that receives learning stimulation)
			core = np.arange(Na)

			# define the subpopulations of the network with resepect to recall
			core_as, core_ans, control = getSubpopulations(Ne, core, p_r)
			Nas = len(core_as)
			Nans = len(core_ans)
				
			# write the values of important parameters to the output file
			fout.write(timestamp + col_sep + \
				       str(Na) + col_sep + \
				       str(w_ei) + col_sep + \
				       str(w_ie) + col_sep + \
				       str(w_ii) + col_sep + \
				       str(readout_time_1) + col_sep)

			# read out *_spikes.txt / *_spike_raster.txt:
			try:
				f = open(full_path)

			except IOError:
				print('Error opening "' + full_path + '"')
				exit()

			# read all spikes from file and count spikes for the different subpopulations and for the individual excitatory neurons and
			time_windows = [ [readout_time_0-time_window_size/2, readout_time_0+time_window_size/2], \
			                 [readout_time_1-time_window_size/2, readout_time_1+time_window_size/2] ]
			as_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in assembly core neurons that receive recall stimulation
			ans_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in assembly core neurons that do not receive recall stimulation
			ctrl_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in the exc. control population
			exc_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in the exc. population
			inh_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in the inh. population
			tot_spikes = np.zeros(len(time_windows), dtype=int) # number of spikes in the whole network
			exc_spikes_individual = np.zeros((len(time_windows), Ne), dtype=int) # number of spikes for each neuron

			for line in f:
				segs = line.split(col_sep)

				if (segs[0] != ""):
					if other_simulator:
						t = float(segs[0].lstrip("\x00")) / 1000 # convert ms to s
					else:
						t = float(segs[0].lstrip("\x00"))

					for j in range(len(time_windows)):
						if t >= time_windows[j][0] and t < time_windows[j][1]: # time window for firing rate
							n = int(segs[1])

							tot_spikes[j] += 1

							if n < Ne: # excitatory population
								exc_spikes[j] += 1
								exc_spikes_individual[j][n] += 1
								if n < Nas: # assembly core, stimulated for recall
									as_spikes[j] += 1
								elif n < Na: # assembly core, not stimulated for recall
									ans_spikes[j] += 1
								elif n < Ne: # control
									ctrl_spikes[j] += 1
							elif n < N: # inhibitory population
								inh_spikes[j] += 1
							else:
								print('Error reading from "' + full_path + '"')

							break # exit the for-loop

			f.close()

			# determine activities for the individual excitatory neurons for each readout time (the activity distributions)
			v_exc_individual = exc_spikes_individual / time_window_size

			# determine mean activities (and standard deviations) for the subpopulations at recall (i.e., at 'readout_time_1')
			v_exc = exc_spikes[1] / Ne / time_window_size # mean firing rate of exc. population
			v_as = as_spikes[1] / Nas / time_window_size # mean firing rate of exc. subpopulation that receives learning and recall stimulation
			v_ans = ans_spikes[1] / Nans / time_window_size # mean firing rate of exc. subpopulation that receives learning but no recall stimulation
			v_ctrl = ctrl_spikes[1] / (Ne-Na) / time_window_size # mean firing rate of exc. subpopulation that receives neither learning nor recall stimulation
			v_inh = inh_spikes[1] / Ni / time_window_size # mean firing rate of inh. population

			v_exc_err = np.sqrt(np.sum(np.power(v_exc_individual, 2)) / Ne - np.power(v_exc, 2))
			v_as_err = np.sqrt(np.sum(np.power(v_exc_individual[1][core_as], 2)) / len(core_as) - np.power(v_as, 2))
			v_ans_err = np.sqrt(np.sum(np.power(v_exc_individual[1][core_ans], 2)) / len(core_ans) - np.power(v_ans, 2))
			v_ctrl_err = np.sqrt(np.sum(np.power(v_exc_individual[1][control], 2)) / len(control) - np.power(v_ctrl, 2))
			v_inh_err = np.nan # no distribution computed to compute this value (TODO?)
					
			# compute the Q value from mean activities at recall
			try:
				Q, Q_err = calculateQ(v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err)
			except OSError as e:
				print(traceback.format_exc())
			except ValueError as e:
				print(traceback.format_exc())
			else:
				print("v_exc = " + str(v_exc) + " +- " + str(v_exc_err))
				print("v_as = " + str(v_as) + " +- " + str(v_as_err))
				print("v_ans = " + str(v_ans) + " +- " + str(v_ans_err))
				print("v_ctrl = " + str(v_ctrl) + " +- " + str(v_ctrl_err))
				print("v_inh = " + str(v_inh) + " +- " + str(v_inh_err))
				print("Q = " + str(Q) + " +- " + str(Q_err))

			# compute the MI value (and entropies) from the activity distributions at learning and recall
			try:
				MI, selfMIL, selfMIR = calculateMIa(v_exc_individual[0], v_exc_individual[1])
			except OSError as e:
				print(traceback.format_exc())
			except ValueError as e:
				print(traceback.format_exc())
			else:
				print("Norm. MIa = " + str(MI / selfMIL)) # MI normalized by the self-information of the reference distribution'''

			# write Q and MI value to the output file
			fout.write(str(v_exc) + col_sep + str(v_exc_err) + col_sep + \
			           str(v_as) + col_sep + str(v_as_err) + col_sep + str(v_ans) + col_sep + str(v_ans_err) + col_sep + str(v_ctrl) + col_sep + str(v_ctrl_err) + col_sep + \
			           str(v_inh) + col_sep + str(v_inh_err) + col_sep)
			fout.write(str(Q) + col_sep + str(Q_err) + col_sep)
			fout.write(str(MI) + col_sep + str(selfMIL) + col_sep + str(selfMIR))

			fout.write("\n")

		elif x.is_dir(): # recurse into the next directory
			ret = extractRecursion(sup_directory + os.sep + path_tail, fout)
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
