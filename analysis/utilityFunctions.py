#####################################################################################################
### Utility functions for different purposes, mainly regarding the processing of simulation data. ###
#####################################################################################################

### Copyright 2017-2023 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
import os
import re
from pathlib import Path
epsilon = 1e-11 # very small number that is counted as zero

######################################
# readWeightMatrixData
# Reads complete early- and late-phase weight matrix and firing rates from network simulation data
# (also cf. loadWeightMatrix in adjacencyFunctions.py and readWeightMatrix in the simulation code)
# filename: name of the file to read the data from
# Nl_exc: number of neurons in one row/column
# return: the adjacency matrix, the early-phase weight matrix, the late-phase weight matrix, the firing rate vector
def readWeightMatrixData(filename, Nl_exc):

	# read weight matrices and firing rates from file
	try:
		with open(filename) as f:
			rawdata = f.read()
	except OSError:
		raise

	rawdata = rawdata.split('\n\n')
	rawmatrix_h = rawdata[0].split('\n')
	rawmatrix_z = rawdata[1].split('\n')
	rawmatrix_v = rawdata[2].split('\n')

	rows = len(rawmatrix_v)

	if (rows != len(rawmatrix_v[0].split('\t\t'))) or (rows != Nl_exc):
		raise ValueError(str(rows) + ' instead of ' + str(Nl_exc) + ' lines in data file "' + filename + '"')
		f.close()
		exit()

	v = np.zeros((Nl_exc,Nl_exc))
	h = np.zeros((Nl_exc**2,Nl_exc**2))
	z = np.zeros((Nl_exc**2,Nl_exc**2))

	for i in range(Nl_exc**2):
		if i < Nl_exc:
			value0 = rawmatrix_v[i].split('\t\t')
		value1 = rawmatrix_h[i].split('\t\t')
		value2 = rawmatrix_z[i].split('\t\t')

		for j in range(Nl_exc**2):
			if i < Nl_exc and j < Nl_exc:
				v[i][j] = float(value0[j])
			h[i][j] = float(value1[j])
			z[i][j] = float(value2[j])

	f.close()
	connections = (h > epsilon)

	return connections, h, z, v

######################################
# readWeightVectorData
# Reads complete weight vector data from a file
# filename: name of the file to read the data from
# N: the number of synapses
# return: the early-phase weight vector and the late-phase weight vector
def readWeightVectorData(filename, N):

	try:
		with open(filename) as f:
			rawdata = f.read()
	except OSError:
		raise

	rawdata = rawdata.split('\n')

	if (len(rawdata) - 1 != N):
		raise ValueError(str(len(rawdata) - 1) + ' instead of ' + str(N) + ' lines in data file "' + filename + '"')
		f.close()
		exit()

	h = np.zeros(N)
	z = np.zeros(N)

	for i in range(N):
		values = rawdata[i].split('\t\t')
		h[i] = float(values[0])
		z[i] = float(values[1])

	f.close()

	return h, z

######################################
# ftos
# Converts a float value to a string with pre-defined number of decimal places
# f: a floating-point number
# output_places [optional]: number of decimal places for output numbers
# return: a string
def ftos(f, output_places = 7):
	return str(round(f, output_places))

######################################
# cond_print
# Prints an argument only if another argument is true
# pr: printing if True, not printing if False
# arg: the argument to be printed
def cond_print(pr, *args):
	if pr:
		print(*args)

######################################
# extractRandomPattern
# Extracts from a provided log file the list of neurons of the randomly drawn pattern
# - path: path to the log file
# - return: list of neurons forming a random pattern (or None if no random pattern has been found)
def extractRandomPattern(path):
	neuron_list = None
	with open(path, 'r') as f:
		for line in f:
			if line.startswith("Randomly drawn neurons for stimulation:"):
				# extract content between the curly braces using regular expression
				match = re.search(r'{(.*?)}', line)
				if match:
					# split the content by commas and strip spaces to get a list
					neuron_list = [int(neuron.strip()) for neuron in match.group(1).split(',') if neuron.strip().isdigit()]
				break
	return neuron_list

######################################
# readParams
# Reads some parameters from a "[timestamp]_PARAMS.txt" file
# path: path to the parameter file to read the data from
# return: a list containing the values of all the parameters read
def readParams(path):

	try:
		f = open(path)

	except IOError:
		print('Error opening "' + path + '"')
		exit()

	r_CA = -1 # radius of the stimulated core
	s_CA = -1

	# read activities from file and determine mean activities for different regions
	rawdata = f.read()
	rawdata = rawdata.split('\n')
	nn = len(rawdata)
	f.close()

	for i in range(nn):
		segs = rawdata[i].split(' ')

		if segs[0] == "Ca_pre":
			Ca_pre = float(segs[2])
		elif segs[0] == "Ca_post":
			Ca_post = float(segs[2])
		elif segs[0] == "theta_p":
			theta_p = float(segs[2])
		elif segs[0] == "theta_d":
			theta_d = float(segs[2])
		elif segs[0] == "R_mem":
			R_mem = float(segs[2])
		elif segs[0] == "learning":
			if segs[3] == "":
				lprot = "none"
			else:
				lprot = segs[3]
		elif segs[0] == "stimulus": # old version of previous condition
			lprot = segs[2]
		elif segs[0] == "recall" and segs[1] == "stimulus":
			rprot = segs[3]
		elif segs[0] == "recall" and segs[1] != "fraction": # old version of previous condition
			rprot = segs[2]
		elif segs[0] == "recall" and segs[1] == "fraction":
			recall_fraction = segs[3]
		elif segs[0] == "pc" or segs[0] == "p_c":
			p_c = float(segs[2])
		elif segs[0] == "w_ei":
			w_ei = float(segs[2])
		elif segs[0] == "w_ie":
			w_ie = float(segs[2])
		elif segs[0] == "w_ii":
			w_ii = float(segs[2])
		elif segs[0] == "theta_pro_c":
			theta_pro_c = float(segs[2])
		elif segs[0] == "N_stim":
			N_stim = int(segs[2])
		elif segs[0] == "I_const" or segs[0] == "I_0":
			I_0 = float(segs[2])
		elif segs[0] == "dt":
			dt = float(segs[2])
		elif segs[0] == "core" and (segs[2] == "first" or segs[2] == "random"):
			s_CA = int(segs[3])

	if s_CA == -1:
		print('Warning: the size of the stimulated core could not be determined.')
		s_CA = int(input('Enter it now to continue: '))

	return [w_ei, w_ie, w_ii, p_c, Ca_pre, Ca_post, theta_p, theta_d, lprot, rprot, dt, theta_pro_c, s_CA, N_stim, I_0, R_mem, recall_fraction]

######################################
# hasTimestamp
# Checks if the name of the given file starts with a 17-character timestamp of the format XX-XX-XX_XX-XX-XX
# path: path pointing to the file
# return: True if probably there is a timestamp, False if not
def hasTimestamp(path):

	try:
		offset = 0 # number of characters after which the timestamp begins
		filename = os.path.split(path)[1] # only consider the tail of the path
		for prefix in ["data_", "tmp_data_"]:
			if filename[0:len(prefix)] == prefix:
				offset = len(prefix)

		if filename[offset+2] == "-" and filename[offset+5] == "-" and filename[offset+8] == "_" and \
		   filename[offset+11] == "-" and filename[offset+14] == "-":
			return True
	except:
		pass

	return False

######################################
# getTimestamp
# Gets the 17-character timestamp from a given filename
# filename: a string
# return: the timestamp as a string, None if an exception occurs
def getTimestamp(filename):

	try:
		offset = 0 # number of characters after which the timestamp begins
		for prefix in ["data_", "tmp_data_"]:
			if filename[0:len(prefix)] == prefix:
				offset = len(prefix)

		return filename[offset:offset+17] 
	except:
		pass

	return None

######################################
# mergeRawData
# To merge data from multiple files into a single file - looks in a specified directory
# for files with a certain substring in the filename, and merges them (merging the content of the lines)
# rootpath: relative path to the output directory
# substr: string that the filename of files to be merged has to contain
# output_file: name of the output file
# remove_raw [optional]: removes the raw data files
# sep_str [optional]: the character or string by which to separate the lines in the output file
def mergeRawData(rootpath, substr, output_file, remove_raw=False, sep_str='\t\t'):

	path = Path(rootpath)
	num_rows = -1
	all_data = []

	for x in sorted(path.iterdir()): # loop through files in the output directory
		x_str = str(x)
		if not x.is_dir() and substr in x_str:

			f = open(x_str)
			single_trial_data = f.read()
			f.close()

			single_trial_data = single_trial_data.split('\n')

			if single_trial_data[-1] == "":
				del single_trial_data[-1] # delete empty line

			if len(single_trial_data) != num_rows:
				if num_rows == -1:
					num_rows = len(single_trial_data)
					all_data = single_trial_data
				else:
					raise Exception("Wrong number of rows encountered in: " + x_str)
			else:
				for i in range(num_rows):
					all_data[i] += sep_str + single_trial_data[i]

			if remove_raw:
				os.remove(x_str)

	fout = open(os.path.join(rootpath, output_file), "w")
	for i in range(num_rows):
		fout.write(all_data[i] + '\n')
	fout.close()

######################################
# countNeuronData
# Counts the datasets of neuron observables in a data file
# filename: name of the file from which to read
# return: number of the datasets or -1 if failed
def countNeuronData(filename):
	
	f = open(filename, 'r')
	line = f.readline() # read header line

	if not line: # end of file is reached
		return -1
	f.close()

	num_mem = len(re.findall(r'[ \t]V\(', line))
	num_curr = len(re.findall(r'[ \t]I\_tot\(', line))
	if num_curr == 0:
		num_curr = len(re.findall(r'[ \t]I\(', line))
	num_p = len(re.findall(r'[ \t]p\^C\(', line))

	#print("num_mem =", num_mem, "num_curr=", num_curr, "num_p =", num_p)

	if num_mem == num_curr == num_p: # check for consistency
		return num_mem
	else:
		return -1

######################################
# countSynapseData
# Counts the datasets of synapse observables in a data file
# filename: name of the file from which to read
# return: number of the datasets or -1 if failed
def countSynapseData(filename):
	
	f = open(filename, 'r')
	line = f.readline() # read header line

	if not line: # end of file is reached
		return -1
	f.close()

	num_h = line.count("h(")
	num_z = line.count("z(")
	num_Ca = line.count("Ca(")

	if num_h == num_z == num_Ca: # check for consistency
		return num_h
	else:
		return -1
