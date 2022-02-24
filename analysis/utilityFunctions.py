########################################################################################
### Utility function for different purposes regarding the reading of simulation data ###
########################################################################################

### Copyright 2017-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import os
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
	return str(round(f, output_precision))

######################################
# cond_print
# Prints an argument only if another argument is true
# pr: printing if True, not printing if False
# arg: the argument to be printed
def cond_print(pr, *args):
	if pr:
		print(*args)

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
# Checks if the given filename starts with a timestamp
# filename: a string
# return: true if presumably there is a timestamp, false if not
def hasTimestamp(filename):
	try:
		if filename[2] == "-" and filename[5] == "-" and filename[8] == "_" and \
		   filename[11] == "-" and filename[14] == "-":
			return True
	except:
		pass

	return False

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
