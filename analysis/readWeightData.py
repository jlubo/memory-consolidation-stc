################################################################################################
### Functions to read the connections, early- and late-phase weight matrix, and firing rates ###
###                               from network simulation data                               ###
################################################################################################

### Copyright 2017-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np

# readWeightMatrixData
# Reads complete weight matrix data from a file (modified from plotFunctions.py)
# filename: name of the file to read the data from
# Nl: number of neurons in one row/column
# return: the adjacency matrix, the early-phase weight matrix, the late-phase weight matrix, the firing rate vector
def readWeightMatrixData(filename, Nl):

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

	if (rows != len(rawmatrix_v[0].split('\t\t'))) or (rows != Nl):
		raise ValueError(str(rows) + ' instead of ' + str(Nl) + ' lines in data file "' + filename + '"')
		f.close()
		exit()

	v = np.zeros((Nl,Nl))
	h = np.zeros((Nl**2,Nl**2))
	z = np.zeros((Nl**2,Nl**2))

	for i in range(Nl**2):
		if i < Nl:
			value0 = rawmatrix_v[i].split('\t\t')
		value1 = rawmatrix_h[i].split('\t\t')
		value2 = rawmatrix_z[i].split('\t\t')

		for j in range(Nl**2):
			if i < Nl and j < Nl:
				v[i][j] = float(value0[j])
			h[i][j] = float(value1[j])
			z[i][j] = float(value2[j])

	f.close()
	connections = (h > 0)

	return connections, h, z, v

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
