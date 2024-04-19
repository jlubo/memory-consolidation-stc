######################################################
### Averages early- and late-phase weight matrices ###
######################################################

### Copyright 2019-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
from pathlib import Path

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays

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

# averageWeights
# Averages early- and late-phase weight matrices over all trials that are available
# (trial data must be given as *_net_* file in a given directory)
# nppath: path to the directory to read the data from
# Nl: the number of excitatory neurons in one line of a quadratic grid
# time: the time that at which the weights shall be read out
def averageWeights(nppath, Nl, time):


	c = np.zeros((Nl**2,Nl**2))
	h = np.zeros((Nl**2,Nl**2))
	z = np.zeros((Nl**2,Nl**2))
	v = np.zeros((Nl,Nl))

	counter = 0

	# Looking for data files ("*_net_*")
	rawpaths = Path(nppath)

	for x in rawpaths.iterdir():

		path = str(x)

		if ("_net_" + time + ".txt") in path:

			try:
				ctmp, htmp, ztmp, vtmp = readWeightMatrixData(path, Nl)
			except ValueError:
				raise
			except OSError:
				raise
			except:
				print("Error in " + path)
				exit()

			c = c + ctmp
			h = h + htmp
			z = z + ztmp
			v = v + vtmp

			counter += 1


	print("Averaged over " + str(counter) + " trials for t = " + str(time) + " s.")
	c /= counter
	h /= counter
	z /= counter
	v /= counter

	# write averaged _net_ file containing averaged early-/late-phase weights and firing rates
	f = open('net_' + time + '_averaged.txt','wb')
	np.savetxt(f, h, fmt='%.6f', delimiter='\t\t')
	f.write(b'\x0a')
	np.savetxt(f, z, fmt='%.6f', delimiter='\t\t')
	f.write(b'\x0a')
	np.savetxt(f, v, fmt='%.6f', delimiter='\t\t')
	f.write(b'\x0a\x0a')
	f.close()

	# write averaged connectivity matrix (just for sanity check)
	#f = open('conn_' + time + '_averaged.txt','wb')
	#np.savetxt(f, c, fmt='%.0f', delimiter=' ')
	#f.close()

#averageWeights(".", 40, "10.0")
averageWeights(".", 40, "20.0")
#averageWeights(".", 40, "28810.0")
