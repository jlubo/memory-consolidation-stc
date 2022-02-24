#######################################################################################
### Functions that serve to analyze the connectivity, early- and late-phase weights ###
###              in a network that contains multiple cell assemblies                ###
#######################################################################################

### Copyright 2019-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import sys
from utilityFunctions import cond_print

np.set_printoptions(precision=8, threshold=1e10, linewidth=200)
epsilon = 1e-9

# loadWeightMatrix
# Loads complete weight matrix from a file (only for excitatory neurons, though)
# filename: name of the file to read the data from
# N_pop: the number of neurons in the considered population
# h_0: initial synaptic weight and normalization factor for z
# return: the adjacency matrix, the early-phase weight matrix, the late-phase weight matrix, the firing rate vector
def loadWeightMatrix(filename, N_pop, h_0):
    global h
    global z
    global adj
    global v

    # read weight matrices and firing rates from file
    with open(filename) as f:
        rawdata = f.read()

    rawdata = rawdata.split('\n\n')
    rawmatrix_h = rawdata[0].split('\n')
    rawmatrix_z = rawdata[1].split('\n')
    rawmatrix_v = rawdata[2].split('\n')

    rows = len(rawmatrix_v)
    N_pop_row = int(round(np.sqrt(N_pop))) # one row of neurons when the population is aligned in a square

    if (rows != len(rawmatrix_v[0].split('\t\t'))) or (rows != N_pop_row):
        print('Data file error in "' + filename + '"')
        f.close()
        exit()

    v = np.zeros((N_pop_row,N_pop_row))
    h = np.zeros((N_pop,N_pop))
    z = np.zeros((N_pop,N_pop))

    for i in range(N_pop):
        if i < N_pop_row:
            value0 = rawmatrix_v[i].split('\t\t')
        value1 = rawmatrix_h[i].split('\t\t')
        value2 = rawmatrix_z[i].split('\t\t')

        for j in range(N_pop):
            if i < N_pop_row and j < N_pop_row:
                v[i][j] = float(value0[j])
            h[i][j] = float(value1[j])
            z[i][j] = h_0*float(value2[j])

    f.close()
    adj = (h > epsilon)

# loadAdjacencyMatrix
# Loads an adjacency matrix from a file into a numpy array
# filename: name of the file containing the adjacency matrix
def loadAdjacencyMatrix(filename):
	global adj

	adj = np.loadtxt(filename)

# getConnectivity
# Computes and prints the connectivity within the adjacency matrix
# pr [optional]: specifies if result is printed
# return: the connectivity as a value between 0 and 1
def getConnectivity(pr = True):
	nconn = np.sum(adj > 0)
	line = adj.shape[0]
	ret = nconn / (line**2 - line)
	cond_print(pr, "Connectivity of adjacency matrix: " + str(ret))
	return ret

# areConnected
# Returns true if neuron i and neuron j are connected (indices have to begin with 0) and prints the result
# i: neuron index
# j: neuron index
# pr [optional]: specifies if result is printed
# return: true or false
def areConnected(i, j, pr = True):
	global adj

	if adj[i][j] == 1:
		cond_print(pr, "Connection " + str(i) + "->" + str(j) + " does exist!")
		return True
	elif adj[i][j] == 0:
		cond_print(pr, "Connection " + str(i) + "->" + str(j) + " does NOT exist!")
		return False

# incomingConnections
# Prints and returns all the neurons from which neuron i receives inputs
# i: neuron index
# pr [optional]: specifies if result is printed
# return: array of neuron numbers
def incomingConnections(i, pr = True):
	global adj
	inc = np.where(adj[:, i] == 1)[0]
	cond_print(pr, "Incoming connections to " + str(i) + " (" + str(len(inc)) + "): \n" + str(inc))

	return inc

# outgoingConnections
# Prints and returns all the neurons to which neuron i provides input
# i: neuron index
# pr [optional]: specifies if result is printed
# return: array of neuron numbers
def outgoingConnections(i, pr = True):
	global adj
	out = np.where(adj[i, :] == 1)[0]
	cond_print(pr, "Outgoing connections from " + str(i) + " (" + str(len(out)) + "): \n" + str(out))

	return out

# incomingEarlyPhaseWeights
# Prints and returns all the early-phase synaptic weights incoming to neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
# return: array of early-phase weights
def incomingEarlyPhaseWeights(i, pr = True):
	global adj
	global h
	hi = h[:, i]
	adj_connections = (adj[:, i] == 1)
	inc = hi[adj_connections]
	cond_print(pr, "Incoming early-phase weights to " + str(i) + " (" + str(len(inc)) + "): \n" + str(inc))

	return inc

# outgoingEarlyPhaseWeights
# Prints and returns all the early-phase synaptic weights outgoing from neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
# return: array of early-phase weights
def outgoingEarlyPhaseWeights(i, pr = True):
	global adj
	global h
	hi = h[i, :]
	adj_connections = (adj[i, :] == 1)
	out = hi[adj_connections]
	cond_print(pr, "Outgoing early-phase weights from " + str(i) + " (" + str(len(out)) + "): \n" + str(out))

	return out

# earlyPhaseWeightsFromSet
# Prints and returns all the early-phase synaptic weights incoming to neuron i from a given set of neurons
# i: neuron index
# set: the set of presynaptic neurons
# return: array of early-phase weights
def earlyPhaseWeightsFromSet(i, set):
	global adj
	global h
	hi = h[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(hi)), set))
	inc_h = hi[adj_connections]

	return inc_h

# latePhaseWeightsFromSet
# Prints and returns all the late-phase synaptic weights incoming to neuron i from a given set of neurons
# i: neuron index
# set: the set of presynaptic neurons
# return: array of late-phase weights
def latePhaseWeightsFromSet(i, set):
	global adj
	global z
	zi = z[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(zi)), set))
	inc_z = zi[adj_connections]

	return inc_z

# meanEarlyPhaseWeight
# Returns the mean early-phase synaptic weight between two sets of neurons; prints the connectivity
# and the mean early-phase weight
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
# return: early-phase weight
def meanEarlyPhaseWeight(set, set2 = None, pr = True):
	summed_weight = 0
	connection_num = 0

	if set2 is None:
		set2 = set # consider internal connections within "set" if "set2" is not specified

	for n in set2:
		inc_weights = earlyPhaseWeightsFromSet(n, set)
		summed_weight += np.sum(inc_weights)
		connection_num += len(inc_weights)

	if connection_num > 0:
		ret = summed_weight / connection_num

		if set2 is None:
			cond_print(pr, "Self-connectivity: " + str(connection_num / (len(set)**2 - len(set))))
		else:
			cond_print(pr, "Connectivity: " + str(connection_num / (len(set)*len(set2))))
	else:
		ret = np.nan
		cond_print(pr, "Connectivity: none")

	cond_print(pr, "Mean early-phase weight: " + str(ret))

	return ret

# sdEarlyPhaseWeight
# Prints and returns the standard deviation of the early-phase synaptic weights between two sets of neurons
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
# return: early-phase weight
def sdEarlyPhaseWeight(set, set2 = None, pr = True):
	mean = meanEarlyPhaseWeight(set, set2, False)
	summed_qu_dev = 0
	connection_num = 0

	if set2 is None:
		set2 = set # consider internal connections within "set" if "set2" is not specified

	for n in set2:
		qu_devs = np.power(mean - earlyPhaseWeightsFromSet(n, set), 2) # quadratic deviations from mean
		summed_qu_dev += np.sum(qu_devs)
		connection_num += len(qu_devs)

	if connection_num > 1:
		ret = np.sqrt(summed_qu_dev / (connection_num-1))
	else:
		ret = np.nan

	cond_print(pr, "Std. dev. of early-phase weight: " + str(ret))

	return ret


# meanLatePhaseWeight
# Returns the mean late-phase synaptic weight between two sets of neurons; prints the connectivity
# and the mean late-phase weight
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
# return: late-phase weight
def meanLatePhaseWeight(set, set2 = None, pr = True):
	summed_weight = 0
	connection_num = 0

	if set2 is None:
		set2 = set # consider internal connections within "set" if "set2" is not specified

	for n in set2:
		inc_weights = latePhaseWeightsFromSet(n, set)
		summed_weight += np.sum(inc_weights)
		connection_num += len(inc_weights)

	if connection_num > 0:
		ret = summed_weight / connection_num
	else:
		ret = np.nan

	if connection_num > 0:
		ret = summed_weight / connection_num

		if set2 is None:
			cond_print(pr, "Self-connectivity: " + str(connection_num / (len(set)**2 - len(set))))
		else:
			cond_print(pr, "Connectivity: " + str(connection_num / (len(set)*len(set2))))
	else:
		ret = np.nan
		cond_print(pr, "Connectivity: none")

	cond_print(pr, "Mean late-phase weight: " + str(ret))

	return ret

# sdLatePhaseWeight
# Prints and returns the standard deviation of the late-phase synaptic weights between two sets of neurons
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
# return: late-phase weight
def sdLatePhaseWeight(set, set2 = [], pr = True):
	mean = meanLatePhaseWeight(set, set2, False)
	summed_qu_dev = 0
	connection_num = 0

	if len(set2) <= 0:
		set2 = set # consider internal connections within "set" if "set2" is not specified

	for n in set2:
		qu_devs = np.power(mean - latePhaseWeightsFromSet(n, set), 2) # quadratic deviations from mean
		summed_qu_dev += np.sum(qu_devs)
		connection_num += len(qu_devs)

	if connection_num > 1:
		ret = np.sqrt(summed_qu_dev / (connection_num-1))
	else:
		ret = np.nan

	cond_print(pr, "Std. dev. of late-phase weight: " + str(ret))

	return ret

# meanCoreWeights
# Computes the mean weights in cores (including intersections) and appends them, together with
# the time for readout, to a file
# ts: timestamp of the data file to read
# time_for_readout: the time of the data file to read
# coreA: array of indices of the first cell assembly (core) neurons
# coreB: array of indices of the second cell assembly (core) neurons
# coreC: array of indices of the third cell assembly (core) neurons
# N_pop: the number of neurons in the considered population
# h_0: initial synaptic weight and normalization factor for z
# pr [optional]: specifies if result shall be printed
def meanCoreWeights(ts, time_for_readout, coreA, coreB, coreC, N_pop, h_0, pr = True):
	
	cond_print(pr, "##############################################")
	cond_print(pr, "At time", time_for_readout)
	loadWeightMatrix(ts + "_net_" + time_for_readout + ".txt", N_pop, h_0)
	
	# early-phase weights
	f = open("cores_mean_early_weights.txt", "a")
	f.write(time_for_readout + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "A -> A:")
	hm_A = meanEarlyPhaseWeight(coreA, coreA, pr)
	hsd_A = sdEarlyPhaseWeight(coreA, coreA, pr)
	f.write(str(hm_A) + "\t\t")
	f.write(str(hsd_A) + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "B -> B:")
	hm_B = meanEarlyPhaseWeight(coreB, coreB, pr)
	hsd_B = sdEarlyPhaseWeight(coreB, coreB, pr)
	f.write(str(hm_B) + "\t\t")
	f.write(str(hsd_B) + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "C -> C:")
	hm_C = meanEarlyPhaseWeight(coreC, coreC, pr)
	hsd_C = sdEarlyPhaseWeight(coreC, coreC, pr)
	f.write(str(hm_C) + "\t\t")
	f.write(str(hsd_C) + "\n")
	f.close()
	
	# late-phase weights
	f = open("cores_mean_late_weights.txt", "a")
	f.write(time_for_readout + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "A -> A:")
	zm_A = meanLatePhaseWeight(coreA, coreA, pr)
	zsd_A = sdLatePhaseWeight(coreA, coreA, pr)
	f.write(str(zm_A) + "\t\t")
	f.write(str(zsd_A) + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "B -> B:")
	zm_B = meanLatePhaseWeight(coreB, coreB, pr)
	zsd_B = sdLatePhaseWeight(coreB, coreB, pr)
	f.write(str(zm_B) + "\t\t")
	f.write(str(zsd_B) + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "C -> C:")
	zm_C = meanLatePhaseWeight(coreC, coreC, pr)
	zsd_C = sdLatePhaseWeight(coreC, coreC, pr)
	f.write(str(zm_C) + "\t\t")
	f.write(str(zsd_C) + "\n")
	f.close()
	
	# total weights
	f = open("cores_mean_tot_weights.txt", "a")
	f.write(time_for_readout + "\t\t")
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "A -> A:")
	wm_A = hm_A + zm_A
	wsd_A = np.sqrt(hsd_A**2 + zsd_A**2)
	f.write(str(wm_A) + "\t\t")
	f.write(str(wsd_A) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(wm_A))
	cond_print(pr, "Std. dev. of total weight: " + str(wsd_A))
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "B -> B:")
	wm_B = hm_B + zm_B
	wsd_B = np.sqrt(hsd_B**2 + zsd_B**2)
	f.write(str(wm_B) + "\t\t")
	f.write(str(wsd_B) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(wm_B))
	cond_print(pr, "Std. dev. of total weight: " + str(wsd_B))
	
	cond_print(pr, "--------------------------------")
	cond_print(pr, "C -> C:")
	wm_C = hm_C + zm_C
	wsd_C = np.sqrt(hsd_C**2 + zsd_C**2)
	f.write(str(wm_C) + "\t\t")
	f.write(str(wsd_C) + "\n")
	cond_print(pr, "Mean total weight: " + str(wm_C))
	cond_print(pr, "Std. dev. of total weight: " + str(wsd_C))
	f.close()

# meanWeightMatrix
# Computes the abstract mean weight matrix and write the outcome to a file
# ts: timestamp of the data file to read
# time_for_readout: the time of the data file to read
# coreA: array of indices of the first cell assembly (core) neurons
# coreB: array of indices of the second cell assembly (core) neurons
# coreC: array of indices of the third cell assembly (core) neurons
# N_pop: the number of neurons in the considered population
# h_0: initial synaptic weight and normalization factor for z
# pr [optional]: specifies if result shall be printed
def meanWeightMatrix(ts, time_for_readout, coreA, coreB, coreC, N_pop, h_0, pr = True):
	# define the whole considered population
	all = np.arange(N_pop)

	# determine masks of whole cores
	mask_coreA = np.in1d(all, coreA) # boolean matrix of neurons in whole core A
	mask_coreB = np.in1d(all, coreB) # boolean matrix of neurons in whole core B
	mask_coreC = np.in1d(all, coreC) # boolean matrix of neurons in whole core C

	# determine exclusive intersections
	mask_I_AB = np.logical_and( np.logical_and(mask_coreA, mask_coreB), np.logical_not(mask_coreC) )
	mask_I_AC = np.logical_and( np.logical_and(mask_coreA, mask_coreC), np.logical_not(mask_coreB) )
	mask_I_BC = np.logical_and( np.logical_and(mask_coreB, mask_coreC), np.logical_not(mask_coreA) )
	mask_I_ABC = np.logical_and( mask_coreA, np.logical_and(mask_coreB, mask_coreC) )
	I_AB = all[mask_I_AB]
	I_AC = all[mask_I_AC]
	I_BC = all[mask_I_BC]
	I_ABC = all[mask_I_ABC]

	# determine exclusive cores by removing exclusive intersections from whole cores
	exA = all[np.logical_and(mask_coreA, \
	                         np.logical_and(np.logical_not(mask_I_AB), \
	                                        np.logical_and(np.logical_not(mask_I_AC), np.logical_not(mask_I_ABC))))]
	exB = all[np.logical_and(mask_coreB, \
	                         np.logical_and(np.logical_not(mask_I_AB), \
	                                        np.logical_and(np.logical_not(mask_I_BC), np.logical_not(mask_I_ABC))))]
	exC = all[np.logical_and(mask_coreC, \
	                         np.logical_and(np.logical_not(mask_I_AC), \
	                                        np.logical_and(np.logical_not(mask_I_BC), np.logical_not(mask_I_ABC))))]

	# determine control subpopulation
	control = all[np.logical_not(np.in1d(all, np.concatenate([I_ABC, I_AB, I_AC, I_BC,
	                                                          coreA, coreB, coreC])))]  # all neurons that are not part of an assembly

	cond_print(pr, "##############################################")
	cond_print(pr, "At time", time_for_readout)
	loadWeightMatrix(ts + "_net_" + time_for_readout + ".txt", N_pop, h_0)
	f = open("mean_tot_weights_" + time_for_readout + ".txt", "w")
	fsd = open("sd_tot_weights_" + time_for_readout + ".txt", "w")

	### OUTGOING FROM I_ABC ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> I_ABC:")
	hm = meanEarlyPhaseWeight(I_ABC, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, I_ABC, pr)
	zm = meanLatePhaseWeight(I_ABC, I_ABC, pr)
	zsd = sdLatePhaseWeight(I_ABC, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> I_AC:")
	hm = meanEarlyPhaseWeight(I_ABC, I_AC, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, I_AC, pr)
	zm = meanLatePhaseWeight(I_ABC, I_AC, pr)
	zsd = sdLatePhaseWeight(I_ABC, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> ~A:")
	hm = meanEarlyPhaseWeight(I_ABC, exA, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, exA, pr)
	zm = meanLatePhaseWeight(I_ABC, exA, pr)
	zsd = sdLatePhaseWeight(I_ABC, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> I_AB:")
	hm = meanEarlyPhaseWeight(I_ABC, I_AB, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, I_AB, pr)
	zm = meanLatePhaseWeight(I_ABC, I_AB, pr)
	zsd = sdLatePhaseWeight(I_ABC, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> ~B:")
	hm = meanEarlyPhaseWeight(I_ABC, exB, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, exB, pr)
	zm = meanLatePhaseWeight(I_ABC, exB, pr)
	zsd = sdLatePhaseWeight(I_ABC, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> I_BC:")
	hm = meanEarlyPhaseWeight(I_ABC, I_BC, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, I_BC, pr)
	zm = meanLatePhaseWeight(I_ABC, I_BC, pr)
	zsd = sdLatePhaseWeight(I_ABC, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> ~C:")
	hm = meanEarlyPhaseWeight(I_ABC, exC, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, exC, pr)
	zm = meanLatePhaseWeight(I_ABC, exC, pr)
	zsd = sdLatePhaseWeight(I_ABC, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_ABC -> Control:")
	hm = meanEarlyPhaseWeight(I_ABC, control, pr)
	hsd = sdEarlyPhaseWeight(I_ABC, control, pr)
	zm = meanLatePhaseWeight(I_ABC, control, pr)
	zsd = sdLatePhaseWeight(I_ABC, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM I_AC ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> I_ABC:")
	hm = meanEarlyPhaseWeight(I_AC, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(I_AC, I_ABC, pr)
	zm = meanLatePhaseWeight(I_AC, I_ABC, pr)
	zsd = sdLatePhaseWeight(I_AC, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> I_AC:")
	hm = meanEarlyPhaseWeight(I_AC, I_AC, pr)
	hsd = sdEarlyPhaseWeight(I_AC, I_AC, pr)
	zm = meanLatePhaseWeight(I_AC, I_AC, pr)
	zsd = sdLatePhaseWeight(I_AC, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> ~A:")
	hm = meanEarlyPhaseWeight(I_AC, exA, pr)
	hsd = sdEarlyPhaseWeight(I_AC, exA, pr)
	zm = meanLatePhaseWeight(I_AC, exA, pr)
	zsd = sdLatePhaseWeight(I_AC, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> I_AB:")
	hm = meanEarlyPhaseWeight(I_AC, I_AB, pr)
	hsd = sdEarlyPhaseWeight(I_AC, I_AB, pr)
	zm = meanLatePhaseWeight(I_AC, I_AB, pr)
	zsd = sdLatePhaseWeight(I_AC, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> ~B:")
	hm = meanEarlyPhaseWeight(I_AC, exB, pr)
	hsd = sdEarlyPhaseWeight(I_AC, exB, pr)
	zm = meanLatePhaseWeight(I_AC, exB, pr)
	zsd = sdLatePhaseWeight(I_AC, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> I_BC:")
	hm = meanEarlyPhaseWeight(I_AC, I_BC, pr)
	hsd = sdEarlyPhaseWeight(I_AC, I_BC, pr)
	zm = meanLatePhaseWeight(I_AC, I_BC, pr)
	zsd = sdLatePhaseWeight(I_AC, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> ~C:")
	hm = meanEarlyPhaseWeight(I_AC, exC, pr)
	hsd = sdEarlyPhaseWeight(I_AC, exC, pr)
	zm = meanLatePhaseWeight(I_AC, exC, pr)
	zsd = sdLatePhaseWeight(I_AC, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AC -> Control:")
	hm = meanEarlyPhaseWeight(I_AC, control, pr)
	hsd = sdEarlyPhaseWeight(I_AC, control, pr)
	zm = meanLatePhaseWeight(I_AC, control, pr)
	zsd = sdLatePhaseWeight(I_AC, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM ~A ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> I_ABC:")
	hm = meanEarlyPhaseWeight(exA, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(exA, I_ABC, pr)
	zm = meanLatePhaseWeight(exA, I_ABC, pr)
	zsd = sdLatePhaseWeight(exA, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> I_AC:")
	hm = meanEarlyPhaseWeight(exA, I_AC, pr)
	hsd = sdEarlyPhaseWeight(exA, I_AC, pr)
	zm = meanLatePhaseWeight(exA, I_AC, pr)
	zsd = sdLatePhaseWeight(exA, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> ~A:")
	hm = meanEarlyPhaseWeight(exA, exA, pr)
	hsd = sdEarlyPhaseWeight(exA, exA, pr)
	zm = meanLatePhaseWeight(exA, exA, pr)
	zsd = sdLatePhaseWeight(exA, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> I_AB:")
	hm = meanEarlyPhaseWeight(exA, I_AB, pr)
	hsd = sdEarlyPhaseWeight(exA, I_AB, pr)
	zm = meanLatePhaseWeight(exA, I_AB, pr)
	zsd = sdLatePhaseWeight(exA, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> ~B:")
	hm = meanEarlyPhaseWeight(exA, exB, pr)
	hsd = sdEarlyPhaseWeight(exA, exB, pr)
	zm = meanLatePhaseWeight(exA, exB, pr)
	zsd = sdLatePhaseWeight(exA, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> I_BC:")
	hm = meanEarlyPhaseWeight(exA, I_BC, pr)
	hsd = sdEarlyPhaseWeight(exA, I_BC, pr)
	zm = meanLatePhaseWeight(exA, I_BC, pr)
	zsd = sdLatePhaseWeight(exA, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> ~C:")
	hm = meanEarlyPhaseWeight(exA, exC, pr)
	hsd = sdEarlyPhaseWeight(exA, exC, pr)
	zm = meanLatePhaseWeight(exA, exC, pr)
	zsd = sdLatePhaseWeight(exA, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~A -> Control:")
	hm = meanEarlyPhaseWeight(exA, control, pr)
	hsd = sdEarlyPhaseWeight(exA, control, pr)
	zm = meanLatePhaseWeight(exA, control, pr)
	zsd = sdLatePhaseWeight(exA, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM I_AB ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> I_ABC:")
	hm = meanEarlyPhaseWeight(I_AB, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(I_AB, I_ABC, pr)
	zm = meanLatePhaseWeight(I_AB, I_ABC, pr)
	zsd = sdLatePhaseWeight(I_AB, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> I_AC:")
	hm = meanEarlyPhaseWeight(I_AB, I_AC, pr)
	hsd = sdEarlyPhaseWeight(I_AB, I_AC, pr)
	zm = meanLatePhaseWeight(I_AB, I_AC, pr)
	zsd = sdLatePhaseWeight(I_AB, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> ~A:")
	hm = meanEarlyPhaseWeight(I_AB, exA, pr)
	hsd = sdEarlyPhaseWeight(I_AB, exA, pr)
	zm = meanLatePhaseWeight(I_AB, exA, pr)
	zsd = sdLatePhaseWeight(I_AB, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> I_AB:")
	hm = meanEarlyPhaseWeight(I_AB, I_AB, pr)
	hsd = sdEarlyPhaseWeight(I_AB, I_AB, pr)
	zm = meanLatePhaseWeight(I_AB, I_AB, pr)
	zsd = sdLatePhaseWeight(I_AB, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> ~B:")
	hm = meanEarlyPhaseWeight(I_AB, exB, pr)
	hsd = sdEarlyPhaseWeight(I_AB, exB, pr)
	zm = meanLatePhaseWeight(I_AB, exB, pr)
	zsd = sdLatePhaseWeight(I_AB, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> I_BC:")
	hm = meanEarlyPhaseWeight(I_AB, I_BC, pr)
	hsd = sdEarlyPhaseWeight(I_AB, I_BC, pr)
	zm = meanLatePhaseWeight(I_AB, I_BC, pr)
	zsd = sdLatePhaseWeight(I_AB, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> ~C:")
	hm = meanEarlyPhaseWeight(I_AB, exC, pr)
	hsd = sdEarlyPhaseWeight(I_AB, exC, pr)
	zm = meanLatePhaseWeight(I_AB, exC, pr)
	zsd = sdLatePhaseWeight(I_AB, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_AB -> Control:")
	hm = meanEarlyPhaseWeight(I_AB, control, pr)
	hsd = sdEarlyPhaseWeight(I_AB, control, pr)
	zm = meanLatePhaseWeight(I_AB, control, pr)
	zsd = sdLatePhaseWeight(I_AB, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM ~B ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> I_ABC:")
	hm = meanEarlyPhaseWeight(exB, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(exB, I_ABC, pr)
	zm = meanLatePhaseWeight(exB, I_ABC, pr)
	zsd = sdLatePhaseWeight(exB, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> I_AC:")
	hm = meanEarlyPhaseWeight(exB, I_AC, pr)
	hsd = sdEarlyPhaseWeight(exB, I_AC, pr)
	zm = meanLatePhaseWeight(exB, I_AC, pr)
	zsd = sdLatePhaseWeight(exB, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> ~A:")
	hm = meanEarlyPhaseWeight(exB, exA, pr)
	hsd = sdEarlyPhaseWeight(exB, exA, pr)
	zm = meanLatePhaseWeight(exB, exA, pr)
	zsd = sdLatePhaseWeight(exB, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> I_AB:")
	hm = meanEarlyPhaseWeight(exB, I_AB, pr)
	hsd = sdEarlyPhaseWeight(exB, I_AB, pr)
	zm = meanLatePhaseWeight(exB, I_AB, pr)
	zsd = sdLatePhaseWeight(exB, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> ~B:")
	hm = meanEarlyPhaseWeight(exB, exB, pr)
	hsd = sdEarlyPhaseWeight(exB, exB, pr)
	zm = meanLatePhaseWeight(exB, exB, pr)
	zsd = sdLatePhaseWeight(exB, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> I_BC:")
	hm = meanEarlyPhaseWeight(exB, I_BC, pr)
	hsd = sdEarlyPhaseWeight(exB, I_BC, pr)
	zm = meanLatePhaseWeight(exB, I_BC, pr)
	zsd = sdLatePhaseWeight(exB, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> ~C:")
	hm = meanEarlyPhaseWeight(exB, exC, pr)
	hsd = sdEarlyPhaseWeight(exB, exC, pr)
	zm = meanLatePhaseWeight(exB, exC, pr)
	zsd = sdLatePhaseWeight(exB, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~B -> Control:")
	hm = meanEarlyPhaseWeight(exB, control, pr)
	hsd = sdEarlyPhaseWeight(exB, control, pr)
	zm = meanLatePhaseWeight(exB, control, pr)
	zsd = sdLatePhaseWeight(exB, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM I_BC ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> I_ABC:")
	hm = meanEarlyPhaseWeight(I_BC, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(I_BC, I_ABC, pr)
	zm = meanLatePhaseWeight(I_BC, I_ABC, pr)
	zsd = sdLatePhaseWeight(I_BC, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> I_AC:")
	hm = meanEarlyPhaseWeight(I_BC, I_AC, pr)
	hsd = sdEarlyPhaseWeight(I_BC, I_AC, pr)
	zm = meanLatePhaseWeight(I_BC, I_AC, pr)
	zsd = sdLatePhaseWeight(I_BC, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> ~A:")
	hm = meanEarlyPhaseWeight(I_BC, exA, pr)
	hsd = sdEarlyPhaseWeight(I_BC, exA, pr)
	zm = meanLatePhaseWeight(I_BC, exA, pr)
	zsd = sdLatePhaseWeight(I_BC, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> I_AB:")
	hm = meanEarlyPhaseWeight(I_BC, I_AB, pr)
	hsd = sdEarlyPhaseWeight(I_BC, I_AB, pr)
	zm = meanLatePhaseWeight(I_BC, I_AB, pr)
	zsd = sdLatePhaseWeight(I_BC, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> ~B:")
	hm = meanEarlyPhaseWeight(I_BC, exB, pr)
	hsd = sdEarlyPhaseWeight(I_BC, exB, pr)
	zm = meanLatePhaseWeight(I_BC, exB, pr)
	zsd = sdLatePhaseWeight(I_BC, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> I_BC:")
	hm = meanEarlyPhaseWeight(I_BC, I_BC, pr)
	hsd = sdEarlyPhaseWeight(I_BC, I_BC, pr)
	zm = meanLatePhaseWeight(I_BC, I_BC, pr)
	zsd = sdLatePhaseWeight(I_BC, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> ~C:")
	hm = meanEarlyPhaseWeight(I_BC, exC, pr)
	hsd = sdEarlyPhaseWeight(I_BC, exC, pr)
	zm = meanLatePhaseWeight(I_BC, exC, pr)
	zsd = sdLatePhaseWeight(I_BC, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "I_BC -> Control:")
	hm = meanEarlyPhaseWeight(I_BC, control, pr)
	hsd = sdEarlyPhaseWeight(I_BC, control, pr)
	zm = meanLatePhaseWeight(I_BC, control, pr)
	zsd = sdLatePhaseWeight(I_BC, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))


	### OUTGOING FROM ~C ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> I_ABC:")
	hm = meanEarlyPhaseWeight(exC, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(exC, I_ABC, pr)
	zm = meanLatePhaseWeight(exC, I_ABC, pr)
	zsd = sdLatePhaseWeight(exC, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> I_AC:")
	hm = meanEarlyPhaseWeight(exC, I_AC, pr)
	hsd = sdEarlyPhaseWeight(exC, I_AC, pr)
	zm = meanLatePhaseWeight(exC, I_AC, pr)
	zsd = sdLatePhaseWeight(exC, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> ~A:")
	hm = meanEarlyPhaseWeight(exC, exA, pr)
	hsd = sdEarlyPhaseWeight(exC, exA, pr)
	zm = meanLatePhaseWeight(exC, exA, pr)
	zsd = sdLatePhaseWeight(exC, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> I_AB:")
	hm = meanEarlyPhaseWeight(exC, I_AB, pr)
	hsd = sdEarlyPhaseWeight(exC, I_AB, pr)
	zm = meanLatePhaseWeight(exC, I_AB, pr)
	zsd = sdLatePhaseWeight(exC, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> ~B:")
	hm = meanEarlyPhaseWeight(exC, exB, pr)
	hsd = sdEarlyPhaseWeight(exC, exB, pr)
	zm = meanLatePhaseWeight(exC, exB, pr)
	zsd = sdLatePhaseWeight(exC, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> I_BC:")
	hm = meanEarlyPhaseWeight(exC, I_BC, pr)
	hsd = sdEarlyPhaseWeight(exC, I_BC, pr)
	zm = meanLatePhaseWeight(exC, I_BC, pr)
	zsd = sdLatePhaseWeight(exC, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> ~C:")
	hm = meanEarlyPhaseWeight(exC, exC, pr)
	hsd = sdEarlyPhaseWeight(exC, exC, pr)
	zm = meanLatePhaseWeight(exC, exC, pr)
	zsd = sdLatePhaseWeight(exC, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "~C -> Control:")
	hm = meanEarlyPhaseWeight(exC, control, pr)
	hsd = sdEarlyPhaseWeight(exC, control, pr)
	zm = meanLatePhaseWeight(exC, control, pr)
	zsd = sdLatePhaseWeight(exC, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	### OUTGOING FROM Control ###########################
	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> I_ABC:")
	hm = meanEarlyPhaseWeight(control, I_ABC, pr)
	hsd = sdEarlyPhaseWeight(control, I_ABC, pr)
	zm = meanLatePhaseWeight(control, I_ABC, pr)
	zsd = sdLatePhaseWeight(control, I_ABC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> I_AC:")
	hm = meanEarlyPhaseWeight(control, I_AC, pr)
	hsd = sdEarlyPhaseWeight(control, I_AC, pr)
	zm = meanLatePhaseWeight(control, I_AC, pr)
	zsd = sdLatePhaseWeight(control, I_AC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> ~A:")
	hm = meanEarlyPhaseWeight(control, exA, pr)
	hsd = sdEarlyPhaseWeight(control, exA, pr)
	zm = meanLatePhaseWeight(control, exA, pr)
	zsd = sdLatePhaseWeight(control, exA, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> I_AB:")
	hm = meanEarlyPhaseWeight(control, I_AB, pr)
	hsd = sdEarlyPhaseWeight(control, I_AB, pr)
	zm = meanLatePhaseWeight(control, I_AB, pr)
	zsd = sdLatePhaseWeight(control, I_AB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> ~B:")
	hm = meanEarlyPhaseWeight(control, exB, pr)
	hsd = sdEarlyPhaseWeight(control, exB, pr)
	zm = meanLatePhaseWeight(control, exB, pr)
	zsd = sdLatePhaseWeight(control, exB, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> I_BC:")
	hm = meanEarlyPhaseWeight(control, I_BC, pr)
	hsd = sdEarlyPhaseWeight(control, I_BC, pr)
	zm = meanLatePhaseWeight(control, I_BC, pr)
	zsd = sdLatePhaseWeight(control, I_BC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> ~C:")
	hm = meanEarlyPhaseWeight(control, exC, pr)
	hsd = sdEarlyPhaseWeight(control, exC, pr)
	zm = meanLatePhaseWeight(control, exC, pr)
	zsd = sdLatePhaseWeight(control, exC, pr)
	f.write(str(hm + zm) + "\t\t")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\t\t")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	cond_print(pr, "--------------------------------")
	cond_print(pr, "Control -> Control:")
	hm = meanEarlyPhaseWeight(control, control, pr)
	hsd = sdEarlyPhaseWeight(control, control, pr)
	zm = meanLatePhaseWeight(control, control, pr)
	zsd = sdLatePhaseWeight(control, control, pr)
	f.write(str(hm + zm) + "\n")
	fsd.write(str(np.sqrt(hsd**2 + zsd**2)) + "\n")
	cond_print(pr, "Mean total weight: " + str(hm + zm))

	f.close()
	fsd.close()


# printMeanWeightsSingleCA
# Computes and prints mean and standard deviation of CA, outgoing, incoming, and control weights
# ts: timestamp of the data file to read
# time_for_readout: the time of the data file to read
# core: array of indices of the cell assembly (core) neurons
# N_pop: the number of neurons in the considered population
# h_0: initial synaptic weight and normalization factor for z
def printMeanWeightsSingleCA(ts, time_for_readout, core, N_pop, h_0):

	all = np.arange(N_pop) # the whole considered population
	noncore = all[np.logical_not(np.in1d(all, core))] # the neurons outside the core

	print("##############################################")
	print("At time", time_for_readout)
	loadWeightMatrix(ts + "_net_" + time_for_readout + ".txt", N_pop, h_0)

	print("--------------------------------")
	print("Core -> core ('CA'):")
	hm = meanEarlyPhaseWeight(core)
	sdEarlyPhaseWeight(core)
	zm = meanLatePhaseWeight(core)
	sdLatePhaseWeight(core)
	print("Mean total weight: " + str(hm + zm))

	print("--------------------------------")
	print("Core -> non-core ('outgoing'):")
	hm = meanEarlyPhaseWeight(core, noncore)
	sdEarlyPhaseWeight(core, noncore)
	zm = meanLatePhaseWeight(core, noncore)
	sdLatePhaseWeight(core, noncore)
	print("Mean total weight: " + str(hm + zm))

	print("--------------------------------")
	print("Non-core -> core ('incoming'):")
	hm = meanEarlyPhaseWeight(noncore, core)
	sdEarlyPhaseWeight(noncore, core)
	zm = meanLatePhaseWeight(noncore, core)
	sdLatePhaseWeight(noncore, core)
	print("Mean total weight: " + str(hm + zm))

	print("--------------------------------")
	print("Non-core -> non-core ('non-CA'):")
	hm = meanEarlyPhaseWeight(noncore)
	sdEarlyPhaseWeight(noncore)
	zm = meanLatePhaseWeight(noncore)
	sdLatePhaseWeight(noncore)
	print("Mean total weight: " + str(hm + zm))


##############################################################################################
# main
# Reads datasets from two simulations and computes mean CA, outgoing, incoming and control weights (early- and late-phase)
# as for Luboeinski and Tetzlaff, 2021 (https://doi.org/10.1038/s42003-021-01778-y)
# argv[]: timestamps of two simulations

# example call from shell: python3 adjacencyFunctions.py "19-11-28_21-07-55" "19-11-28_22-10-17"
if __name__ == "__main__":

	if len(sys.argv) < 3:
		print("Not enough arguments provided!")
		exit()
	else:
		ts1 = str(sys.argv[1]) # timestamp for simulation data before consolidation
		ts2 = str(sys.argv[2]) # timestamp for simulation data after consolidation

	core = np.arange(150) # define the cell assembly core
	N_pop = 1600 # number of neurons in the population
	h_0 = 4.20075 # initial/median synaptic weight

	print("##############################################")
	print("Before 10s-recall:")
	printMeanWeightsSingleCA(ts1, "20.0", core, N_pop, h_0)

	print("##############################################")
	print("After 10s-recall:")
	printMeanWeightsSingleCA(ts1, "20.1", core, N_pop, h_0)

	print("##############################################")
	print("Before 8h-recall:")
	printMeanWeightsSingleCA(ts2, "28810.0", core, N_pop, h_0)

	print("##############################################")
	print("After 8h-recall:")
	printMeanWeightsSingleCA(ts2, "28810.1", core, N_pop, h_0)
