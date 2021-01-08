#######################################################################################
### Functions that serve to analyze the connectivity, early- and late-phase weights ###
###                in a network that contains a cell assembly                       ###
#######################################################################################

### Copyright 2019-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
import sys

# the assembly core (the neurons that receive learning stimulation) and the recall neurons (the neurons that receive recall stimulation)
core = np.arange(350)
recall_fraction = 0.5
core_recall = core[0:int(np.floor(recall_fraction*core.shape[0]))]
core_norecall = core[np.logical_not(np.in1d(core, core_recall))]

Nl = 40 # number of excitatory neurons in one line (Nl^2 is the total number of excitatory neurons)

all = np.arange(Nl**2)
noncore = all[np.logical_not(np.in1d(all, core))]

np.set_printoptions(precision=8, threshold=1e10, linewidth=200)
epsilon = 1e-9

# loadWeightMatrix
# Loads complete weight matrix from a file (only for excitatory neurons, though)
# filename: name of the file to read the data from
# return: the adjacency matrix, the early-phase weight matrix, the late-phase weight matrix, the firing rate vector
def loadWeightMatrix(filename):

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

    if (rows != len(rawmatrix_v[0].split('\t\t'))) or (rows != Nl):
        print('Data file error in "' + filename + '"')
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
def getConnectivity(pr = True):
	nconn = np.sum(adj > 0)
	line = adj.shape[0]
	ret = nconn / (line**2 - line)
	if (pr):
		print("Connectivity of adjacency matrix: " + str(ret))
	return ret

# getConnectivityInCore
# Computes and prints the connectivity within the core
# pr [optional]: specifies if result is printed
def getConnectivityInCore(pr = True):
	connections_within_core = 0
	N_core = len(core)

	for n in core:
		connections_within_core += len(connectionsFromCore(n, False))

	ret = connections_within_core / (N_core**2 - N_core)

	if pr:
		print("Connectivity within core: " + str(ret))

	return ret

# areConnected
# Returns true if neuron i and neuron j are connected (indices have to begin with 0) and prints the result
# i: neuron index
# j: neuron index
# pr [optional]: specifies if result is printed
def areConnected(i, j, pr = True):
	global adj

	if adj[i][j] == 1:
		if pr:
			print("Connection " + str(i) + "->" + str(j) + " does exist!")
		return True
	elif adj[i][j] == 0:
		if pr:
			print("Connection " + str(i) + "->" + str(j) + " does NOT exist!")
		return False

# incomingConnections
# Prints and returns all the neurons from which neuron i receives inputs
# i: neuron index
# pr [optional]: specifies if result is printed
def incomingConnections(i, pr = True):
	global adj
	inc = np.where(adj[:, i] == 1)[0]
	if pr:
		print("Incoming connections to " + str(i) + " (" + str(len(inc)) + "): \n" + str(inc))

	return inc

# outgoingConnections
# Prints and returns all the neurons to which neuron i provides input
# i: neuron index
# pr [optional]: specifies if result is printed
def outgoingConnections(i, pr = True):
	global adj
	out = np.where(adj[i, :] == 1)[0]
	if pr:
		print("Outgoing connections from " + str(i) + " (" + str(len(out)) + "): \n" + str(out))

	return out

# incomingEarlyPhaseWeights
# Prints and returns all the early-phase synaptic weights incoming to neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
def incomingEarlyPhaseWeights(i, pr = True):
	global adj
	global h
	hi = h[:, i]
	adj_connections = (adj[:, i] == 1)
	inc = hi[adj_connections]
	if pr:
		print("Incoming early-phase weights to " + str(i) + " (" + str(len(inc)) + "): \n" + str(inc))

	return inc

# outgoingEarlyPhaseWeights
# Prints and returns all the early-phase synaptic weights outgoing from neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
def outgoingEarlyPhaseWeights(i, pr = True):
	global adj
	global h
	hi = h[i, :]
	adj_connections = (adj[i, :] == 1)
	out = hi[adj_connections]
	if pr:
		print("Outgoing early-phase weights from " + str(i) + " (" + str(len(out)) + "): \n" + str(out))

	return out

# connectionsFromCore
# Prints and returns all the core neurons from which neuron i receives inputs
# i: neuron index
# pr [optional]: specifies if result is printed
def connectionsFromCore(i, pr = True):
	global adj
	cfc = core[np.in1d(core, incomingConnections(i, False))].flatten()
	if pr:
		print("Connections from core to neuron " + str(i) + " (" + str(len(cfc)) + "): \n" + str(cfc))

	return cfc

# connectionsFromRSCore
# Prints and returns all the recall-stimulated core neurons from which neuron i receives inputs
# i: neuron index
# pr [optional]: specifies if result is printed
def connectionsFromRSCore(i, pr = True):
	global adj
	cfrsc = core_recall[np.in1d(core_recall, incomingConnections(i, False))].flatten()
	if pr:
		print("Connections from recall-stim. core to neuron " + str(i) + " (" + str(len(cfrsc)) + "): \n" + str(cfrsc))

	return cfrsc

# earlyPhaseWeightsFromSet
# Prints and returns all the early-phase synaptic weights incoming to neuron i from a given set of neurons
# i: neuron index
# set: the set of presynaptic neurons
def earlyPhaseWeightsFromSet(i, set):
	global adj
	global h
	hi = h[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(hi)), set))
	inc_h = hi[adj_connections]

	return inc_h

# earlyPhaseWeightsFromCore
# Prints and returns all the early-phase synaptic weights incoming to neuron i from core neurons
# i: neuron index
# pr [optional]: specifies if result is printed
def earlyPhaseWeightsFromCore(i, pr = True):
	inc_h = earlyPhaseWeightsFromSet(i, core)
	if pr:
		print("Incoming early-phase weights from core to neuron " + str(i) + " (" + str(len(inc_h)) + "): \n" + str(inc_h))

	return inc_h

# earlyPhaseWeightsFromRSCore
# Prints and returns all the early-phase synaptic weights incoming to neuron i from recall-stimulated core neurons
# i: neuron index
# pr [optional]: specifies if result is printed
def earlyPhaseWeightsFromRSCore(i, pr = True):
	inc_h = earlyPhaseWeightsFromSet(i, core_recall)
	if pr:
		print("Incoming early-phase weights from recall-stim. core to neuron " + str(i) + " (" + str(len(inc_h)) + "): \n" + str(inc_h))

	return inc_h

# earlyPhaseWeightsToCore
# Prints and returns all the early-phase synaptic weights outgoing to core neurons from neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
def earlyPhaseWeightsToCore(i, pr = True):
	global adj
	global h
	hi = h[i, :]
	adj_connections = np.logical_and(adj[i, :] == 1, np.in1d(np.arange(len(hi)), core))
	out_h = hi[adj_connections]
	if pr:
		print("Outgoing early-phase weights from neuron " + str(i) + " to core (" + str(len(out_h)) + "): \n" + str(out_h))

	return out_h

# latePhaseWeightsFromSet
# Prints and returns all the late-phase synaptic weights incoming to neuron i from a given set of neurons
# i: neuron index
# set: the set of presynaptic neurons
def latePhaseWeightsFromSet(i, set):
	global adj
	global z
	zi = z[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(zi)), set))
	inc_z = zi[adj_connections]

	return inc_z

# latePhaseWeightsFromCore
# Prints and returns all the late-phase synaptic weights incoming to neuron i from core neurons
# i: neuron index
# pr [optional]: specifies if result is printed
def latePhaseWeightsFromCore(i, pr = True):
	inc_z = latePhaseWeightsFromSet(i, core)
	if pr:
		print("Incoming late-phase weights from core to neuron " + str(i) + " (" + str(len(inc_z)) + "): \n" + str(inc_z))

	return inc_z

# latePhaseWeightsToCore
# Prints and returns all the late-phase synaptic weights outgoing to core neurons from neuron i
# i: neuron index
# pr [optional]: specifies if result is printed
def latePhaseWeightsToCore(i, pr = True):
	global adj
	global z
	zi = z[i, :]
	adj_connections = np.logical_and(adj[i, :] == 1, np.in1d(np.arange(len(zi)), core))
	out_z = zi[adj_connections]
	if pr:
		print("Outgoing late-phase weights from neuron " + str(i) + " to core (" + str(len(out_z)) + "): \n" + str(out_z))

	return out_z

# meanEarlyPhaseWeight
# Returns the mean early-phase synaptic weight between two sets of neurons; prints the connectivity
# and the mean early-phase weight
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
def meanEarlyPhaseWeight(set, set2 = [], pr = True):
	summed_weight = 0
	connection_num = 0
	self_set = False

	if len(set2) <= 0:
		set2 = set # consider internal connections within "set" if "set2" is not specified
		self_set = True

	for n in set2:
		inc_weights = earlyPhaseWeightsFromSet(n, set)
		summed_weight += np.sum(inc_weights)
		connection_num += len(inc_weights)

	ret = summed_weight / connection_num

	if pr:
		if self_set:
			print("Self-connectivity: " + str(connection_num / (len(set)**2 - len(set))))
		else:
			print("Connectivity: " + str(connection_num / (len(set)*len(set2))))

		print("Mean early-phase weight: " + str(ret))

	return ret

# sdEarlyPhaseWeight
# Prints and returns the standard deviation of the early-phase synaptic weights between two sets of neurons
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
def sdEarlyPhaseWeight(set, set2 = [], pr = True):
	mean = meanEarlyPhaseWeight(set, set2, False)
	summed_qu_dev = 0
	connection_num = 0

	if len(set2) <= 0:
		set2 = set # consider internal connections within "set" if "set2" is not specified

	for n in set2:
		qu_devs = np.power(mean - earlyPhaseWeightsFromSet(n, set), 2) # quadratic deviations from mean
		summed_qu_dev += np.sum(qu_devs)
		connection_num += len(qu_devs)

	ret = np.sqrt(summed_qu_dev / (connection_num-1))

	if pr:
		print("Std. dev. of early-phase weight: " + str(ret))

	return ret

# meanLatePhaseWeight
# Returns the mean late-phase synaptic weight between two sets of neurons; prints the connectivity
# and the mean late-phase weight
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
def meanLatePhaseWeight(set, set2 = [], pr = True):
	summed_weight = 0
	connection_num = 0
	self_set = False

	if len(set2) <= 0:
		set2 = set # consider internal connections within "set" if "set2" is not specified
		self_set = True

	for n in set2:
		inc_weights = latePhaseWeightsFromSet(n, set)
		summed_weight += np.sum(inc_weights)
		connection_num += len(inc_weights)

	ret = summed_weight / connection_num

	if pr:
		if self_set:
			print("Self-connectivity: " + str(connection_num / (len(set)**2 - len(set))))
		else:
			print("Connectivity: " + str(connection_num / (len(set)*len(set2))))

		print("Mean late-phase weight: " + str(ret))

	return ret

# sdLatePhaseWeight
# Prints and returns the standard deviation of the late-phase synaptic weights between two sets of neurons
# set: the first set of neurons (presynaptic)
# set2 [optional]: the seconds set of neurons (postsynaptic); if not specified, connections within "set" are considered
# pr [optional]: specifies if result shall be printed
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

	ret = np.sqrt(summed_qu_dev / (connection_num-1))

	if pr:
		print("Std. dev. of late-phase weight: " + str(ret))

	return ret


# connectionsToCore
# Prints and returns all the core neurons to which neuron i provides input
# i: neuron index
# pr [optional]: specifies if result is printed
def connectionsToCore(i, pr = True):
	global adj
	ctc = core[np.in1d(core, outgoingConnections(i, False))].flatten()
	if pr:
		print("Connections from neuron " + str(i) + " to core (" + str(len(ctc)) + "): \n" + str(ctc))

	return ctc

# connectionsToRSCore
# Prints and returns all the recall-stimulated core neurons to which neuron i provides input
# i: neuron index
# pr [optional]: specifies if result is printed
def connectionsToRSCore(i, pr = True):
	global adj
	ctrsc = core_recall[np.in1d(core_recall, outgoingConnections(i, False))].flatten()
	if pr:
		print("Connections from neuron " + str(i) + " to recall-stim. core (" + str(len(ctrsc)) + "): \n" + str(ctrsc))

	return ctrsc

# connectionsToNonCore
# Prints and returns all the non-core neurons to which neuron i provides input
# i: neuron index
# pr [optional]: specifies if result is printed
def connectionsToNonCore(i, pr = True):
	global adj
	outg = outgoingConnections(i, False)
	ctnc = outg[np.negative(np.in1d(outg, core))]
	if pr:
		print("Connections from neuron " + str(i) + " to non-core neurons (" + str(len(ctnc)) + "): \n" + str(ctnc))

	return ctnc

# setRhombCore
# Sets a rhomb-shaped core depending on given center and radius
# core_center: the central neuron of the rhomb
# core_radius: the "radius" of the rhomb
def setRhombCore(core_center, core_radius):
	global core

	core = np.array([], dtype=np.int32)
	core_size = 2*core_radius**2 + 2*core_radius + 1

	for i in range(-core_radius, core_radius+1, 1):
		num_cols = (core_radius-abs(i))

		for j in range(-num_cols, num_cols+1, 1):
			core = np.append(core, np.array([core_center+i*Nl+j]))

# printMeanWeights
# Prints mean and standard deviation of CA, outgoing, incoming, and control weight
def printMeanWeights():
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
# argv[]: timestamps of two simulations

if __name__ == "__main__":

	if len(sys.argv) < 3:
		print("Not enough arguments provided!")
		exit()
	else:
		ts1 = str(sys.argv[1]) # timestamp for simulation data before consolidation
		ts2 = str(sys.argv[2]) # timestamp for simulation data after consolidation

	print("##############################################")
	print("Before 10s-recall:")
	loadWeightMatrix(ts1 + "_net_20.0.txt")
	printMeanWeights()

	print("##############################################")
	print("After 10s-recall:")
	loadWeightMatrix(ts1 + "_net_20.1.txt")
	printMeanWeights()

	print("##############################################")
	print("Before 8h-recall:")
	loadWeightMatrix(ts2 + "_net_28810.0.txt")
	printMeanWeights()

	print("##############################################")
	print("After 8h-recall:")
	loadWeightMatrix(ts2 + "_net_28810.1.txt")
	printMeanWeights()
