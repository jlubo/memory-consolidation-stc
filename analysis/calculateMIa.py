#######################################################################################
### Functions to calculate the mutual information from the firing rate distribution ###
###                    in a network at two different times                          ###
#######################################################################################

### Copyright 2018-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

from readWeightData import *
from probDistributions import *
import numpy as np
from pathlib import Path

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays

# calculateMIa
# Calculates mutual information of the activity distribution of the network at two timesteps and returns it along with the self-information of
# the reference distribution
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# Nl: the number of excitatory neurons in one line of a quadratic grid
# time_for_activity: the time that at which the activites shall be read out (some time during recall)
# time_ref: the reference time (for getting the activity distribution during learning)
# core: the neurons in the cell assembly (for stipulation; only required if no activity distribution during learning is available)
def calculateMIa(nppath, timestamp, Nl, time_for_activity, time_ref = "11.0", core = np.array([])):

	if time_ref: # use reference firing rate distribution from data (for learned cell assembly)
		times_for_readout_list = [time_ref, time_for_activity] # the simulation times at which the activities shall be read out
	else: # use model firing rate distribution (for stipulated cell assembly)
		times_for_readout_list = [time_for_activity]
		v_model = np.zeros(Nl**2)
		v_model[core] = 1 # entropy/MI of this distribution for Nl=40: 0.4689956
		v_model = np.reshape(v_model, (Nl,Nl))

	connections = [np.zeros((Nl**2,Nl**2)) for x in times_for_readout_list]
	h = [np.zeros((Nl**2,Nl**2)) for x in times_for_readout_list]
	z = [np.zeros((Nl**2,Nl**2)) for x in times_for_readout_list]
	v = [np.zeros((Nl,Nl)) for x in times_for_readout_list]

	v_array = np.zeros(Nl*Nl) # data array

	rawpaths = Path(nppath)

	for i in range(len(times_for_readout_list)):

		time_for_readout = times_for_readout_list[i]
		path = ""

		# look for data file [timestamp]_net_[time_for_readout].txt
		for x in rawpaths.iterdir():

			tmppath = str(x)

			if (timestamp + "_net_" + time_for_readout + ".txt") in tmppath:
				path = tmppath

		if path == "":
			raise FileNotFoundError('"' + timestamp + '_net_' + time_for_readout + '.txt" was not found')

		try:
			connections[i], h[i], z[i], v[i] = readWeightMatrixData(path, Nl)

		except ValueError:
			raise
		except OSError:
			raise

	### Activity ###
	if time_ref: # use reference firing rate distribution from data (for learned cell assembly)
		margEntropyActL = marginalEntropy(v[0])
		print("margEntropyActL = " + str(margEntropyActL))

		margEntropyActR = marginalEntropy(v[1])
		print("margEntropyActR = " + str(margEntropyActR))

		jointEntropyAct = jointEntropy(v[0],v[1])
		print("jointEntropyAct = " + str(jointEntropyAct))

	else: # use model firing rate distribution (for stipulated cell assembly)
		margEntropyActL = marginalEntropy(v_model)
		print("margEntropyActL = " + str(margEntropyActL))

		margEntropyActR = marginalEntropy(v[0])
		print("margEntropyActR = " + str(margEntropyActR))

		jointEntropyAct = jointEntropy(v_model,v[0])
		print("jointEntropyAct = " + str(jointEntropyAct))

	### Results and Output ###
	MIa = mutualInformation(margEntropyActL, margEntropyActR, jointEntropyAct)
	print("MIa = " + str(MIa))

	return MIa, margEntropyActL

# marginalEntropy
# Computes the marginal entropy of an array
# a: array of outcomes (e.g., array of activities of all neurons in a network or array of weights of all synapses in a network)
# return: the Shannon entropy of a
def marginalEntropy(a):

	val, p = marginalProbDist(a)
	p = p[p > 0]
	entropy = - np.sum(p * np.log2(p))

	return entropy


# jointEntropy
# Computes the joint entropy of two arrays
# a1: first array of outcomes (e.g., array of activities of all neurons in a network or array of weights of all synapses in a network)
# a2: second array of outcomes
# return: the joint Shannon entropy of (a1, a2)
def jointEntropy(a1, a2):

	val, p = jointProbDist(a1, a2)
	p = p[p > 0]
	entropy = - np.sum(p * np.log2(p))

	return entropy

# mutualInformation
# Computes the mutual information of two arrays
# margEntropy1: marginal entropy of the first array of outcomes
# margEntropy2: marginal entropy of the second array of outcomes
# jEntropy: joint entropy of both arrays of outcomes
# return: the mutual information of (a1, a2)
def mutualInformation(margEntropy1, margEntropy2, jEntropy):

	#margEntropy1 = marginalEntropy(a1)
	#margEntropy2 = marginalEntropy(a2)
	#jEntropy = jointEntropy(a1, a2)

	MI = margEntropy1 + margEntropy2 - jEntropy

	return MI
