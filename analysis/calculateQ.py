##########################################################################################
### Functions to calculate the quality coefficient Q from the firing rate distribution ###
###                   during recall of an input-defined cell assembly                  ###
##########################################################################################

### Copyright 2018-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
from pathlib import Path
from subprocess import call

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays

# calculateQ
# From given mean firing rates of subpopulations, calculates the Q value 
# that measures how well recall of an input-defined pattern works
# v_as: mean firing rate of the neuronal subpopulation that receives learning and recall stimulation
# v_as_err: uncertainty of v_as
# v_ans: mean firing rate of the neuronal subpopulation that receives learning but no recall stimulation
# v_ans_err: uncertainty of v_ans
# v_ctrl: mean firing rate of the neuronal subpopulation that receives neither learning nor recall stimulation
# v_ctrl_err: uncertainty of v_ctrl
def calculateQ(v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err):

	Q = (v_ans - v_ctrl) / v_as

	Q_err = np.sqrt( np.power( v_ans_err / v_as, 2 ) + np.power( v_ctrl_err / v_as, 2 ) \
	                 + np.power (v_as_err * (v_ans - v_ctrl) / np.power(v_as, 2), 2 ) ) # error propagated from v_as, v_ans, v_ctrl

	return Q, Q_err

# getSubpopulations
# Retrieves the indices of the neurons belonging to the different network subpopulations
# N_pop: the number of neurons in the whole population
# core: array of the neurons belonging to the stimulated core
# p_r: the fraction of core neurons that are stimulated for recall
# return: arrays of neuron indices of network subpopulations that receive 1. learning and recall stimulation, 2. learning but no recall stimulation, 3. neither learning nor recall stimulation
def getSubpopulations(N_pop, core, p_r):

	core_recall = core[0:int(np.floor(float(p_r)*core.shape[0]))]
	core_norecall = core[np.logical_not(np.in1d(core, core_recall))]
	control = np.delete(np.arange(N_pop), core)

	return core_recall, core_norecall, control

# calculateMeanRatesAndQ
# Calculates the mean firing rates of subpopulations from given firing rate distribution across a network,
# then calculates the Q value that measures how well recall of an input-defined pattern works
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# core: array of the neurons belonging to the stimulated core
# Nl_exc: the number of excitatory neurons in one line of a quadratic grid
# time_for_activity: the time that at which the activities shall be read out (some time during recall)
# recall_fraction: the fraction of core neurons that are activated for recall
def calculateMeanRatesAndQ(nppath, timestamp, core, Nl_exc, time_for_activity, recall_fraction):

	path = ""

	core_recall, core_norecall, control = getSubpopulations(Nl_exc*Nl_exc, core, recall_fraction)

	v_array = np.zeros(Nl_exc*Nl_exc) # data array

	# look for data file [timestamp]_net_[time_for_activity].txt
	path = ""
	rawpaths = Path(nppath)

	for x in rawpaths.iterdir():

		tmppath = str(x)

		if (timestamp + "_net_" + time_for_activity + ".txt") in tmppath:
			path = tmppath

	if path == "":
		raise FileNotFoundError('"' + timestamp + '_net_' + time_for_activity + '.txt" was not found')

	try:
		f = open(path)

	except OSError:
		raise

	# read activities from file and determine mean activities for different regions
	rawdata = f.read()
	rawdata = rawdata.split('\n')
	nn = len(rawdata)-1
	f.close()

	if nn != 2*Nl_exc*Nl_exc+Nl_exc+3:
		raise ValueError(str(nn) + ' instead of ' + str(2*Nl_exc*Nl_exc+Nl_exc+3) + ' lines in data file "' + path + '"')

	offset = 2*Nl_exc*Nl_exc+2

	for n in range(nn-1): # loop over lines

		if n >= offset:
			n2 = (n - offset) * Nl_exc
			line_values = rawdata[n].split("\t\t")

			for p in range(len(line_values)):
				v_array[n2+p] = float(line_values[p])

	v_as = np.sum(v_array[core_recall]) / len(core_recall)
	v_ans = np.sum(v_array[core_norecall]) / len(core_norecall)
	v_ctrl = np.sum(v_array[control]) / len(control)

	v_as_err = np.sqrt(np.sum(np.power(v_array[core_recall], 2)) / len(core_recall) - np.power(v_as, 2)) # standard deviation
	v_ans_err = np.sqrt(np.sum(np.power(v_array[core_norecall], 2)) / len(core_norecall) - np.power(v_ans, 2))
	v_ctrl_err = np.sqrt(np.sum(np.power(v_array[control], 2)) / len(control) - np.power(v_ctrl, 2))

	Q, Q_err = calculateQ(v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err)

	return Q, Q_err, v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err
