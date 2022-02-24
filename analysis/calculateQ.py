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
# Calculates Q measure for how well recall works
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# core: array of the neurons belonging to the stimulated core
# Nl_exc: the number of excitatory neurons in one line of a quadratic grid
# time_for_activity: the time that at whcih the activites shall be read out (some time during recall)
# recall_fraction: the fraction of core neurons that are activated for recall
def calculateQ(nppath, timestamp, core, Nl_exc, time_for_activity, recall_fraction):

	core_recall = core[0:int(np.floor(float(recall_fraction)*core.shape[0]))]
	core_norecall = core[np.logical_not(np.in1d(core, core_recall))]
	control = np.delete(np.arange(Nl_exc*Nl_exc), core)

	path = ""

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

	Q = (v_ans - v_ctrl) / v_as

	Q_err = np.sqrt( np.power( v_ans_err / v_as, 2 ) + np.power( v_ctrl_err / v_as, 2 ) \
	                 + np.power (v_as_err * (v_ans - v_ctrl) / np.power(v_as, 2), 2 ) ) # error propagated from v_as, v_ans, v_ctrl

	return Q, Q_err, v_as, v_as_err, v_ans, v_ans_err, v_ctrl, v_ctrl_err
