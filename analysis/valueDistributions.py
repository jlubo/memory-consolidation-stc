############################################################################################
### Functions to analyze and plot weight and activity distributions from simulation data ###
############################################################################################

### Copyright 2019-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

from utilityFunctions import *
import sys
import warnings
import os.path
import numpy as np
from pathlib import Path
from subprocess import call

# findOverallMinMax
# Determines the minimum and maximum values across all data files that are located somewhere in a given directory
# and that have the same readout time
# nppath: path to the directory to read the data from
# Nl_exc: the number of excitatory neurons in one line of a quadratic grid
# time_for_readout: the time that at which the weights shall be read out (as a string)
# h_0: the initial weight, and normalization factor for z
# return: two-dimensional array containing the minimum and maximum values for the four different data types
def findOverallMinMax(nppath, Nl_exc, time_for_readout, h_0):

	sysmin, sysmax = sys.float_info.min, sys.float_info.max
	(h_min, z_min, w_min, v_min) = (sysmax, sysmax, sysmax, sysmax) # initially, assign the maximum possible value
	(h_max, z_max, w_max, v_max) = (sysmin, sysmin, sysmin, sysmin) # initially, assign the minimum possible value

	# recurseFind
	# Function to recursively move through directories and look for data to find their minima/maxima
	# path: the directory to iterate through
	def recurseFindMinMax(path):
		nonlocal h_min, z_min, w_min, v_min
		nonlocal h_max, z_max, w_max, v_max

		rawpaths = Path(path)
		for x in rawpaths.iterdir():
			if x.is_dir():
				recurseFindMinMax(x) # if the found file is a directory, recurse into it

			tmppath = str(x)
			if ("_net_" + time_for_readout + ".txt") in tmppath: # file containing network simulation data found

				# read data from file
				try:
					connections, h, z, v = readWeightMatrixData(tmppath, Nl_exc)
					h = h[connections] # reduce h (leave out non-existent synapses)
					z = h_0*z[connections] # reduce and normalize z
					w = h + z # compute total synaptic weight
				except ValueError:
					raise
				except OSError:
					raise

				# checkAndAdjust
				# Compares two numbers and returns the larger/lower one, depending on the operator
				# a: a floating point number
				# b: a floating point number
				# op [optional]: the operator to be used
				# return: the larger/lower one of the two numbers
				def checkAndAdjust(a, b, op=">"):
					if b > a:
						return b if op == ">" else a
					else:
						return a if op == ">" else b

				# adjust maxima
				h_max = checkAndAdjust(h_max, np.max(h), ">")
				z_max = checkAndAdjust(z_max, np.max(z), ">")
				w_max = checkAndAdjust(w_max, np.max(w), ">")
				v_max = checkAndAdjust(v_max, np.max(v), ">")

				# adjust minima
				h_min = checkAndAdjust(h_min, np.min(h), "<")
				z_min = checkAndAdjust(z_min, np.min(z), "<")
				w_min = checkAndAdjust(w_min, np.min(w), "<")
				v_min = checkAndAdjust(v_min, np.min(v), "<")

	# iterate across files in the directory
	recurseFindMinMax(nppath)

	return np.array([[h_min, h_max], [z_min, z_max], [w_min, w_max], [v_min, v_max]])

# plotDistributions
# Creates data and plot files of the weight and activity distribution at a given time
# nppath: path to the directory to read the data from
# timestamp: a string containing date and time (to access correct paths) OR equal to "any"
# add: additional descriptor
# Nl_exc: the number of excitatory neurons in one line of a quadratic grid
# h_0: the initial weight, and normalization factor for z
# time_for_readout: the time that at which the weights shall be read out (as a string)
# core: array of indices of the cell assembly (core) neurons
# norm_all [optional]: specifies whether to normalize across all subpopulations (True) or across each subpop. individually (False)
#                      - the first is recommendable if samples of different subpopulations are compared against each other,
#                        the latter is recommendable if different samples of the same subpopulation are compared
# bins [optional]: list of four arrays, each containing the bins for one of the four quantities
def plotDistributions(nppath, timestamp, add, Nl_exc, h_0, time_for_readout, core, norm_all=False, bins=None):

	orgdir = os.getcwd() # store the current working directory

	 # "any" case: not looking for a specific timestamp, but for any data with a certain time_for_readout in the given directory
	if timestamp == "any":
		if bins is None:
			warnings.warn("Warning: timestamp=\"any\": bins should be provided by the calling function to compare across trials.")
		rawpaths = Path(nppath)
		for x in rawpaths.iterdir():
			tmppath = os.path.split(str(x))[1] # remove head from path
			if ("_net_" + time_for_readout + ".txt") in tmppath:
				timestamp = tmppath.split("_net_")[0]
				plotDistributions(nppath, timestamp, add, Nl_exc, h_0, time_for_readout, core, norm_all, bins) # call this function again, now with specific timestamp
		return

	# read data from file [timestamp]_net_[time_for_readout].txt
	os.chdir(nppath) # change to data directory
	try:
		connections, h, z, v = readWeightMatrixData(timestamp + "_net_" + time_for_readout + ".txt", Nl_exc)
		z = h_0*z # normalize z
		w = h + z # compute total synaptic weight
	except ValueError:
		raise
	except OSError:
		raise

	# determine subpopulations
	N_tot = Nl_exc**2 # total number of neurons
	N_CA = len(core) # number of neurons in the cell assembly
	N_control = N_tot - N_CA # number of neurons in the control subpopulation
	all = np.arange(N_tot)
	noncore = all[np.logical_not(np.in1d(all, core))] # array of indices of the neurons not in the cell assembly (core)

	block_CA_within = np.ones((N_CA, N_CA), dtype=bool) # array of ones for the synapses within the cell assembly
	block_CA_outgoing = np.ones((N_CA, N_control), dtype=bool) # array of ones for the synapses outgoing from the cell assembly
	block_CA_incoming = np.ones((N_control, N_CA), dtype=bool) # array of ones for the synapses incoming to the cell assembly
	block_control = np.ones((N_control, N_control), dtype=bool) # array of ones for the synapses within the control subpopulation

	mask_CA_within = np.append(np.append(block_CA_within, np.logical_not(block_CA_outgoing), axis=1), \
	                           np.logical_not(np.append(block_CA_incoming, block_control, axis=1)),
	                           axis=0) # binary mask defining the synapses within the cell assembly
	mask_CA_outgoing = np.append(np.append(np.logical_not(block_CA_within), block_CA_outgoing, axis=1), \
	                             np.logical_not(np.append(block_CA_incoming, block_control, axis=1)),
	                             axis=0) # binary mask defining the synapses outgoing from the cell assembly
	mask_CA_incoming = np.append(np.logical_not(np.append(block_CA_within, block_CA_outgoing, axis=1)), \
	                             np.append(block_CA_incoming, np.logical_not(block_control), axis=1),
	                             axis=0) # binary mask defining the synapses incoming to the cell assembly
	mask_control = np.append(np.logical_not(np.append(block_CA_within, block_CA_outgoing, axis=1)), \
	                         np.append(np.logical_not(block_CA_incoming), block_control, axis=1),
	                         axis=0) # binary mask defining the synapses within the control subpopulation

	# early-phase weights
	'''h_CA_within = h[mask_CA_within]
	h_CA_outgoing = h[mask_CA_outgoing]
	h_CA_incoming = h[mask_CA_incoming]
	h_control = h[mask_control]'''
	h_CA_within = h[np.logical_and(connections, mask_CA_within)]
	h_CA_outgoing = h[np.logical_and(connections, mask_CA_outgoing)]
	h_CA_incoming = h[np.logical_and(connections, mask_CA_incoming)]
	h_control = h[np.logical_and(connections, mask_control)]

	# late-phase weights
	'''z_CA_within = z[mask_CA_within]
	z_CA_outgoing = z[mask_CA_outgoing]
	z_CA_incoming = z[mask_CA_incoming]
	z_control = z[mask_control]'''
	z_CA_within = z[np.logical_and(connections, mask_CA_within)]
	z_CA_outgoing = z[np.logical_and(connections, mask_CA_outgoing)]
	z_CA_incoming = z[np.logical_and(connections, mask_CA_incoming)]
	z_control = z[np.logical_and(connections, mask_control)]

	# total synaptic weights
	w_CA_within = h_CA_within + z_CA_within
	w_CA_outgoing = h_CA_outgoing + z_CA_outgoing
	w_CA_incoming = h_CA_incoming + z_CA_incoming
	w_control = h_control + z_control

	# firing rates
	v_CA = v.flatten()[np.in1d(all, core)]
	v_control = v.flatten()[np.logical_not(np.in1d(all, core))]

	# discretization of the distribution
	if bins is None:
		binh = np.linspace(np.min(h), np.max(h), 101, endpoint=True) # create range of bins for marginalProbDist(h...)
		binz = np.linspace(np.min(z), np.max(z), 101, endpoint=True) # create range of bins for marginalProbDist(z...)
		binw = np.linspace(np.min(w), np.max(w), 101, endpoint=True) # create range of bins for marginalProbDist(w...)
		binv = np.linspace(np.min(v), np.max(v), 101, endpoint=True) # create range of bins for marginalProbDist(v...)
	else:
		[binh, binz, binw, binv] = bins # use pre-defined bins

	hstep = binh[1]-binh[0]
	zstep = binz[1]-binz[0]
	wstep = binw[1]-binw[0]
	vstep = binv[1]-binv[0]
	valh = np.delete(binh, -1) + hstep/2 # use mean values instead of lower bounds of the bins as values
	valz = np.delete(binz, -1) + zstep/2 # use mean values instead of lower bounds of the bins as values
	valw = np.delete(binw, -1) + wstep/2 # use mean values instead of lower bounds of the bins as values
	valv = np.delete(binv, -1) + vstep/2 # use mean values instead of lower bounds of the bins as values

	# normalization of the distribution
	if norm_all:
		norm_value_w = np.sum(connections) # normalization factor for weights (number of all connections)
		norm_value_v = N_CA + N_control # normalization factor for activities (number of all neurons)
	else:
		norm_value_w = None # use default (normalization across each subpopulation individually)
		norm_value_v = None # use default (normalization across CA and control individually)

	buf, ph_CA_within = marginalProbDist(h_CA_within, binning = True, bin_edges = binh, norm = norm_value_w)
	buf, ph_CA_outgoing = marginalProbDist(h_CA_outgoing, binning = True, bin_edges = binh, norm = norm_value_w)
	buf, ph_CA_incoming = marginalProbDist(h_CA_incoming, binning = True, bin_edges = binh, norm = norm_value_w)
	buf, ph_control = marginalProbDist(h_control, binning = True, bin_edges = binh, norm = norm_value_w)
	buf, pz_CA_within = marginalProbDist(z_CA_within, binning = True, bin_edges = binz, norm = norm_value_w)
	buf, pz_CA_outgoing = marginalProbDist(z_CA_outgoing, binning = True, bin_edges = binz, norm = norm_value_w)
	buf, pz_CA_incoming = marginalProbDist(z_CA_incoming, binning = True, bin_edges = binz, norm = norm_value_w)
	buf, pz_control = marginalProbDist(z_control, binning = True, bin_edges = binz, norm = norm_value_w)
	buf, pw_CA_within = marginalProbDist(w_CA_within, binning = True, bin_edges = binw, norm = norm_value_w)
	buf, pw_CA_outgoing = marginalProbDist(w_CA_outgoing, binning = True, bin_edges = binw, norm = norm_value_w)
	buf, pw_CA_incoming = marginalProbDist(w_CA_incoming, binning = True, bin_edges = binw, norm = norm_value_w)
	buf, pw_control = marginalProbDist(w_control, binning = True, bin_edges = binw, norm = norm_value_w)
	buf, pv_CA = marginalProbDist(v_CA, binning = True, bin_edges = binv, norm = norm_value_v)
	buf, pv_control = marginalProbDist(v_control, binning = True, bin_edges = binv, norm = norm_value_v)

	# write early-phase weight distribution to file
	f = open(timestamp + "_eweight_dist_" + time_for_readout + add + ".txt", "w")
	for i in range(len(valh)):
		f.write(str(valh[i]) + "\t\t" + str(ph_CA_within[i]) + "\t\t" + str(ph_CA_outgoing[i]) + "\t\t" + \
		        str(ph_CA_incoming[i]) + "\t\t" + str(ph_control[i]) + "\n")
	f.close()

	# write late-phase weight distribution to file
	f = open(timestamp + "_lweight_dist_" + time_for_readout + add + ".txt", "w")
	for i in range(len(valz)):
		f.write(str(valz[i]) + "\t\t" + str(pz_CA_within[i]) + "\t\t" + str(pz_CA_outgoing[i]) + "\t\t" + \
		        str(pz_CA_incoming[i]) + "\t\t" + str(pz_control[i]) + "\n")
	f.close()

	# write distribution of total synaptic weights to file
	f = open(timestamp + "_totweight_dist_" + time_for_readout + add + ".txt", "w")
	for i in range(len(valw)):
		f.write(str(valw[i]) + "\t\t" + str(pw_CA_within[i]) + "\t\t" + str(pw_CA_outgoing[i]) + "\t\t" + \
		        str(pw_CA_incoming[i]) + "\t\t" + str(pw_control[i]) + "\n")
	f.close()

	# write activity distribution to file
	f = open(timestamp + "_act_dist_" + time_for_readout + add + ".txt", "w")
	for i in range(len(valv)):
		f.write(str(valv[i]) + "\t\t" + str(pv_CA[i]) + "\t\t" + str(pv_control[i]) + "\n")
	f.close()

	# write gnuplot script
	f = open(timestamp + "_plot_dist.gpl", "w")
	f.write("### DO NOT EDIT THIS FILE! IT WILL BE OVERWRITTEN. ###\n\n" + \
	        "#set terminal png size 1024,640 enhanced\nset terminal pdf enhanced\n\n" + \
	        "#set style fill transparent solid 0.8 noborder # for 'boxes' style\n" + \
			"#set style fill transparent pattern 4 bo # for 'boxes' style\n" + \
	        "set log y\nset format y \"%.0e\"\nset yrange [3e-06:1]\nset key outside\n\n" + \
			"h_0 = " + str(h_0) + "\n" +\
			"epsilon = " + str(epsilon) + "\n\n")

	# plotting of early-phase weight distribution
	f.write("set output \"" + timestamp + "_eweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binh[0]-10*hstep) + "/h_0:" + str(binh[-1]+10*hstep) + "/h_0]\n")
	f.write("set xlabel \"Early-phase weight / h_0\"\nset ylabel \"Relative frequency\"\n")
	f.write("plot \"" + timestamp + "_eweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0):($2 > 0 ? $2 : epsilon) t \"CA\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($3 > 0 ? $3 : epsilon) t \"outgoing\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($4 > 0 ? $4 : epsilon) t \"incoming\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($5 > 0 ? $5 : epsilon) t \"control\" with histeps\n")

	# plotting of late-phase weight distribution
	f.write("\nset output \"" + timestamp + "_lweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binz[0]-10*zstep) + "/h_0:" + str(binz[-1]+10*zstep) + "/h_0]\n")
	f.write("set xlabel \"Late-phase weight / h_0\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_lweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0):($2 > 0 ? $2 : epsilon) t \"CA\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($3 > 0 ? $3 : epsilon) t \"outgoing\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($4 > 0 ? $4 : epsilon) t \"incoming\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0):($5 > 0 ? $5 : epsilon) t \"control\" with histeps\n")

	# plotting of total weight distribution
	f.write("\nset output \"" + timestamp + "_totweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binw[0]-10*wstep) + "/h_0*100:" + str(binw[-1]+10*wstep) + "/h_0*100]\n")
	f.write("set xlabel \"Total synaptic weight (%)\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_totweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0*100):($2 > 0 ? $2 : epsilon) t \"CA\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0*100):($3 > 0 ? $3 : epsilon) t \"outgoing\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0*100):($4 > 0 ? $4 : epsilon) t \"incoming\" with histeps, \\\n" + \
	        "     \"\" using ($1/h_0*100):($5 > 0 ? $5 : epsilon) t \"control\" with histeps\n")

	# plotting of activity distribution
	f.write("\nset output \"" + timestamp + "_act_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binv[0]-10*vstep) + ":" + str(binv[-1]+10*vstep) + "]\n")
	f.write("set xlabel \"Neuronal firing rate (Hz)\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_act_dist_" + time_for_readout + add + ".txt\" using 1:($2 > 0 ? $2 : epsilon) t \"CA\" with histeps, \\\n" + \
	        "     \"\" using 1:($3 > 0 ? $3 : epsilon) t \"control\" with histeps\n\n")

	f.close()

	call(["gnuplot", timestamp + "_plot_dist.gpl"])
	os.chdir(orgdir) # change back to original directory

# plotWeightDistributions3CAs
# Creates data and plot files of the weight distribution of a network with 2 or 3, possibly overlapping, assemblies at a given time
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# add: additional descriptor
# Nl_exc: the number of excitatory neurons in one line of a quadratic grid
# h_0: the initial weight, and normalization factor for z
# time_for_readout: the time that at which the weights shall be read out
# coreA: array of indices of the first cell assembly (core) neurons
# coreB [optional]: array of indices of the second cell assembly (core) neurons
# coreC [optional]: array of indices of the third cell assembly (core) neurons
# bins [optional]: list of three arrays, each containing the bins for one of the four quantities
def plotWeightDistributions3CAs(nppath, timestamp, add, Nl_exc, h_0, time_for_readout, coreA, coreB = None, coreC = None, bins = None):

	orgdir = os.getcwd() # store the current working directory

	# "any" case: not looking for a specific timestamp, but for any data with a certain time_for_readout in the given directory
	if timestamp == "any":

		rawpaths = Path(nppath)
		for x in rawpaths.iterdir():
			tmppath = os.path.split(str(x))[1] # remove head from path
			if ("_net_" + time_for_readout + ".txt") in tmppath:
				timestamp = tmppath.split("_net_")[0]
				plotWeightDistributions3CAs(nppath, timestamp, add, Nl_exc, h_0, time_for_readout, coreA, coreB, coreC, bins) # call this function again, now with specific timestamp; bins should be provided by calling function
		return

	# read data from file [timestamp]_net_[time_for_readout].txt
	os.chdir(nppath) # change to data directory
	try:
		connections, h, z, v = readWeightMatrixData(timestamp + "_net_" + time_for_readout + ".txt", Nl_exc)
		z = h_0*z # normalize z
		w = h + z # compute total synaptic weight
	except ValueError:
		raise
	except OSError:
		raise

	# determine synapses within the cell assemblies
	N_tot = Nl_exc**2 # total number of neurons

	mask_coreA = np.zeros((N_tot, N_tot), dtype=bool)
	for syn_pre in coreA:
		for syn_post in coreA:
			if connections[syn_pre,syn_post]: # NEW, TEST
				mask_coreA[syn_pre,syn_post] = True

	mask_coreB = np.zeros((N_tot, N_tot), dtype=bool)
	if coreB is not None:
		for syn_pre in coreB:
			for syn_post in coreB:
				if connections[syn_pre,syn_post]: # NEW, TEST
					mask_coreB[syn_pre,syn_post] = True

	mask_coreC = np.zeros((N_tot, N_tot), dtype=bool)
	if coreC is not None:
		for syn_pre in coreC:
			for syn_post in coreC:
				if connections[syn_pre,syn_post]: # NEW, TEST
					mask_coreC[syn_pre,syn_post] = True

	# find control synapses (all synapses that are not within a cell assembly)
	mask_control = np.logical_and(connections, np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC))))

	# synapses outgoing from A // TODO
	#block_outgoing_A = np.ones((len(coreA), N_tot-len(coreA)), dtype=bool) # array of ones for the synapses outgoing from assembly A
	#mask_A_to_B = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_A_to_C = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_A_to_ctrl = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))

	# synapses outgoing from B // TODO
	#mask_B_to_A = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_B_to_C = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_B_to_ctrl = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))

	# synapses outgoing from C // TODO
	#mask_C_to_A = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_C_to_B = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))
	#mask_C_to_ctrl = np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))

	# synapses incoming... // TODO

	# find exclusive intersections
	mask_I_AB = np.logical_and( np.logical_and(mask_coreA, mask_coreB), np.logical_not(mask_coreC) )
	mask_I_AC = np.logical_and( np.logical_and(mask_coreA, mask_coreC), np.logical_not(mask_coreB) )
	mask_I_BC = np.logical_and( np.logical_and(mask_coreB, mask_coreC), np.logical_not(mask_coreA) )
	mask_I_ABC = np.logical_and( mask_coreA, np.logical_and(mask_coreB, mask_coreC) )

	# remove intersections from exclusive cores
	mask_coreA = np.logical_and(mask_coreA, \
	                            np.logical_and(np.logical_not(mask_I_AB), \
				                               np.logical_and(np.logical_not(mask_I_AC), np.logical_not(mask_I_ABC))))
	mask_coreB = np.logical_and(mask_coreB, \
	                            np.logical_and(np.logical_not(mask_I_AB), \
				                               np.logical_and(np.logical_not(mask_I_BC), np.logical_not(mask_I_ABC))))
	mask_coreC = np.logical_and(mask_coreC, \
	                            np.logical_and(np.logical_not(mask_I_AC), \
				                               np.logical_and(np.logical_not(mask_I_BC), np.logical_not(mask_I_ABC))))

	# tests (each should yield true)
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_coreB)))
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_coreC)))
	#print("Test:", not np.any(np.logical_and(mask_coreB, mask_coreC)))
	#print("Test:", not np.any(np.logical_and(mask_I_AB, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_I_AB, mask_I_AC)))
	#print("Test:", not np.any(np.logical_and(mask_I_AB, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_I_AC, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_I_AC, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_I_BC, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_coreA)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_coreB)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_coreC)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_I_AB)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_I_AC)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_control, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_I_AB)))
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_I_AC)))
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_coreA, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_coreB, mask_I_AB)))
	#print("Test:", not np.any(np.logical_and(mask_coreB, mask_I_AC)))
	#print("Test:", not np.any(np.logical_and(mask_coreB, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_coreB, mask_I_ABC)))
	#print("Test:", not np.any(np.logical_and(mask_coreC, mask_I_AB)))
	#print("Test:", not np.any(np.logical_and(mask_coreC, mask_I_AC)))
	#print("Test:", not np.any(np.logical_and(mask_coreC, mask_I_BC)))
	#print("Test:", not np.any(np.logical_and(mask_coreC, mask_I_ABC)))

	# early-phase weights
	h_coreA = h[mask_coreA]
	h_coreB = h[mask_coreB]
	h_coreC = h[mask_coreC]
	h_I_AB = h[mask_I_AB]
	h_I_AC = h[mask_I_AC]
	h_I_BC = h[mask_I_BC]
	h_I_ABC = h[mask_I_ABC]
	h_control = h[mask_control]

	# late-phase weights
	z_coreA = z[mask_coreA]
	z_coreB = z[mask_coreB]
	z_coreC = z[mask_coreC]
	z_I_AB = z[mask_I_AB]
	z_I_AC = z[mask_I_AC]
	z_I_BC = z[mask_I_BC]
	z_I_ABC = z[mask_I_ABC]
	z_control = z[mask_control]

	# total synaptic weights
	w_coreA = h_coreA + z_coreA
	w_coreB = h_coreB + z_coreB
	w_coreC = h_coreC + z_coreC
	w_I_AB = h_I_AB + z_I_AB
	w_I_AC = h_I_AC + z_I_AC
	w_I_BC = h_I_BC + z_I_BC
	w_I_ABC = h_I_ABC + z_I_ABC
	w_control = h_control + z_control

	# mean and standard deviation of the subpopulations (to compare to values from adjacencyFunctionsAttractors.py)
	#mean_z_coreA = np.mean(z_coreA)
	#mean_z_coreB = np.mean(z_coreB)
	#mean_z_coreC = np.mean(z_coreC)
	#mean_z_I_AB = np.mean(z_I_AB)
	#mean_z_I_AC = np.mean(z_I_AC)
	#mean_z_I_BC = np.mean(z_I_BC)
	#mean_z_I_ABC = np.mean(z_I_ABC)
	#mean_z_control = np.mean(z_control)
	#sd_z_coreA = np.std(z_coreA)
	#sd_z_coreB = np.std(z_coreB)
	#sd_z_coreC = np.std(z_coreC)
	#sd_z_I_AB = np.std(z_I_AB)
	#sd_z_I_AC = np.std(z_I_AC)
	#sd_z_I_BC = np.std(z_I_BC)
	#sd_z_I_ABC = np.std(z_I_ABC)
	#sd_z_control = np.std(z_control)

	# discretization of the distribution
	if bins is None:
		binh = np.linspace(np.min(h[connections]), np.max(h), 101, endpoint=True) # create range of bins for marginalProbDist(h...)
		binz = np.linspace(np.min(z[connections]), np.max(z), 101, endpoint=True) # create range of bins for marginalProbDist(z...)
		binw = np.linspace(np.min(w[connections]), np.max(w), 101, endpoint=True) # create range of bins for marginalProbDist(w...)
	else:
		[binh, binz, binw] = bins # use pre-defined bins

	hstep = binh[1]-binh[0]
	zstep = binz[1]-binz[0]
	wstep = binw[1]-binw[0]
	valh = np.delete(binh, -1) + hstep/2 # use mean values instead of lower bounds of the bins as values
	valz = np.delete(binz, -1) + zstep/2 # use mean values instead of lower bounds of the bins as values
	valw = np.delete(binw, -1) + wstep/2 # use mean values instead of lower bounds of the bins as values

	numconn = len(h[connections]) # normalization constant

	buf, ph_coreA = marginalProbDist(h_coreA, binning = True, bin_edges = binh, norm = numconn)
	buf, pz_coreA = marginalProbDist(z_coreA, binning = True, bin_edges = binz, norm = numconn)
	buf, pw_coreA = marginalProbDist(w_coreA, binning = True, bin_edges = binw, norm = numconn)
	if coreB is not None:
		buf, ph_coreB = marginalProbDist(h_coreB, binning = True, bin_edges = binh, norm = numconn)
		buf, pz_coreB = marginalProbDist(z_coreB, binning = True, bin_edges = binz, norm = numconn)
		buf, pw_coreB = marginalProbDist(w_coreB, binning = True, bin_edges = binw, norm = numconn)
		if h_I_AB.size > 0:
			buf, ph_I_AB = marginalProbDist(h_I_AB, binning = True, bin_edges = binh, norm = numconn)
			buf, pz_I_AB = marginalProbDist(z_I_AB, binning = True, bin_edges = binz, norm = numconn)
			buf, pw_I_AB = marginalProbDist(w_I_AB, binning = True, bin_edges = binw, norm = numconn)
		if coreC is not None:
			buf, ph_coreC = marginalProbDist(h_coreC, binning = True, bin_edges = binh, norm = numconn)
			buf, pz_coreC = marginalProbDist(z_coreC, binning = True, bin_edges = binz, norm = numconn)
			buf, pw_coreC = marginalProbDist(w_coreC, binning = True, bin_edges = binw, norm = numconn)
			if h_I_AC.size > 0:
				buf, ph_I_AC = marginalProbDist(h_I_AC, binning = True, bin_edges = binh, norm = numconn)
				buf, pz_I_AC = marginalProbDist(z_I_AC, binning = True, bin_edges = binz, norm = numconn)
				buf, pw_I_AC = marginalProbDist(w_I_AC, binning = True, bin_edges = binw, norm = numconn)
			if h_I_BC.size > 0:
				buf, ph_I_BC = marginalProbDist(h_I_BC, binning = True, bin_edges = binh, norm = numconn)
				buf, pz_I_BC = marginalProbDist(z_I_BC, binning = True, bin_edges = binz, norm = numconn)
				buf, pw_I_BC = marginalProbDist(w_I_BC, binning = True, bin_edges = binw, norm = numconn)
			if h_I_ABC.size > 0:
				buf, ph_I_ABC = marginalProbDist(h_I_ABC, binning = True, bin_edges = binh, norm = numconn)
				buf, pz_I_ABC = marginalProbDist(z_I_ABC, binning = True, bin_edges = binz, norm = numconn)
				buf, pw_I_ABC = marginalProbDist(w_I_ABC, binning = True, bin_edges = binw, norm = numconn)
	buf, ph_control = marginalProbDist(h_control, binning = True, bin_edges = binh, norm = numconn)
	buf, pz_control = marginalProbDist(z_control, binning = True, bin_edges = binz, norm = numconn)
	buf, pw_control = marginalProbDist(w_control, binning = True, bin_edges = binw, norm = numconn)

	# Write weight distributions to files
	fh = open(timestamp + "_eweight_dist_" + time_for_readout + add + ".txt", "w")
	fz = open(timestamp + "_lweight_dist_" + time_for_readout + add + ".txt", "w")
	fw = open(timestamp + "_totweight_dist_" + time_for_readout + add + ".txt", "w")

	for i in range(len(valh)):
		fh.write(str(valh[i]) + "\t\t" + str(ph_coreA[i]) + "\t\t")
		fz.write(str(valz[i]) + "\t\t" + str(pz_coreA[i]) + "\t\t")
		fw.write(str(valw[i]) + "\t\t" + str(pw_coreA[i]) + "\t\t")

		if coreB is not None:
			fh.write(str(ph_coreB[i]) + "\t\t")
			fz.write(str(pz_coreB[i]) + "\t\t")
			fw.write(str(pw_coreB[i]) + "\t\t")

			if h_I_AB.size > 0:
				fh.write(str(ph_I_AB[i]) + "\t\t")
				fz.write(str(pz_I_AB[i]) + "\t\t")
				fw.write(str(pw_I_AB[i]) + "\t\t")
			else:
				fh.write("nan\t\t")
				fz.write("nan\t\t")
				fw.write("nan\t\t")

			if coreC is not None:
				fh.write(str(ph_coreC[i]) + "\t\t")
				fz.write(str(pz_coreC[i]) + "\t\t")
				fw.write(str(pw_coreC[i]) + "\t\t")

				if h_I_AC.size > 0:
					fh.write(str(ph_I_AC[i]) + "\t\t")
					fz.write(str(pz_I_AC[i]) + "\t\t")
					fw.write(str(pw_I_AC[i]) + "\t\t")
				else:
					fh.write("nan\t\t")
					fz.write("nan\t\t")
					fw.write("nan\t\t")

				if h_I_BC.size > 0:
					fh.write(str(ph_I_BC[i]) + "\t\t")
					fz.write(str(pz_I_BC[i]) + "\t\t")
					fw.write(str(pw_I_BC[i]) + "\t\t")
				else:
					fh.write("nan\t\t")
					fz.write("nan\t\t")
					fw.write("nan\t\t")

				if h_I_ABC.size > 0:
					fh.write(str(ph_I_ABC[i]) + "\t\t")
					fz.write(str(pz_I_ABC[i]) + "\t\t")
					fw.write(str(pw_I_ABC[i]) + "\t\t")
				else:
					fh.write("nan\t\t")
					fz.write("nan\t\t")
					fw.write("nan\t\t")
			else:
				fh.write("nan\t\tnan\t\tnan\t\tnan\t\t")
				fz.write("nan\t\tnan\t\tnan\t\tnan\t\t")
				fw.write("nan\t\tnan\t\tnan\t\tnan\t\t")
		else:
			fh.write("nan\t\tnan\t\tnan\t\tnan\t\tnan\t\tnan\t\t")
			fz.write("nan\t\tnan\t\tnan\t\tnan\t\tnan\t\tnan\t\t")
			fw.write("nan\t\tnan\t\tnan\t\tnan\t\tnan\t\tnan\t\t")

		fh.write(str(ph_control[i]) + "\n")
		fz.write(str(pz_control[i]) + "\n")
		fw.write(str(pw_control[i]) + "\n")

	fh.close()
	fz.close()
	fw.close()

	# write gnuplot script
	f = open(timestamp + "_plot_dist.gpl", "w")
	f.write("### DO NOT EDIT THIS FILE! IT WILL BE OVERWRITTEN. ###\n\n" + \
	        "#set terminal png size 1024,640 enhanced\nset terminal pdf enhanced\n\n" + \
	        "#set style fill transparent solid 0.8 noborder # for 'boxes' style\n" + \
			"#set style fill transparent pattern 4 bo # for 'boxes' style\n" + \
	        "set log y\nset format y \"%.0e\"\nset yrange [3e-06:1]\nset key outside\n\n" + \
			"h_0 = " + str(h_0) + "\n" +\
			"epsilon = " + str(epsilon) + "\n\n")

	# plotting of early-phase weight distribution
	f.write("set output \"" + timestamp + "_eweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binh[0]-10*hstep) + "/h_0:" + str(binh[-1]+10*hstep) + "/h_0]\n")
	f.write("set xlabel \"Early-phase weight / h_0\"\nset ylabel \"Relative frequency\"\n")
	f.write("plot \"" + timestamp + "_eweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0):($2 > 0 ? $2 : epsilon) t \"Within ~A{.8\\\\~}\" lc 6 with histeps, \\\n")
	if coreB is not None:
		f.write("\"\" using ($1/h_0):($3 > 0 ? $3 : epsilon) t \"Within ~B{.8\\\\~}\" lc 7 with histeps, \\\n")
		if coreC is not None:
			f.write("\"\" using ($1/h_0):($5 > 0 ? $5 : epsilon) t \"Within ~C{.8\\\\~}\" lc 5 with histeps, \\\n")
		if h_I_AB.size > 0:
			f.write("\"\" using ($1/h_0):($4 > 0 ? $4 : epsilon) t \"Within I_{AB}\" lc 1 with histeps, \\\n")
		if coreC is not None:
			if h_I_AC.size > 0:
				f.write("\"\" using ($1/h_0):($6 > 0 ? $6 : epsilon) t \"Within I_{AC}\" lc 2 with histeps, \\\n")
			if h_I_BC.size > 0:
				f.write("\"\" using ($1/h_0):($7 > 0 ? $7 : epsilon) t \"Within I_{BC}\" lc 4 with histeps, \\\n")
			if h_I_ABC.size > 0:
				f.write("\"\" using ($1/h_0):($8 > 0 ? $8 : epsilon) t \"Within I_{ABC}\" lc 3 with histeps, \\\n")
	f.write("\"\" using ($1/h_0):($9 > 0 ? $9 : epsilon) t \"Others\" lc rgb \"#eeeeee\" with histeps\n")

	# plotting of late-phase weight distribution
	f.write("\nset output \"" + timestamp + "_lweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binz[0]-10*zstep) + "/h_0:" + str(binz[-1]+10*zstep) + "/h_0]\n")
	f.write("set xlabel \"Late-phase weight / h_0\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_lweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0):($2 > 0 ? $2 : epsilon) t \"Within ~A{.8\\\\~}\" lc 6 with histeps, \\\n")
	if coreB is not None:
		f.write("\"\" using ($1/h_0):($3 > 0 ? $3 : epsilon) t \"Within ~B{.8\\\\~}\" lc 7 with histeps, \\\n")
		if coreC is not None:
			f.write("\"\" using ($1/h_0):($5 > 0 ? $5 : epsilon) t \"Within ~C{.8\\\\~}\" lc 5 with histeps, \\\n")
		if h_I_AB.size > 0:
			f.write("\"\" using ($1/h_0):($4 > 0 ? $4 : epsilon) t \"Within I_{AB}\" lc 1 with histeps, \\\n")
		if coreC is not None:
			if h_I_AC.size > 0:
				f.write("\"\" using ($1/h_0):($6 > 0 ? $6 : epsilon) t \"Within I_{AC}\" lc 2 with histeps, \\\n")
			if h_I_BC.size > 0:
				f.write("\"\" using ($1/h_0):($7 > 0 ? $7 : epsilon) t \"Within I_{BC}\" lc 4 with histeps, \\\n")
			if h_I_ABC.size > 0:
				f.write("\"\" using ($1/h_0):($8 > 0 ? $8 : epsilon) t \"Within I_{ABC}\" lc 3 with histeps, \\\n")
	f.write("\"\" using ($1/h_0):($9 > 0 ? $9 : epsilon) t \"Others\" lc rgb \"#eeeeee\" with histeps\n")

	# plotting of total weight distribution
	f.write("\nset output \"" + timestamp + "_totweight_dist_" + time_for_readout + add + ".pdf\"\n")
	f.write("set xrange [" + str(binw[0]-10*wstep) + "/h_0*100:" + str(binw[-1]+10*wstep) + "/h_0*100]\n")
	f.write("set xlabel \"Total synaptic weight (%)\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_totweight_dist_" + time_for_readout + add + ".txt\" using ($1/h_0*100):($2 > 0 ? $2 : epsilon) t \"Within ~A{.8\\\\~}\" lc 6 with histeps, \\\n")
	if coreB is not None:
		f.write("\"\" using ($1/h_0*100):($3 > 0 ? $3 : epsilon) t \"Within ~B{.8\\\\~}\" lc 7 with histeps, \\\n")
		if coreC is not None:
			f.write("\"\" using ($1/h_0*100):($5 > 0 ? $5 : epsilon) t \"Within ~C{.8\\\\~}\" lc 5 with histeps, \\\n")
		if h_I_AB.size > 0:
			f.write("\"\" using ($1/h_0*100):($4 > 0 ? $4 : epsilon) t \"Within I_{AB}\" lc 1 with histeps, \\\n")
		if coreC is not None:
			if h_I_AC.size > 0:
				f.write("\"\" using ($1/h_0*100):($6 > 0 ? $6 : epsilon) t \"Within I_{AC}\" lc 2 with histeps, \\\n")
			if h_I_BC.size > 0:
				f.write("\"\" using ($1/h_0*100):($7 > 0 ? $7 : epsilon) t \"Within I_{BC}\" lc 4 with histeps, \\\n")
			if h_I_ABC.size > 0:
				f.write("\"\" using ($1/h_0*100):($8 > 0 ? $8 : epsilon) t \"Within I_{ABC}\" lc 3 with histeps, \\\n")
	f.write("\"\" using ($1/h_0*100):($9 > 0 ? $9 : epsilon) t \"Others\" lc rgb \"#eeeeee\" with histeps\n")

	f.close()

	call(["gnuplot", timestamp + "_plot_dist.gpl"])


	os.chdir(orgdir) # change back to original directory

# mapToBins
# Adjusts the values of an array such that they match bins defined by another array
# a: array of values
# b: array of bin edges (including the terminal one)
# return: array of the shape of a with values discretized according to b
def mapToBins(a, b):
	a = a.flatten()

	for i in range(len(a)):

		if a[i] < b[0] or a[i] > b[len(b)-1]:
			raise ValueError("Value " + a[i] + " at index " + i + " is out of range.")

		for j in reversed(range(len(b)-1)):
			if a[i] > b[j]:
				a[i] = (b[j+1] + b[j]) / 2
				break
	return a

# marginalProbDist
# Computes a marginal probability distribution
# a: array of outcomes (e.g., array of activities of all neurons in a network or array of weights of all synapses in a network)
# binning [optional]: specifies if bins are used to discretize the values of the distribution
# bin_edges [optional]: pre-specified range of binning edges; only applicable if binning is used
# norm [optional]: value by which to normalize
# return: array of values (unless bin_edges were pre-defined) and corresponding array of their probabilities
def marginalProbDist(a, binning = False, bin_edges = None, norm = None):

	if binning == True:
		if bin_edges is None: # generate bin_edges for 100 bins
			bin_edges = np.linspace(np.min(h), np.max(h), 101, endpoint=True) # create array of bins (defined by their edges)
			values_mean = np.delete(binh, -1) + (np.max(a)-np.min(a)) / 100 / 2 # use mean values instead of lower bounds as values of the bins
		else:
			values_mean = None

		freq_of_values = np.histogram(a, bin_edges)[0] # create histogram of activity value occurrences
		#values_mean, freq_of_values = np.unique(mapToBins(a, bin_edges), return_counts=True, axis=0) # yields the same result but without empty bins
	else:
		values_mean, freq_of_values = np.unique(a, return_counts=True) # determine and count occuring activities

	if norm is None:
		freq_of_values = freq_of_values / np.sum(freq_of_values) # normalize only over this quantity
	else:
		freq_of_values = freq_of_values / norm # normalize over given number (i.e., typically, over all values)

	return values_mean, freq_of_values

# jointProbDist
# Computes a joint probability distribution
# a1: first array of outcomes (e.g., array of activities of all neurons in a network or array of weights of all synapses in a network)
# a2: second array of outcomes
# binning [optional]: specifies if bins are used to discretize the values of the distribution
# return: array of value pairs and corresponding array of their probabilities
def jointProbDist(a1, a2, binning = False):

	if binning == True:
		ab = np.concatenate((a1,a2))

		try:
			values = np.arange(np.min(ab), np.max(ab), (np.max(ab)-np.min(ab)) / 100) # create array of activity value bins
		except ValueError:
			values = np.array([np.min(ab)])
		bin_edges = np.concatenate((values, np.max(ab)), axis=None) # add last edge
		a1 = mapToBins(a1, bin_edges)
		a2 = mapToBins(a2, bin_edges)

		ja = np.array([a1.flatten(), a2.flatten()]).T # array of pairs of the two outcomes (for one neuron or synapse or ... at two times)

		value_pairs, freq_of_values = np.unique(ja, return_counts=True, axis=0) # determine and count occurring activities
		freq_of_values = freq_of_values / np.sum(freq_of_values) # normalize
	else:
		ja = np.array([a1.flatten(),a2.flatten()]).T # array of pairs of the two outcomes (for one neuron or synapse or ... at two times)

		value_pairs, freq_of_values = np.unique(ja, return_counts=True, axis=0) # determine and count occurring activities
		freq_of_values = freq_of_values / np.sum(freq_of_values) # normalize
		#dist = np.asarray(value_pairs, freq_of_values).T # create distribution

	return value_pairs, freq_of_values

# main
# argv[]: comandline arguments
'''### examples:
if __name__ == "__main__":

	######################################################################################################
	## early- and late-phase weight distributions, four each, with overall normalization,
	## as used in Luboeinski and Tetzlaff, Commun. Biol., 2021

	Nl_exc = 40 # number of excitatory neurons in one line of a square
	h_0 = 4.20075 # initial/median synaptic weight
	core = np.arange(150) # size of the assembly

	if len(sys.argv) == 3:
		
		ts1 = str(sys.argv[1]) # timestamp for simulation data before consolidation
		ts2 = str(sys.argv[2]) # timestamp for simulation data after consolidation

		# define bins (can also be done automatically, but then bins are different for each sample)
		binh = np.linspace(0.00212467, 0.74913949, 101, endpoint=True) # bins for early-phase weights
		binz = np.linspace(0.00000000, 0.33602527, 101, endpoint=True) # bins for late-phase weights
		binw = np.linspace(0.00212467, 0.76000000, 101, endpoint=True) # bins for total synaptic weights
		binv = np.linspace(0.0, 102.0, 51, endpoint=True) # bins for activities
		bins = [binh, binz, binw, binv]

		# compute and plot distributions
		plotDistributions("10s", ts1, "_150default", Nl_exc, "20.0", core, h_0, bins=bins, norm_all=True) # before 10s-recall
		plotDistributions("10s", ts1, "_150default", Nl_exc, "20.1", core, h_0, bins=bins, norm_all=True) # after 10s-recall
		plotDistributions("8h", ts2, "_150default", Nl_exc, "28810.0", core, h_0, bins=bins, norm_all=True) # before 8h-recall
		plotDistributions("8h", ts2, "_150default", Nl_exc, "28810.1", core, h_0, bins=bins, norm_all=True) # after 8h-recall

	######################################################################################################
	## weight and activity distributions with automatic determination of bins and without overall 
	## normalization
	
	Nl_exc = 40 # number of excitatory neurons in one line of a square
	h_0 = 4.20075 # initial/median synaptic weight
	core = np.arange(150) # size of the assembly

	### 28810.0 ###
	minmax = findOverallMinMax(".", Nl_exc, "28810.0", h_0)
	binh = np.linspace(minmax[0][0], minmax[0][1], 101, endpoint=True) # bins for early-phase weights
	binz = np.linspace(minmax[1][0], minmax[1][1], 101, endpoint=True) # bins for late-phase weights
	binw = np.linspace(minmax[2][0], minmax[2][1], 101, endpoint=True) # bins for total synaptic weights
	binv = np.linspace(minmax[3][0], minmax[3][1], 51, endpoint=True) # bins for activities
	bins = [binh, binz, binw, binv]
	np.save("bins_28810.0.npy", bins)
	np.savetxt("minmax_28810.0.txt", minmax)
	#bins = np.load("bins_28810.0.npy", allow_pickle=True)

	plotDistributions("nm=0.02", "any", "", Nl_exc, h_0, "28810.0", core, bins=bins, norm_all=False) # before 8h-recall
	plotDistributions("nm=0.06", "any", "", Nl_exc, h_0, "28810.0", core, bins=bins, norm_all=False) # before 8h-recall
	plotDistributions("nm=0.10", "any", "", Nl_exc, h_0, "28810.0", core, bins=bins, norm_all=False) # before 8h-recall
	plotDistributions("nm=0.20", "any", "", Nl_exc, h_0, "28810.0", core, bins=bins, norm_all=False) # before 8h-recall

	### 28810.1 ###
	minmax = findOverallMinMax(".", Nl_exc, "28810.1", h_0)
	binh = np.linspace(minmax[0][0], minmax[0][1], 101, endpoint=True) # bins for early-phase weights
	binz = np.linspace(minmax[1][0], minmax[1][1], 101, endpoint=True) # bins for late-phase weights
	binw = np.linspace(minmax[2][0], minmax[2][1], 101, endpoint=True) # bins for total synaptic weights
	binv = np.linspace(minmax[3][0], minmax[3][1], 51, endpoint=True) # bins for activities
	bins = [binh, binz, binw, binv]
	np.save("bins_28810.1.npy", bins)
	np.savetxt("minmax_28810.1.txt", minmax)
	#bins = np.load("bins_28810.1.npy", allow_pickle=True)

	plotDistributions("nm=0.02", "any", "", Nl_exc, h_0, "28810.1", core, bins=bins, norm_all=False) # after 8h-recall
	plotDistributions("nm=0.06", "any", "", Nl_exc, h_0, "28810.1", core, bins=bins, norm_all=False) # after 8h-recall
	plotDistributions("nm=0.10", "any", "", Nl_exc, h_0, "28810.1", core, bins=bins, norm_all=False) # after 8h-recall
	plotDistributions("nm=0.20", "any", "", Nl_exc, h_0, "28810.1", core, bins=bins, norm_all=False) # after 8h-recall

	### 20.0 ###
	minmax = findOverallMinMax(".", Nl_exc, "20.0", h_0)
	binh = np.linspace(minmax[0][0], minmax[0][1], 101, endpoint=True) # bins for early-phase weights
	binz = np.linspace(minmax[1][0], minmax[1][1], 101, endpoint=True) # bins for late-phase weights
	binw = np.linspace(minmax[2][0], minmax[2][1], 101, endpoint=True) # bins for total synaptic weights
	binv = np.linspace(minmax[3][0], minmax[3][1], 51, endpoint=True) # bins for activities
	bins = [binh, binz, binw, binv]
	np.save("bins_20.0.npy", bins)
	np.savetxt("minmax_20.0.txt", minmax)
	#bins = np.load("bins_20.0.npy", allow_pickle=True)

	plotDistributions("nm=0.02", "any", "", Nl_exc, h_0, "20.0", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.06", "any", "", Nl_exc, h_0, "20.0", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.10", "any", "", Nl_exc, h_0, "20.0", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.20", "any", "", Nl_exc, h_0, "20.0", core, bins=bins, norm_all=False) # before 10s-recall

	### 20.1 ###
	minmax = findOverallMinMax(".", Nl_exc, "20.1", h_0)
	binh = np.linspace(minmax[0][0], minmax[0][1], 101, endpoint=True) # bins for early-phase weights
	binz = np.linspace(minmax[1][0], minmax[1][1], 101, endpoint=True) # bins for late-phase weights
	binw = np.linspace(minmax[2][0], minmax[2][1], 101, endpoint=True) # bins for total synaptic weights
	binv = np.linspace(minmax[3][0], minmax[3][1], 51, endpoint=True) # bins for activities
	bins = [binh, binz, binw, binv]
	np.save("bins_20.1.npy", bins)
	np.savetxt("minmax_20.1.txt", minmax)
	#bins = np.load("bins_20.1.npy", allow_pickle=True)

	plotDistributions("nm=0.02", "any", "", Nl_exc, h_0, "20.1", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.06", "any", "", Nl_exc, h_0, "20.1", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.10", "any", "", Nl_exc, h_0, "20.1", core, bins=bins, norm_all=False) # before 10s-recall
	plotDistributions("nm=0.20", "any", "", Nl_exc, h_0, "20.1", core, bins=bins, norm_all=False) # before 10s-recall

	######################################################################################################
	## Network with 3 cell assemblies - weight distributions with automatic determination of bins 
	## and without overall normalization

	Nl_exc = 50 # number of excitatory neurons in one line of a square
	h_0 = 4.20075 # initial/median synaptic weight
	core_size = 600 # size of the assembly
	core1 = np.arange(core_size)
	core2 = np.arange(core_size, 2*core_size)
	core3 = np.arange(2*core_size, 3*core_size)

	minmax = findOverallMinMax(".", Nl_exc, "28810.0", h_0)
	binh = np.linspace(minmax[0][0], minmax[0][1], 101, endpoint=True) # bins for early-phase weights
	binz = np.linspace(minmax[1][0], minmax[1][1], 101, endpoint=True) # bins for late-phase weights
	binw = np.linspace(minmax[2][0], minmax[2][1], 101, endpoint=True) # bins for total synaptic weights
	bins = [binh, binz, binw]
	np.save("bins_28810.0.npy", bins)
	np.savetxt("minmax_28810.0.txt", minmax)
	#bins = np.load("bins_28810.0.npy", allow_pickle=True)

	plotWeightDistributions3CAs("stdconn", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn1", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn2", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn3", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn4", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn5", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn6", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn7", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn8", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
	plotWeightDistributions3CAs("altconn9", "any", "", Nl_exc, h_0, "28810.0", core1, core2, core3, bins=bins)
'''
