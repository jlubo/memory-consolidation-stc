############################################################################################
### Functions to analyze and plot weight and activity distributions from simulation data ###
############################################################################################

### Copyright 2019-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

from readWeightData import *
import sys
import os.path
import numpy as np
from pathlib import Path
from subprocess import call

# plotDistributions
# Creates data and plot files of the weight and activity distribution at a given time
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# add: additional descriptor
# Nl: the number of excitatory neurons in one line of a quadratic grid
# time: the time that at which the weights shall be read out
# core: array of indices of the cell assembly (core) neurons
def plotDistributions(nppath, timestamp, add, Nl, time, core):

	# look for data file [timestamp]_net_[time].txt
	path = ""
	rawpaths = Path(nppath)

	for x in rawpaths.iterdir():

		tmppath = str(x)

		if (timestamp + "_net_" + time + ".txt") in tmppath:
			path = tmppath

	if path == "":
		raise FileNotFoundError('"' + timestamp + '_net_' + time + '.txt" was not found')

	# read data from file
	try:
		connections, h, z, v = readWeightMatrixData(path, Nl)

	except ValueError:
		raise
	except OSError:
		raise

	# determine subpopulations
	N_tot = Nl**2 # total number of neurons
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

	h_CA_within = h[mask_CA_within]
	h_CA_outgoing = h[mask_CA_outgoing]
	h_CA_incoming = h[mask_CA_incoming]
	h_control = h[mask_control]

	z_CA_within = z[mask_CA_within]
	z_CA_outgoing = z[mask_CA_outgoing]
	z_CA_incoming = z[mask_CA_incoming]
	z_control = z[mask_control]

	v_CA = v.flatten()[np.in1d(all, core)]
	v_control = v.flatten()[np.logical_not(np.in1d(all, core))]

	binh = np.concatenate((np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h)) / 100), np.max(h)), axis=None) # create range of bins for marginalProbDist(h...)
	binz = np.concatenate((np.arange(np.min(z), np.max(z), (np.max(z)-np.min(z)) / 100), np.max(z)), axis=None) # create range of bins for marginalProbDist(z...)
	binv = np.concatenate((np.arange(np.min(v), np.max(v), (np.max(v)-np.min(v)) / 100), np.max(v)), axis=None) # create range of bins for marginalProbDist(v...)

	valh = np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h)) / 100) + (np.max(h)-np.min(h)) / 200 # use mean values instead of lower bounds of the bins as values
	valz = np.arange(np.min(z), np.max(z), (np.max(z)-np.min(z)) / 100) + (np.max(z)-np.min(z)) / 200 # use mean values instead of lower bounds of the bins as values
	valv = np.arange(np.min(v), np.max(v), (np.max(v)-np.min(v)) / 100) + (np.max(v)-np.min(v)) / 200 # use mean values instead of lower bounds of the bins as values

	buf, ph_CA_within = marginalProbDist(h_CA_within, binning = True, bin_edges = binh)
	buf, ph_CA_outgoing = marginalProbDist(h_CA_outgoing, binning = True, bin_edges = binh)
	buf, ph_CA_incoming = marginalProbDist(h_CA_incoming, binning = True, bin_edges = binh)
	buf, ph_control = marginalProbDist(h_control, binning = True, bin_edges = binh)
	buf, pz_CA_within = marginalProbDist(z_CA_within, binning = True, bin_edges = binz)
	buf, pz_CA_outgoing = marginalProbDist(z_CA_outgoing, binning = True, bin_edges = binz)
	buf, pz_CA_incoming = marginalProbDist(z_CA_incoming, binning = True, bin_edges = binz)
	buf, pz_control = marginalProbDist(z_control, binning = True, bin_edges = binz)
	buf, pv_CA = marginalProbDist(v_CA, binning = True, bin_edges = binv)
	buf, pv_control = marginalProbDist(v_control, binning = True, bin_edges = binv)

	f = open(timestamp + "_eweight_dist_" + time + add + ".txt", "w")
	for i in range(len(valh)):
		f.write(str(valh[i]) + "\t\t" + str(ph_CA_within[i]) + "\t\t" + str(ph_CA_outgoing[i]) + "\t\t" + \
		        str(ph_CA_incoming[i]) + "\t\t" + str(ph_control[i]) + "\n")
	f.close()

	f = open(timestamp + "_lweight_dist_" + time + add + ".txt", "w")
	for i in range(len(valz)):
		f.write(str(valz[i]) + "\t\t" + str(pz_CA_within[i]) + "\t\t" + str(pz_CA_outgoing[i]) + "\t\t" + \
		        str(pz_CA_incoming[i]) + "\t\t" + str(pz_control[i]) + "\n")
	f.close()

	f = open(timestamp + "_act_dist_" + time + add + ".txt", "w")
	for i in range(len(valv)):
		f.write(str(valv[i]) + "\t\t" + str(pv_CA[i]) + "\t\t" + str(pv_control[i]) + "\n")
	f.close()

	if os.path.exists("plot_dist.gpl"):
		f = open("plot_dist.gpl", "a")
	else:
		f = open("plot_dist.gpl", "w")
		f.write("#set terminal png size 1024,640 enhanced\nset terminal pdf enhanced\n\n" + \
		        "set style fill transparent solid 0.8 noborder\n" + \
		        "set log y\nset format y \"%.0e\"\nset yrange [3e-06:1]\nset key outside\n\n")

	f.write("set output \"" + timestamp + "_eweight_dist_" + time + add + ".pdf\"\n")
	f.write("set xlabel \"Early-phase weight / nC\"\nset ylabel \"Relative frequency\"\n")
	f.write("plot [0.3:0.9] \"" + timestamp + "_eweight_dist_" + time + add + ".txt\" using 1:($1 > 0 ? $2 : $2) t \"CA\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $3 : $3) t \"outgoing\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $4 : $4) t \"incoming\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $5 : $5) t \"control\" with boxes\n")

	f.write("\nset output \"" + timestamp + "_lweight_dist_" + time + add + ".pdf\"\n")
	f.write("set xlabel \"Late-phase weight\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_lweight_dist_" + time + add + ".txt\" using 1:($1 > 0 ? $2 : $2) t \"CA\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $3 : $3) t \"outgoing\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $4 : $4) t \"incoming\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $5 : $5) t \"control\" with boxes\n")

	f.write("\nset output \"" + timestamp + "_act_dist_" + time + add + ".pdf\"\n")
	f.write("set xlabel \"Neuronal firing rate / Hz\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_act_dist_" + time + add + ".txt\" using 1:2 t \"CA\" with boxes, " + \
	        "\"\" using 1:3 t \"control\" with boxes\n\n")

	f.close()

	call(["gnuplot", "plot_dist.gpl"])

# plotDistributions3CAs
# Creates data and plot files of the weight distribution of a network with 3 overlapping assemblies at a given time
# nppath: path to the network_plots directory to read the data from
# timestamp: a string containing date and time (to access correct paths)
# add: additional descriptor
# Nl: the number of excitatory neurons in one line of a quadratic grid
# time: the time that at which the weights shall be read out
# coreA: array of indices of the first cell assembly (core) neurons
# coreB: array of indices of the second cell assembly (core) neurons
# coreC: array of indices of the third cell assembly (core) neurons
def plotDistributions3CAs(nppath, timestamp, add, Nl, time, coreA, coreB, coreC):

	# look for data file [timestamp]_net_[time].txt
	path = ""
	rawpaths = Path(nppath)

	for x in rawpaths.iterdir():

		tmppath = str(x)

		if (timestamp + "_net_" + time + ".txt") in tmppath:
			path = tmppath

	if path == "":
		raise FileNotFoundError('"' + timestamp + '_net_' + time + '.txt" was not found')

	# read data from file
	try:
		connections, h, z, v = readWeightMatrixData(path, Nl)

	except ValueError:
		raise
	except OSError:
		raise

	# determine synapses within the cell assemblies
	N_tot = Nl**2 # total number of neurons

	mask_coreA = np.zeros((N_tot, N_tot), dtype=bool)
	for syn_pre in coreA:
		for syn_post in coreA:
			mask_coreA[syn_pre,syn_post] = True

	mask_coreB = np.zeros((N_tot, N_tot), dtype=bool)
	for syn_pre in coreB:
		for syn_post in coreB:
			mask_coreB[syn_pre,syn_post] = True

	mask_coreC = np.zeros((N_tot, N_tot), dtype=bool)
	for syn_pre in coreC:
		for syn_post in coreC:
			mask_coreC[syn_pre,syn_post] = True

	h_coreA = h[mask_coreA]
	h_coreB = h[mask_coreB]
	h_coreC = h[mask_coreC]

	z_coreA = z[mask_coreA]
	z_coreB = z[mask_coreB]
	z_coreC = z[mask_coreC]

	binh = np.concatenate((np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h)) / 100), np.max(h)), axis=None) # create range of bins for marginalProbDist(h...)
	binz = np.concatenate((np.arange(np.min(z), np.max(z), (np.max(z)-np.min(z)) / 100), np.max(z)), axis=None) # create range of bins for marginalProbDist(z...)

	valh = np.arange(np.min(h), np.max(h), (np.max(h)-np.min(h)) / 100) + (np.max(h)-np.min(h)) / 200 # use mean values instead of lower bounds of the bins as values
	valz = np.arange(np.min(z), np.max(z), (np.max(z)-np.min(z)) / 100) + (np.max(z)-np.min(z)) / 200 # use mean values instead of lower bounds of the bins as values

	buf, ph_coreA = marginalProbDist(h_coreA, binning = True, bin_edges = binh)
	buf, ph_coreB = marginalProbDist(h_coreB, binning = True, bin_edges = binh)
	buf, ph_coreC = marginalProbDist(h_coreC, binning = True, bin_edges = binh)
	buf, pz_coreA = marginalProbDist(z_coreA, binning = True, bin_edges = binz)
	buf, pz_coreB = marginalProbDist(z_coreB, binning = True, bin_edges = binz)
	buf, pz_coreC = marginalProbDist(z_coreC, binning = True, bin_edges = binz)

	f = open(timestamp + "_eweight_dist_" + time + add + ".txt", "w")
	for i in range(len(valh)):
		f.write(str(valh[i]) + "\t\t" + str(ph_coreA[i]) + "\t\t" + str(ph_coreB[i]) + "\t\t" + \
		        str(ph_coreC[i]) + "\n")
	f.close()

	f = open(timestamp + "_lweight_dist_" + time + add + ".txt", "w")
	for i in range(len(valz)):
		f.write(str(valz[i]) + "\t\t" + str(pz_coreA[i]) + "\t\t" + str(pz_coreB[i]) + "\t\t" + \
		        str(pz_coreC[i]) + "\n")
	f.close()

	if os.path.exists("plot_dist.gpl"):
		f = open("plot_dist.gpl", "a")
	else:
		f = open("plot_dist.gpl", "w")
		f.write("#set terminal png size 1024,640 enhanced\nset terminal pdf enhanced\n\n" + \
		        "set style fill transparent solid 0.8 noborder\n" + \
		        "set log y\nset format y \"%.0e\"\nset yrange [3e-06:1]\nset key outside\n\n")

	f.write("set output \"" + timestamp + "_eweight_dist_" + time + add + ".pdf\"\n")
	f.write("set xlabel \"Early-phase weight / nC\"\nset ylabel \"Relative frequency\"\n")
	f.write("plot [0.3:0.9] \"" + timestamp + "_eweight_dist_" + time + add + ".txt\" using 1:($1 > 0 ? $2 : $2) t \"A\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $3 : $3) t \"B\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $4 : $4) t \"C\" with boxes\n")

	f.write("\nset output \"" + timestamp + "_lweight_dist_" + time + add + ".pdf\"\n")
	f.write("set xlabel \"Late-phase weight\"\nset ylabel \"Relative frequency\"\nset format y \"%.0e\"\n")
	f.write("plot \"" + timestamp + "_lweight_dist_" + time + add + ".txt\" using 1:($1 > 0 ? $2 : $2) t \"A\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $3 : $3) t \"B\" with boxes, \\\n" + \
	        "\"\" using 1:($1 > 0 ? $4 : $4) t \"C\" with boxes\n")

	f.close()

	call(["gnuplot", "plot_dist.gpl"])

# mapToBins
# Maps the values of an array to bins defined by another array
# a: array of values
# b: array of bin edges (including the terminal one)
# return: array of the shape of a with values of b
def mapToBins(a, b):
	a = a.flatten()

	for i in range(len(a)):

		if a[i] < b[0] or a[i] > b[len(b)-1]:
			raise ValueError("Value " + a[i] + " at index " + i + " is out of range.")

		for j in reversed(range(len(b)-1)):
			if a[i] >= b[j]:
				a[i] = (b[j+1] + b[j]) / 2
				break
	return a

# marginalProbDist
# Computes a marginal probability distribution
# a: array of outcomes (e.g., array of activities of all neurons in a network or array of weights of all synapses in a network)
# binning [optional]: specifies if bins are used to discretize the values of the distribution
# bin_edges [optional]: pre-specified range of binning edges; only applicable if binning is used
# return: array of values and corresponding array of their probabilities
def marginalProbDist(a, binning = False, bin_edges = None):

	if binning == True:

		if np.max(a) != np.min(a):
			values_low = np.arange(np.min(a), np.max(a), (np.max(a)-np.min(a)) / 100) # create array of activity value bins
			values_mean = values_low + (np.max(a)-np.min(a)) / 200 # use mean values instead of lower bounds of the bins as values
		else:
			values_low = np.array([np.min(a)])
			values_mean = values_low

		if bin_edges is None:
			bin_edges = np.concatenate((values_low, np.max(a)), axis=None) # add last edge

		freq_of_values = np.histogram(a, bin_edges)[0] # create histogram of activity value occurrences

		#values, freq_of_values = np.unique(mapToBins(a, bin_edges), return_counts=True, axis=0) # yields the same result

		freq_of_values = freq_of_values / np.sum(freq_of_values) # normalize

	else:
		#a = np.sort(a)
		values_mean, freq_of_values = np.unique(a, return_counts=True) # determine and count occuring activities
		freq_of_values = freq_of_values / np.sum(freq_of_values) # normalize
		#dist = np.asarray(values, freq_of_values).T # create distribution

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
# Plots distributions from two simulations at different times
# argv[]: timestamps of two simulations

if __name__ == "__main__":

	if len(sys.argv) < 3:
		print("Not enough arguments provided!")
		exit()
	else:
		ts1 = str(sys.argv[1]) # timestamp for simulation data before consolidation
		ts2 = str(sys.argv[2]) # timestamp for simulation data after consolidation

	core = np.arange(150) # size of the assembly
	plotDistributions(".", ts1, "_150default", 40, "20.0", core) # before 10s-recall
	plotDistributions(".", ts1, "_150default", 40, "20.1", core) # after 10s-recall
	plotDistributions(".", ts2, "_150default", 40, "28810.0", core) # before 8h-recall
	plotDistributions(".", ts2, "_150default", 40, "28810.1", core) # after 8h-recall
