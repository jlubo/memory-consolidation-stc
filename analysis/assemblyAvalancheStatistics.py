###########################################################################################################################
### Script to extract the time series of avalanches in cell assemblies and the whole population from spike raster data, ###
###                       and to compute the likelihoods of avalanche occurrence and "transition",                      ###
###                             as well as the mean firing rates in the assemblies                                      ###
###########################################################################################################################

### Copyright 2020-2022 Jannik Luboeinski
### License: Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import os
import time
import sys
from pathlib import Path
from shutil import copy2
from overlapParadigms import *
from utilityFunctions import *

# main properties (most can be adjusted via commandline parameters, see at the end of the script)
paradigm = "NOOVERLAP" # paradigm of overlaps between assemblies
period_duration = 0.01 # binning period (in units of seconds)
n_thresh = 10 # number of spikes in a binning period to consider them an avalanche
new_plots = True # defines if new spike raster plots shall be created using gnuplot
exc_pop_size = 2500 # number of neurons in the excitatory population
core_size = 600 # total size of one cell assembly

# cell assemblies
coreA, coreB, coreC = coreDefinitions(paradigm, core_size)

# control population
mask_coreA = np.in1d(np.arange(exc_pop_size), coreA)
mask_coreB = np.in1d(np.arange(exc_pop_size), coreB)
mask_coreC = np.in1d(np.arange(exc_pop_size), coreC)
ctrl = np.arange(exc_pop_size)[np.logical_not(np.logical_or(mask_coreA, np.logical_or(mask_coreB, mask_coreC)))] # control neurons (neurons that are not within a cell assembly)
ctrl_size = len(ctrl)

####################################
# timeSeries
# Computes the time series of avalanche occurrence and saves it to a file
# timestamp: the timestamp of the simulation data
# spike_raster_file: the name of the spike raster file
# output_dir: relative path to the output directory
# return: the time series as a list of characters
def timeSeries(timestamp, spike_raster_file, output_dir):

	t0 = time.time()

	# read the last line and compute number of periods
	with open(spike_raster_file, 'rb') as f:
		f.seek(-2, os.SEEK_END)
		while f.read(1) != b'\n': # seek last line
			f.seek(-2, os.SEEK_CUR)
		last_line = f.readline().decode()
	num_periods_tot = int(float(last_line.split('\t\t')[0]) / period_duration) + 1

	# count lines
	with open(spike_raster_file) as f:
		num_rows = sum(1 for _ in f)
	print("num_rows =", num_rows)

	# counters per period for the different cell assemblies
	counterA = np.zeros(num_periods_tot, dtype=int)
	counterB = np.zeros(num_periods_tot, dtype=int)
	counterC = np.zeros(num_periods_tot, dtype=int)
	counterOverall = np.zeros(num_periods_tot, dtype=int)
	counterCtrl = np.zeros(num_periods_tot, dtype=int)
	series = ["-" for i in range(num_periods_tot)]

	# read all data
	f = open(spike_raster_file)
	for line in f:
		row = line.split('\t\t')
		t = float(row[0])
		n = int(row[1])

		current_period = int(np.floor(t / period_duration))

		if n < exc_pop_size:
			counterOverall[current_period] += 1

			if n in coreA:
				counterA[current_period] += 1

			if n in coreB:
				counterB[current_period] += 1

			if n in coreC:
				counterC[current_period] += 1

			if n in ctrl:
				counterCtrl[current_period] += 1

	f.close()

	# determine active CAs for each period and write data to file
	fout = open(os.path.join(output_dir, timestamp + '_CA_time_series.txt'), 'w')

	for i in range(num_periods_tot):
		fout.write(str(round((i+0.5)*period_duration,4)) + "\t\t") # write time at 1/2 of a period

		if counterOverall[i] > n_thresh:
			series[i] = "o"
			if counterC[i] > n_thresh:
				series[i] = "C" + series[i]
			if counterB[i] > n_thresh:
				series[i] = "B" + series[i]
			if counterA[i] > n_thresh:
				series[i] = "A" + series[i]
			
		fout.write(series[i] + "\t\t" + str(counterA[i]) + "\t\t" + str(counterB[i]) + "\t\t" + str(counterC[i]) + "\t\t" + str(counterOverall[i]) + "\n")

	fout.close()

	time_el = round(time.time()-t0) # elapsed time in seconds
	time_el_str = "Elapsed time: "
	if time_el < 60:
		time_el_str += str(time_el) + " s"
	else:
		time_el_str += str(time_el // 60) + " m " + str(time_el % 60) + " s"
	print(time_el_str)


	# write firing rates to file
	fout = open(os.path.join(output_dir, timestamp + '_firing_rates.txt'), 'w')
	fout.write("nu(A) = " + str(np.sum(counterA) / (num_periods_tot*period_duration) / core_size) + "\n")
	fout.write("nu(B) = " + str(np.sum(counterB) / (num_periods_tot*period_duration) / core_size) + "\n")
	fout.write("nu(C) = " + str(np.sum(counterC) / (num_periods_tot*period_duration) / core_size) + "\n")
	fout.write("nu(ctrl) = " + str(np.sum(counterCtrl) / (num_periods_tot*period_duration) / ctrl_size) + "\n")
	fout.write("core_size = " + str(core_size) + "\n")
	fout.write("ctrl_size = " + str(ctrl_size) + "\n")
	fout.close()

	return series

####################################
# transitionProbabilities
# Computes the likelihood of avalanche occurrence for each assembly, as well as and the likelihoods of "transitions"/triggering
# of assemblies, and saves the results to a file
# timestamp: the timestamp of the simulation data
# series: the time series as provided by timeSeries(...)
# output_dir: relative path to the output directory
def transitionProbabilities(timestamp, series, output_dir):

	np.seterr(divide='ignore', invalid='ignore')

	num_periods_tot = len(series)

	nA, nB, nC, nO, nN = 0, 0, 0, 0, 0 # frequencies of avalanches in different assemblies/overall
	nAA, nBB, nCC, nAB, nBA, nAC, nCA, nBC, nCB = 0, 0, 0, 0, 0, 0, 0, 0, 0 # frequencies of transitions
	nOO, nOA, nOB, nOC, nAO, nBO, nCO = 0, 0, 0, 0, 0, 0, 0 # frequencies of transitions involving "overall"
	nAN, nBN, nCN, nON, nNA, nNB, nNC, nNO, nNN = 0, 0, 0, 0, 0, 0, 0, 0, 0 # frequencies of transitions into/from void

	for i in range(num_periods_tot-1):

		if series[i] == "-": # there is no avalanche in this period
			nN += 1
			if series[i+1] == "-":
				nNN += 1
			else:
				nNO += 1
				if "A" in series[i+1]:
					nNA += 1
				if "B" in series[i+1]:
					nNB += 1
				if "C" in series[i+1]:
					nNC += 1

		else: # there is an avalanche in this period
			nO += 1
			if series[i+1] == "-":
				nON += 1
			else:
				nOO += 1
				if "A" in series[i+1]:
					nOA += 1
				if "B" in series[i+1]:
					nOB += 1
				if "C" in series[i+1]:
					nOC += 1

			if "A" in series[i]: # there is an avalanche in A
				nA += 1
				if series[i+1] == "-":
					nAN += 1
				else:
					nAO += 1
					if "A" in series[i+1]:
						nAA += 1
					if "B" in series[i+1]:
						nAB += 1
					if "C" in series[i+1]:
						nAC += 1						
			if "B" in series[i]: # there is an avalanche in B
				nB += 1
				if series[i+1] == "-":
					nBN += 1
				else:
					nBO += 1
					if "A" in series[i+1]:
						nBA += 1
					if "B" in series[i+1]:
						nBB += 1
					if "C" in series[i+1]:
						nBC += 1
			if "C" in series[i]: # there is an avalanche in C
				nC += 1
				if series[i+1] == "-":
					nCN += 1
				else:
					nCO += 1
					if "A" in series[i+1]:
						nCA += 1
					if "B" in series[i+1]:
						nCB += 1
					if "C" in series[i+1]:
						nCC += 1			

	# human-readable output
	fout = open(os.path.join(output_dir, timestamp + '_CA_probabilities.txt'), 'w')
	fout.write("Timestamp: " + timestamp)
	fout.write("\nTotal number of periods: " + str(num_periods_tot))
	fout.write("\n\nTotal probabilities:")
	fout.write("\np(A) = " + str(nA / num_periods_tot)) # likelihood of avalanche in A
	fout.write("\np(B) = " + str(nB / num_periods_tot)) # likelihood of avalanche in B
	fout.write("\np(C) = " + str(nC / num_periods_tot)) # likelihood of avalanche in C
	fout.write("\np(overall) = " + str(nO / num_periods_tot)) # probability of avalanche in overall population
	fout.write("\np(-) = " + str(nN / num_periods_tot)) # probability of no avalanche

	fout.write("\n\nTrigger probabilities:")
	fout.write("\np_A(A) = " + str(np.divide(nAA, num_periods_tot))) # likelihood of triggering A
	fout.write("\np_A(B) = " + str(np.divide(nAB, num_periods_tot))) # likelihood of triggering B
	fout.write("\np_A(C) = " + str(np.divide(nAC, num_periods_tot))) # likelihood of triggering C
	fout.write("\np_A(overall) = " + str(np.divide(nAO, nAO+nAN))) # probability of triggering something
	fout.write("\np_A(-) = " + str(np.divide(nAN, nAO+nAN))) # probability of triggering nothing
	
	fout.write("\np_B(A) = " + str(np.divide(nBA, num_periods_tot))) # likelihood of triggering A
	fout.write("\np_B(B) = " + str(np.divide(nBB, num_periods_tot))) # likelihood of triggering B
	fout.write("\np_B(C) = " + str(np.divide(nBC, num_periods_tot))) # likelihood of triggering C
	fout.write("\np_B(overall) = " + str(np.divide(nBO, nBO+nBN))) # probability of triggering something
	fout.write("\np_B(-) = " + str(np.divide(nBN, nBO+nBN))) # probability of triggering nothing

	fout.write("\np_C(A) = " + str(np.divide(nCA, num_periods_tot))) # likelihood of triggering A
	fout.write("\np_C(B) = " + str(np.divide(nCB, num_periods_tot))) # likelihood of triggering B
	fout.write("\np_C(C) = " + str(np.divide(nCC, num_periods_tot))) # likelihood of triggering C
	fout.write("\np_C(overall) = " + str(np.divide(nCO, nCO+nCN))) # probability of triggering something
	fout.write("\np_C(-) = " + str(np.divide(nCN, nCO+nCN))) # probability of triggering nothing

	fout.write("\np_overall(A) = " + str(np.divide(nOA, num_periods_tot))) # likelihood of triggering A
	fout.write("\np_overall(B) = " + str(np.divide(nOB, num_periods_tot))) # likelihood of triggering B
	fout.write("\np_overall(C) = " + str(np.divide(nOC, num_periods_tot))) # likelihood of triggering C
	fout.write("\np_overall(overall) = " + str(np.divide(nOO, nOO+nON))) # probability of triggering something
	fout.write("\np_overall(-) = " + str(np.divide(nON, nOO+nON))) # probability of triggering nothing

	fout.write("\np_-(A) = " + str(np.divide(nNA, num_periods_tot))) # likelihood of triggering A
	fout.write("\np_-(B) = " + str(np.divide(nNB, num_periods_tot))) # likelihood of triggering B
	fout.write("\np_-(C) = " + str(np.divide(nNC, num_periods_tot))) # likelihood of triggering C
	fout.write("\np_-(overall) = " + str(np.divide(nNO, nNO+nNN))) # probability of triggering something
	fout.write("\np_-(-) = " + str(np.divide(nNN, nNO+nNN))) # probability of triggering nothing
	fout.close()

	# output for facilitated machine readability
	fout = open(os.path.join(output_dir, timestamp + '_CA_probabilities_raw.txt'), 'w')
	fout.write(timestamp + "\n\n\n")
	fout.write(str(nA / num_periods_tot) + "\n")
	fout.write(str(nB / num_periods_tot) + "\n")
	fout.write(str(nC / num_periods_tot) + "\n")
	fout.write(str(nO / num_periods_tot) + "\n")
	fout.write(str(nN / num_periods_tot) + "\n")
	fout.close()

	n_transitions = sum((nOO, nON, nNO, nNN))
	if n_transitions == num_periods_tot-1:
		print("Normalization check suceeeded.")
	else:
		print("Normalization check failed:", n_transitions, "triggerings found, as compared to", num_periods_tot-1, "expected.")

####################################
# spikeRasterPlot
# Creates two spike raster plots in the data directory and copies them to the output directory
# timestamp: the timestamp of the data
# data_dir: the directory containing the simulation data
# output_dir: relative path to the output directory
# new_plots: specifies if new plots shall be created
def spikeRasterPlot(timestamp, data_dir, output_dir, new_plots):

	plot_file1 = timestamp + "_spike_raster.png"
	plot_file2 = timestamp + "_spike_raster2.png"

	work_dir = os.getcwd() # get the current working directory
	if data_dir == "":
		os.chdir(".")
	else:
		os.chdir(data_dir) # change to data directory

	if new_plots:
		fout = open("spike_raster.gpl", "w")
		fout.write("set term png enhanced font Sans 20 size 1280,960 lw 2.5\n")
		fout.write("set output '" + plot_file1 + "'\n\n")
		fout.write("Ne = " + str(exc_pop_size) + "\n")
		fout.write("Ni = " + str(inh_pop_size) + "\n")
		fout.write("set xlabel 'Time (s)'\n")
		fout.write("unset ylabel\n")
		fout.write("set yrange [0:1]\n")
		fout.write("set ytics out ('#0' 0.05, '#625' 0.23, '#1250' 0.41, '#1875' 0.59, '#2500' 0.77)\n")
		fout.write("plot [x=100:120] '" + timestamp + "_spike_raster.txt' using 1:($2 < Ne ? (0.9*$2/(Ne+Ni) + 0.05) : 1/0) notitle with dots lc 'blue', \\\n")
		fout.write("                 '" + timestamp + "_spike_raster.txt' using 1:($2 >= Ne ? (0.9*$2/(Ne+Ni) + 0.05) : 1/0) notitle with dots lc 'red'\n\n")
		fout.write("###########################################\n")
		fout.write("set output '" + plot_file2 + "'\n")
		fout.write("plot [x=100:180] '" + timestamp + "_spike_raster.txt' using 1:($2 < Ne ? (0.9*$2/(Ne+Ni) + 0.05) : 1/0) notitle with dots lc 'blue', \\\n")
		fout.write("                 '" + timestamp + "_spike_raster.txt' using 1:($2 >= Ne ? (0.9*$2/(Ne+Ni) + 0.05) : 1/0) notitle with dots lc 'red'\n")
		fout.close()

		os.system("gnuplot spike_raster.gpl")

	if os.path.exists(plot_file1) and os.path.exists(plot_file2):
		copy2(plot_file1, os.path.join(work_dir, output_dir)) # copy spike raster plot #1 to output directory
		copy2(plot_file2, os.path.join(work_dir, output_dir)) # copy spike raster plot #2 to output directory
	else:
		print("Warning: " + data_dir + ": plot files not found.")

	os.chdir(work_dir) # change back to previous working directory


######################################
# dirRecursion
# Walks recursively through a directory looking for spike raster data;
# if data are found, computes time series, avalanche likelihoods (and, if specified, creates spike raster plots)
# directory: the directory to consider
# output_dir: relative path to the output directory
# new_plots: specifies if new plots shall be created
def dirRecursion(directory, output_dir, new_plots):

	rawpaths = Path(directory)

	print("Reading directory " + directory)
	rawpaths = Path(directory)

	for x in sorted(rawpaths.iterdir()):

		dest_file = ""

		full_path = str(x)
		hpath = os.path.split(full_path)[0] # take head
		tpath = os.path.split(full_path)[1] # take tail

		if not x.is_dir():

			if "_spike_raster.txt" in tpath:
				timestamp = tpath.split("_spike_raster.txt")[0]
				series = timeSeries(timestamp, full_path, output_dir)
				transitionProbabilities(timestamp, series, output_dir)
				spikeRasterPlot(timestamp, hpath, output_dir, new_plots)

				params_file = os.path.join(hpath, timestamp + "_PARAMS.txt")
				if os.path.exists(params_file):
					copy2(params_file, output_dir)
				else:
					print("Warning: " + hpath + ": no parameter file found.")

				

		else:
			if hasTimestamp(tpath):
				dirRecursion(directory + os.sep + tpath, output_dir, new_plots)


###############################################
# main:

### example call from shell: python3 assemblyAvalancheStatistics.py "OVERLAP10 no AC, no ABC" 0.01 10 False 

if len(sys.argv) > 1: # if there is at least one additional commandline arguments
	paradigm = sys.argv[1]
	try:
		coreA, coreB, coreC = coreDefinitions(paradigm, core_size) # re-define cell assemblies
	except:
		raise
if len(sys.argv) > 2: # if there are at least 2 additional commandline arguments
	period_duration = float(sys.argv[2])
if len(sys.argv) > 3: # if there are at least 3 additional commandline arguments
	n_thresh = int(sys.argv[3])
if len(sys.argv) > 4: # if there are at least 4 additional commandline arguments
	if sys.argv[4] == "0" or sys.argv[4] == "False":
		new_plots = False
		print("Creation of new plots switched off")
	else:
		new_plots = True

output_dir = "./avalanche_statistics_" + str(period_duration) + "_" + str(n_thresh) # output directory for analysis results

if not os.path.exists(output_dir):
	os.mkdir(output_dir)

print("Output directory:", output_dir)
print("Paradigm:", paradigm)
print("Bin size:", str(period_duration), "s")
print("Detection threshold:", str(n_thresh))

dirRecursion('.', output_dir, new_plots) # walk through directories and analyze data
mergeRawData(output_dir, "_CA_probabilities_raw.txt", "all_trials_raw.txt", remove_raw=True) # merge machine-readable output

