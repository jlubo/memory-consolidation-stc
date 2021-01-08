###########################################################################################
### Script to extract the mean correlation between neuronal subpopulations of a network ###
###                             from spike raster data                                  ###
###########################################################################################

### Copyright 2020-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
from pathlib import Path

# size of the excitatory population
Ne = 1600

 # size of one cell assembly
core_size = 150

# sizes of the (sub)populations
Nas = core_size // 2 # assembly neurons receiving recall stimulation
Nans = core_size // 2 # assembly neurons not receiving recall stimulation
Nc = Ne - Nas - Nans # control neurons

# firingRateMatrix
# timestamp: timestamp string that labels the data set
# spike_raster_file: name of the file containg the spike raster data
# t_prime: time for which the firing rate matrix shall be computed
# t_wfr: duration of the sliding window for computing the firing rate
def firingRateMatrix(timestamp, spike_raster_file, t_prime, t_wfr):

	# determine the start and end of the effective sliding window
	if t_prime < t_wfr/2.:
		t_wfr_start = 0
		t_wfr_end = t_prime + t_wfr/2.
	else:
		t_wfr_start = t_prime - t_wfr/2.
		t_wfr_end = t_prime + t_wfr/2.

	# read data
	f = open(spike_raster_file)
	rawdata = f.read()
	f.close()

	rawdata = rawdata.split('\n')
	num_rows = len(rawdata)-1
	t_max = np.int(np.double(rawdata[num_rows-1].split('\t\t')[0]))

	count_e = np.zeros(Ne, dtype=np.int) # neuronal spike counts of the whole excitatory population

	# go through the spike raster data and find for every neuron the spikes occurring in the effective window
	for i in range(num_rows):
		row = rawdata[i].split('\t\t')
		t = np.double(row[0])
		n = np.int(row[1])

		if t >= t_wfr_start and t < t_wfr_end and n < Ne: # check if effective time window is being hit and only consider exc. neurons

			count_e[n] += 1

	# average over time
	fr = count_e / (t_wfr_end - t_wfr_start)

	# write FR data to file
	filename = timestamp + '_fr_' + ('%.1f' % t_prime) + '.txt'
	Nl = int(round(np.sqrt(Ne)))
	#np.savetxt(filename, np.reshape(fr, (Nl, Nl)))

	return fr

# corrMatrix
# timestamp: timestamp string that labels the data set
# fr_mat1: firing rate matrix of the first pattern
# fr_mat2: firing rate matrix of the second pattern
# suffix: suffix for the filename
def corrMatrix(timestamp, fr_mat1, fr_mat2, suffix):

	cm = np.zeros((Ne,Ne), dtype=np.float32)

	for i in range(Ne):
		for j in range(Ne):
			cm[i][j] = (fr_mat1[i] * fr_mat2[j])

	# write correlation data to file
	#filename = timestamp + '_corr_' + suffix + '.txt'
	#np.savetxt(filename, cm)

	return cm


####################################
# corrBetweenSubPop
# timestamp: timestamp string that labels the data set
# cm: the correlation matrix
# suffix: suffix for the filename
def corrBetweenSubPop(timestamp, cm, suffix):

	# integrated correlation means and standard deviations for the subpopulations
	as_within, as_within_sd = 0, 0
	as_to_ans, as_to_ans_sd = 0, 0
	as_to_c, as_to_c_sd = 0, 0
	ans_within, ans_within_sd = 0, 0
	ans_to_as, ans_to_as_sd = 0, 0
	ans_to_c, ans_to_c_sd = 0, 0
	c_within, c_within_sd = 0, 0
	c_to_as, c_to_as_sd = 0, 0
	c_to_ans, c_to_ans_sd = 0, 0

	# compute mean correlation values
	for i in range(Nas): # loop over 'as' neurons

		for j in range(Nas):
			as_within += cm[i][j]

		for j in range(Nas, Nas+Nans):
			as_to_ans += cm[i][j]

		for j in range(Nas+Nans, Ne):
			as_to_c += cm[i][j]

	for i in range(Nas, Nas+Nans): # loop over 'ans' neurons

		for j in range(Nas):
			ans_to_as += cm[i][j]

		for j in range(Nas, Nas+Nans):
			ans_within += cm[i][j]

		for j in range(Nas+Nans, Ne):
			ans_to_c += cm[i][j]

	for i in range(Nas+Nans, Ne): # loop over 'ctrl' neurons

		for j in range(Nas):
			c_to_as += cm[i][j]

		for j in range(Nas, Nas+Nans):
			c_to_ans += cm[i][j]

		for j in range(Nas+Nans, Ne):
			c_within += cm[i][j]

	as_within /= Nas**2
	as_to_ans /= (Nas*Nans)
	as_to_c /= (Nas*Nc)
	ans_within /= Nans**2
	ans_to_as /= (Nans*Nas)
	ans_to_c /= (Nans*Nc)
	c_within /= Nc**2
	c_to_as /= (Nc*Nas)
	c_to_ans /= (Nc*Nans)

	# compute standard devations of the correlation values
	for i in range(Nas): # loop over 'as' neurons

		for j in range(Nas):
			as_within_sd += (as_within - cm[i][j])**2

		for j in range(Nas, Nas+Nans):
			as_to_ans_sd += (as_to_ans - cm[i][j])**2

		for j in range(Nas+Nans, Ne):
			as_to_c_sd += (as_to_c - cm[i][j])**2

	for i in range(Nas, Nas+Nans): # loop over 'ans' neurons

		for j in range(Nas):
			ans_to_as_sd += (ans_to_as - cm[i][j])**2

		for j in range(Nas, Nas+Nans):
			ans_within_sd += (ans_within - cm[i][j])**2

		for j in range(Nas+Nans, Ne):
			ans_to_c_sd += (ans_to_c - cm[i][j])**2

	for i in range(Nas+Nans, Ne): # loop over 'ctrl' neurons

		for j in range(Nas):
			c_to_as_sd += (c_to_as - cm[i][j])**2

		for j in range(Nas, Nas+Nans):
			c_to_ans_sd += (c_to_ans - cm[i][j])**2

		for j in range(Nas+Nans, Ne):
			c_within_sd += (c_within - cm[i][j])**2

	as_within_sd = np.sqrt(as_within_sd / Nas**2)
	as_to_ans_sd = np.sqrt(as_to_ans_sd / (Nas*Nans))
	as_to_c_sd = np.sqrt(as_to_c_sd / (Nas*Nc))
	ans_within_sd = np.sqrt(ans_within_sd / Nans**2)
	ans_to_as_sd = np.sqrt(ans_to_as_sd / (Nans*Nas))
	ans_to_c_sd = np.sqrt(ans_to_c_sd / (Nans*Nc))
	c_within_sd = np.sqrt(c_within_sd / Nc**2)
	c_to_as_sd = np.sqrt(c_to_as_sd / (Nc*Nas))
	c_to_ans_sd = np.sqrt(c_to_ans_sd / (Nc*Nans))

	# write averaged correlation data to file
	filename = 'corr_subpop_' + suffix + '.txt'
	#"\"as→as\"\t\terr\t\t\"as→ans\"\t\terr\t\t\"as→c\"\t\terr\t\t\"ans→ans\"\t\terr\t\t\"ans→as\"\t\terr\t\t\"ans→c\"\t\terr\t\t\"c→c\"\t\terr\t\t\"c→as\"\t\terr\t\t\"c→ans\"\t\terr"

	f = open(filename, 'a', encoding="utf-8")
	f.write(str(as_within) + "\t\t" + str(as_within_sd) + "\t\t" + \
	        str(as_to_ans) + "\t\t" + str(as_to_ans_sd) + "\t\t" + \
	        str(as_to_c) + "\t\t" + str(as_to_c_sd) + "\t\t" + \
	        str(ans_within) + "\t\t" + str(ans_within_sd) + "\t\t" + \
	        str(ans_to_as) + "\t\t" + str(ans_to_as_sd) + "\t\t" + \
	        str(ans_to_c) + "\t\t" + str(ans_to_c_sd) + "\t\t" + \
	        str(c_within) + "\t\t" + str(c_within_sd) + "\t\t" + \
	        str(c_to_as) + "\t\t" + str(c_to_as_sd) + "\t\t" + \
	        str(c_to_ans) + "\t\t" + str(c_to_ans_sd) + "\n")
	f.close()

# main
rawpaths = Path('Standby')
for x in rawpaths.iterdir():
	if "_spike_raster.txt" in str(x) and not x.is_dir():
		spike_raster_file = str(x)
		timestamp = spike_raster_file.split("_spike_raster.txt")[0]

		fr = firingRateMatrix(timestamp, spike_raster_file, 7.0, 0.5)
		cm = corrMatrix(timestamp, fr, fr, "Standby")
		corrBetweenSubPop(timestamp, cm, "Standby")

rawpaths = Path('Control stim')
for x in rawpaths.iterdir():
	if "_spike_raster.txt" in str(x) and not x.is_dir():
		spike_raster_file = str(x)
		timestamp = spike_raster_file.split("_spike_raster.txt")[0]

		fr = firingRateMatrix(timestamp, spike_raster_file, 7.0, 0.5)
		cm = corrMatrix(timestamp, fr, fr, "Control stim")
		corrBetweenSubPop(timestamp, cm, "Control stim")

rawpaths = Path('10s-recall')
for x in rawpaths.iterdir():
	if "_spike_raster.txt" in str(x) and not x.is_dir():
		spike_raster_file = str(x)
		timestamp = spike_raster_file.split("_spike_raster.txt")[0]

		fr = firingRateMatrix(timestamp, spike_raster_file, 20.0, 0.5)
		cm = corrMatrix(timestamp, fr, fr, "10s-recall")
		corrBetweenSubPop(timestamp, cm, "10s-recall")

rawpaths = Path('8h-recall')
for x in rawpaths.iterdir():
	if "_spike_raster.txt" in str(x) and not x.is_dir():
		spike_raster_file = str(x)
		timestamp = spike_raster_file.split("_spike_raster.txt")[0]

		fr = firingRateMatrix(timestamp, spike_raster_file, 28810.0, 0.5)
		cm = corrMatrix(timestamp, fr, fr, "8h-recall")
		corrBetweenSubPop(timestamp, cm, "8h-recall")
