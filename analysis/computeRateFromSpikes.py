######################################################################################################
###   Script to compute the firing rate over time via fixed time windows from spike raster data;   ###
###       computes the firing rates for all subpopulations in an "equal overlap" paradigm          ###
######################################################################################################

### Copyright 2020-2022 Jannik Luboeinski
### License: Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
from pathlib import Path
import os

np.set_printoptions(threshold=1e10, linewidth=200) # extend console print range for numpy arrays

# main parameters
Ne = 2500 # number of neurons in the excitatory population
Ni = 625 # number of neurons in the inhibitory population
Na = 600 # number of neurons in one cell assembly
overlap = 0.1 # relative size of the overlap
period_duration = 0.01 # binning period (in units of seconds)

# derived parameters
Nl = int(round(np.sqrt(Ne))) # the number of excitatory neurons in one line of a quadratic grid
No = int(np.round(overlap*Na)) # size of the overlap
No_half = int(np.round(No/2)) # half the size of the overlap
Na_wo_overlaps = Na-No-No_half
N = Ne + Ni # the number of neurons in the whole network

# extractRecursion
# Recursively looks for a data files, extracts spiking data from them and computes mean firing rates
# directory: the directory to look into
# fout: file handle to output file
def extractRecursion(directory, fout):

	data_found = False # specifies if any data has been found
	rawpaths = Path(directory)
	MI = []

	print("Contents of directory " + directory)
	print([str(x) for x in rawpaths.iterdir()])
	rawpaths = Path(directory)

	for x in sorted(rawpaths.iterdir()):

		full_path = str(x)
		(full_path_head, filename) = os.path.split(full_path)

		if x.is_file() and "_spike_raster.txt" in filename:

			data_found = True

			print("========================")
			print(filename + " in " + directory)
			print("------------------------")

			# read the last line and compute number of periods
			try:
				with open(full_path, 'rb') as f:
					f.seek(-2, os.SEEK_END)
					while f.read(1) != b'\n': # seek last line
						f.seek(-2, os.SEEK_CUR)
					last_line = f.readline().decode()
				num_periods_tot = int(float(last_line.split('\t\t')[0]) / period_duration) + 1
			except IOError:
				print('Error opening "' + filename + '"')
				exit()

			# allocate arrays for counting spikes
			tot_spikes = np.zeros(num_periods_tot, dtype=int) # number of total spikes for each timestep
			A_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in exclusive cell assembly ~A
			I_AB_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the intersection AB
			B_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in exclusive cell assembly ~B
			I_AC_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the intersection AC
			I_BC_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the intersection BC
			I_ABC_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the intersection ABC
			C_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in exclusive cell assembly ~C
			ctrl_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the exc. control population
			inh_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the inh. population
			tot_spikes = np.zeros(num_periods_tot, dtype=int) # number of spikes in the whole network

			# compute mean firing rates for different subpopulations
			print("I_AC: <", No_half, ", ~A: <", No_half+Na_wo_overlaps, ", I_ABC: <", No_half+Na_wo_overlaps+No_half, \
			      ", I_AB: <", Na, ", ~B: <", Na+Na_wo_overlaps, ", I_BC: <", Na+Na_wo_overlaps+No_half, \
			      ", ~C: <", Na+Na_wo_overlaps+No_half+Na_wo_overlaps)

			f = open(full_path)
			for line in f:
				segs = line.split('\t\t')

				if (segs[0] != ""):
					t = float(segs[0])
					n = int(segs[1])
					current_period = int(np.floor(t / period_duration))
					tot_spikes[current_period] += 1

					if n < No_half:
						I_AC_spikes[current_period] += 1
					elif n < No_half+Na_wo_overlaps:
						A_spikes[current_period] += 1
					elif n < No_half+Na_wo_overlaps+No_half:
						I_ABC_spikes[current_period] += 1
					elif n < Na:
						I_AB_spikes[current_period] += 1
					elif n < Na+Na_wo_overlaps:
						B_spikes[current_period] += 1
					elif n < Na+Na_wo_overlaps+No_half:
						I_BC_spikes[current_period] += 1
					elif n < Na+Na_wo_overlaps+No_half+Na_wo_overlaps:
						C_spikes[current_period] += 1
					elif n < Ne:
						ctrl_spikes[current_period] += 1
					elif n < N:
						inh_spikes[current_period] += 1
					else:
						print('Error reading neuron from "' + filename + '"')
				else:
					print('Error reading time from "' + filename + '"')

			f.close()

			# write mean firing rate data to file
			for i in range(num_periods_tot):
				mean_A = A_spikes[i] / Na_wo_overlaps / period_duration
				mean_I_AB = I_AB_spikes[i] / No_half / period_duration
				mean_B = B_spikes[i] / Na_wo_overlaps / period_duration
				mean_I_AC = I_AC_spikes[i] / No_half / period_duration
				mean_I_BC = I_BC_spikes[i] / No_half / period_duration
				mean_I_ABC = I_ABC_spikes[i] / No_half / period_duration
				mean_C = C_spikes[i] / Na_wo_overlaps / period_duration
				mean_ctrl = ctrl_spikes[i] / (Ne-Na-Na_wo_overlaps-Na_wo_overlaps-No_half) / period_duration
				mean_inh = inh_spikes[i] / Ni / period_duration
				mean_tot = tot_spikes[i] / N / period_duration

				# write time (at 1/2 of a period) and values for this period
				fout.write(str(round((i+0.5)*period_duration,4)) + "\t\t" + \
				           str(mean_A) + "\t\t" + str(mean_I_AB) + "\t\t" + str(mean_B) + "\t\t" + \
						   str(mean_I_AC) + "\t\t" + str(mean_I_BC) + "\t\t" + str(mean_I_ABC) + "\t\t" + str(mean_C) + "\t\t" + \
				           str(mean_ctrl) + "\t\t" + str(mean_inh) + "\t\t" + str(mean_tot) + "\n")


		elif x.is_dir():
			ret = extractRecursion(full_path, fout)
			data_found = data_found or ret

	return data_found

try:
	fout = open("rates_over_time.txt", "w")

except IOError:
	print('Error opening "rates_over_time.txt"')
	exit()

if extractRecursion('.', fout):
	print("========================")

fout.close()
