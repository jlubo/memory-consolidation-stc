###############################################################################################
### Reads out Q and MI values (mean and standard deviation) from data files and plots them, ###
### individually for 10s- and 8h-recall.                                                    ###
###############################################################################################

### Copyright 2023-2024 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import scipy.stats as st

Q_COL = 1 # data column in which the mean of Q is found
MI_COL = 3 # data column in which the mean of the mutual information is found
TRIALS_COL = 7 # data column in which the number of trials is found

# plotQMIGraphs
# Plots Q and MI values over assembly core size. Creates figures with up to two graphs,
# one for 10s- and one for 8h-recall.
# data_file_10s: file containing the data of 10s-recall
# data_file_8h: file containing the data of 8h-recall
# show [optional]: specifies whether to show the plots immediately or not (mainly intended for Jupyter)
# store [optional]: specifies whether to store the plots
def plotQMIGraphs(data_file_10s, data_file_8h, show = False, store = True):

	# Loop over measures
	for msr, unit, col in [("Pattern completion coefficient", "", Q_COL), ("Mutual information", "bits", MI_COL)]:

		# Prepare for the plotting
		plt.rcParams.update({'font.size': 16})
		fig, ax = plt.subplots()
		data_files = [data_file_10s, data_file_8h]
		data_labels = ["10s", "8h"]
		colors = ["#bb1111", "#1111bb"]
		markers = ["o", "s"]

		# Loop over paradigms
		for data_file, data_label, color, marker in zip(data_files, data_labels, colors, markers):

			# Test if a data file has been provided, otherwise, skip this paradigm
			if not data_file:
				continue

			# Load NumPy array from data file
			data = np.loadtxt(data_file)

			# Create a DataFrame, and sort
			data_mean_sd = pd.DataFrame(data).sort_values(by=[0,1])

			# Get the set of available assembly core sizes
			core_sizes = np.int32(sorted(list(set(data_mean_sd.iloc[:,0]))))
			print(f"Creating graph plot of {data_label}-recall over the following core sizes: {core_sizes}.")

			# Determine the number of trials and the factor for 95% compatibility interval
			n_trials = int(data[0,TRIALS_COL])
			print(f"Assumed (read out) number of trials: {n_trials}.")
			ci95_factor = st.t.ppf(1-0.05/2, df=n_trials-1) / np.sqrt(n_trials)
			#ci95_factor = st.t.interval(0.95, n_trials-1, loc=0, scale=1/np.sqrt(n_trials))
			#ci95_factor = 1.98 / np.sqrt(100)  # 100 trials
			#ci95_factor = 2.01 / np.sqrt(50)  # 50 trials

			# Plot the data
			selected_data_mean = data_mean_sd.iloc[:,col]
			selected_data_ci95 = data_mean_sd.iloc[:,col+1]*ci95_factor
			x_range = np.arange(len(core_sizes))
			ax.plot(x_range, selected_data_mean, marker=marker, label=f"{data_label}-recall", linestyle='None', color=color)
			ax.errorbar(x_range, selected_data_mean, yerr=selected_data_ci95, fmt='none', ecolor=color, capsize=4, capthick=1.5)

		# Finalize the figure
		ax.set_xlabel("Pattern size (#cells)")
		ax.set_xticks(x_range)
		ax.set_xticklabels(core_sizes)
		if msr == "Pattern completion coefficient":
			y_min = 0.00
			y_max = 0.10
			y_range = np.append(np.arange(y_min, y_max, 0.01), y_max)
		elif msr == "Mutual information":
			y_min = 0.0
			y_max = 1.4
			y_range = np.append(np.arange(y_min, y_max, 0.2), y_max)
		#ax.set_ylim(y_min, y_max)
		#ax.set_yticks(y_range)
		if unit:
			ax.set_ylabel(f'{msr} ({unit})')
		else:
			ax.set_ylabel(msr)
		#ax.minorticks_on()
		ax.grid(True, axis="y", linestyle = '-', alpha=0.6, which="major")
		ax.grid(True, axis="y", linestyle = '--', alpha=0.3, which="minor")
		ax.legend()

		plt.tight_layout()
		if show:
			plt.show()
		if store:
			plt.savefig(f"{msr} (graph plot).svg")

# plotQMIBarGroups
# Plots Q and MI values over assembly core size. Creates figures with groups of bar plots,
# each group showing the data for 10s- and 8h-recall.
# data_file_10s: file containing the data of 10s-recall
# data_file_8h: file containing the data of 8h-recall
# show [optional]: specifies whether to show the plots immediately or not (mainly intended for Jupyter)
# store [optional]: specifies whether to store the plots
def plotQMIBarGroups(data_file_10s, data_file_8h, show = False, store = True):

	# Labels of paradigms to consider
	prdgms = ["10s", "8h"]

	# Color scheme for paradigms
	colors = ["#bb1111", "#1111bb"]

	# Loop over measures
	for msr, unit, col in [("Pattern completion coefficient", "", Q_COL), ("Mutual information", "bits", MI_COL)]:

		# Load NumPy array from data file, create a DataFrame, and sort
		data_10s = pd.read_table(data_file_10s, header=None, sep=" ")
		data_10s = data_10s.assign(prdgm = pd.Series([prdgms[0]]*data_10s.shape[0]).values)
		data_8h = pd.read_table(data_file_8h, header=None, sep=" ")
		data_8h = data_8h.assign(prdgm = pd.Series([prdgms[1]]*data_8h.shape[0]).values)
		data = pd.concat((data_10s, data_8h))
		data_mean_sd = pd.DataFrame(data).sort_values(by=[0, "prdgm"])

		# Remove the rows corresponding core size < 100 and > 250
		data_mean_sd.drop(data_mean_sd[data_mean_sd.iloc[:,0] < 100].index, inplace=True)
		data_mean_sd.drop(data_mean_sd[data_mean_sd.iloc[:,0] > 250].index, inplace=True)

		# Get the set of available assembly core sizes
		core_sizes = np.int32(sorted(list(set(data_mean_sd.iloc[:,0]))))
		print(f"Creating bar group plot of {prdgms[0]}- and {prdgms[1]}-recall "
		      f"over the following core sizes: {core_sizes}.")

		# Determine the number of trials and the factor for 95% compatibility interval
		n_trials = int(data.iloc[0,TRIALS_COL])
		print(f"Assumed (read out) number of trials: {n_trials}.")
		ci95_factor = st.t.ppf(1-0.05/2, df=n_trials-1) / np.sqrt(n_trials)
		#ci95_factor = 1.98 / np.sqrt(100)  # 100 trials

		# Prepare for the plotting
		plt.rcParams.update({'font.size': 16})
		x_range = np.arange(len(core_sizes))
		x2_range = np.arange(len(prdgms))
		width = 0.2
		fig, ax = plt.subplots()

		# Do the plotting (inspired by https://matplotlib.org/stable/gallery/lines_bars_and_markers/barchart.html)
		for x2, prdgm, color in zip(x2_range, prdgms, colors):
			offset = width * x2
			selected_data_mean = data_mean_sd[data_mean_sd["prdgm"] == prdgm].iloc[:,col]
			selected_data_ci95 = data_mean_sd[data_mean_sd["prdgm"] == prdgm].iloc[:,col+1]*ci95_factor
			bars = ax.bar(x_range + offset, 
				          selected_data_mean,
				          yerr=selected_data_ci95,
				          width=width, 
				          label=f"{prdgm}-recall", 
				          color=color)
			#ax.bar_label(bars, padding=2.5)
			# Separately plot errorbars
			ax.errorbar(x_range + offset, selected_data_mean, yerr=selected_data_ci95, fmt='none', ecolor='darkgray', capsize=4, capthick=1.5)
		ax.set_xlabel("Pattern size (#cells)")
		ax.set_xticks(x_range + width/2)
		ax.set_xticklabels(core_sizes)
		if msr == "Pattern completion coefficient":
			y_min = 0.00
			y_max = 0.09
			y_range = np.append(np.arange(y_min, y_max, 0.01), y_max)
		elif msr == "Mutual information":
			y_min = 0.0
			y_max = 1.4
			y_range = np.append(np.arange(y_min, y_max, 0.2), y_max)
		#ax.set_ylim(y_min, y_max)
		#ax.set_yticks(y_range)
		if unit:
			ax.set_ylabel(f'{msr} ({unit})')
		else:
			ax.set_ylabel(msr)
		#ax.minorticks_on()
		ax.grid(True, axis="y", linestyle = '-', alpha=0.6, which="major")
		ax.grid(True, axis="y", linestyle = '--', alpha=0.3, which="minor")
		ax.legend()

		plt.tight_layout()
		if show:
			plt.show()
		if store:
			plt.savefig(f"{msr} (bar plot).svg")

# main
if __name__ == '__main__':
	plotQMIGraphs("./data_Q_MI_10s.txt",
	              "./data_Q_MI_8h.txt",
	              show = False)
	plotQMIBarGroups("./data_Q_MI_10s.txt",
	                 "./data_Q_MI_8h.txt",
	                 show = False)
