################################################################################################################
### The main purpose of this script is to facilitate the comparison of data from two different simulators.   ###
### It allows to plot the mean and the s.e.m. across `M` batches of the mean _and_ the standard deviation of ###
### the synaptic weight, membrane potential, membrane current, and calcium/protein concentration across `N`  ###
### trials. Nonetheless, the script also provides functions to plot data from deterministic processes.       ###
################################################################################################################

### Copyright 2022-2024 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import matplotlib.pyplot as plt
import json
import os
from pathlib import Path
from utilityFunctions import hasTimestamp

#####################################
# plotMeanAndSEMofX
# Plots the mean and the s.e.m. across `n_batches` for given values of synaptic weight, membrane voltage, and calcium/protein concentration. Deterministic data
# can be plotted as well, corresponding to `n_batches = 1`, with no s.e.m. considered.
# - X_name: the name of the measure whose mean and s.e.m. is to be plotted
# - n_batches: number of batches (= 1 implies there is no s.e.m.)
# - config: configuration parameters in JSON format
# - framework_name_1: name of the first framework that is to be compared (may be shown alone as well)
# - data_stacked_1: two-dimensional array containing the values of the membrane potential, weights, calcium and protein concentration, etc. over time
# - spikes_stacked_1: two-dimensional array containing the neuron and time of each spike
# - framework_name_2 [optional]: name of the second framework that is to be compared (may be left out)
# - data_stacked_2 [optional]: two-dimensional array containing the values of the membrane potential, weights, calcium and protein concentration, etc. over time (may be left out)
# - spikes_stacked_2 [optional]: two-dimensional array containing the neuron and time of each spike (may be left out)
# - X_cols [optional]: dictionary specifying the data columns of the measures whose mean and s.e.m. is to be plotted (must point to the mean, the std.dev. is sought in the next column;
#                      if a value is <= 0, the related variable is not plotted)
# - X_cols_2 [optional]: like `X_cols`, but for data from second framework
# - spike_neuron [optional]: index of the neuron the spikes of which shall be plotted
# - t_max [optional]: maximum time for plotting the data points (in ms, `None` means all data will be plotted)
# - store_path [optional]: path to store resulting graphics file
# - sampling [optional]: sampling (1: best resolution, >1: samples are left out) -- can be useful to reduce file sizes
# - wls [optional]: shift for better visualization of late-phase weight
# - time_unit_1 [optional]: unit of time of the first framework (1: millisecond, 1000: second, ...)
# - time_unit_2 [optional]: unit of time of the second framework (1: millisecond, 1000: second, ...)
def plotMeanAndSEMofX(X_name, n_batches, config, 
					  framework_name_1, data_stacked_1, spikes_stacked_1, 
					  framework_name_2 = "", data_stacked_2 = None, spikes_stacked_2 = None,
                      X_cols = {"voltage": 1, "current": -1, "weight": 3, "weight_l": 4, "calcium": 5, "protein": 7}, 
					  X_cols_2 = {"voltage": 1, "current": -1, "weight": 3, "weight_l": 4, "calcium": 5, "protein": 7},
					  spike_neuron = None, t_max = None, store_path = 'figure.svg', sampling = 1, wls = 1,
					  time_unit_1 = 1, time_unit_2 = 1):

	FIGSIZE = 8 # fonts inversely scale with this 
	h_0 = config["synapses"]["h_0"]
	fig, axes = plt.subplots(nrows=3, ncols=1, sharex=False, figsize=(FIGSIZE, FIGSIZE))

	# downsampling with integer sampling variable
	if type(sampling) == int and sampling > 1:
		data_stacked_1 = data_stacked_1[np.arange(len(data_stacked_1[:,0])) % sampling == 0] # sample only every n-th value
		if data_stacked_2 is not None:
			data_stacked_2 = data_stacked_2[np.arange(len(data_stacked_2[:,0])) % sampling == 0] # sample only every n-th value

	# time unit scaling of certain data (if necessary) and adjustment of time unit description
	def check_timescale(timescale):
		if t_max and t_max > timescale:
			return True
		elif data_stacked_1[:,0][-1]*time_unit_1 > timescale:
			return True
		return False
	def time_unit_conversion(factor = 1):
		nonlocal data_stacked_1, spikes_stacked_1
		nonlocal data_stacked_2, spikes_stacked_2
		nonlocal t_max
		# only do the rescaling of data from framework #1 if the effective factor is != 1
		if time_unit_1 * factor != 1:
			if data_stacked_1 is not None:
				data_stacked_1 = np.copy(data_stacked_1)
				data_stacked_1[:,0] *= time_unit_1 * factor
			if spikes_stacked_1 is not None:
				spikes_stacked_1 = np.copy(spikes_stacked_1)
				spikes_stacked_1[0] *= time_unit_1 * factor
		# only do the rescaling of data from framework #2 if the effective factor is != 1
		if time_unit_2 * factor != 1:
			if data_stacked_2 is not None:
				data_stacked_2 = np.copy(data_stacked_2)
				data_stacked_2[:,0] *= time_unit_2 * factor
			if spikes_stacked_2 is not None:
				spikes_stacked_2 = np.copy(spikes_stacked_2)
				spikes_stacked_2[0] *= time_unit_2 * factor
		# adjust the maximum time
		if t_max:
			t_max *= factor
	if check_timescale(7.2e6):  # more than 2 h are given
		time_unit_str = "h"
		time_unit_conversion(1/3.6e6)
	elif check_timescale(6e4):  # more than 1 min is given
		time_unit_str = "min"
		time_unit_conversion(1/6e4)
		t_max 
	elif check_timescale(2e3):  # more than 2 s are given
		time_unit_str = "s"
		time_unit_conversion(1/1e3)
	else:
		time_unit_str = "ms"
		time_unit_conversion()

	#----------------------------------
	# Weight dynamics
	if X_cols.get("weight", 0) > 0:
		# set axis labels for axes[0]
		axes[0].set_xlabel(f"Time ({time_unit_str})")
		if X_name != "":
			axes[0].set_ylabel(f"{X_name} of synaptic weight (%)")
		else:
			axes[0].set_ylabel("Synaptic weight (%)")
		
		# plot data for axes[0]
		axes[0].plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["weight"]]/h_0*100, color="#800000", label=framework_name_1, marker='None', zorder=8) # early phase
		if n_batches > 1:
			axes[0].fill_between(data_stacked_1[:,0], (data_stacked_1[:,X_cols["weight"]]-data_stacked_1[:,X_cols["weight"]+1]/np.sqrt(n_batches))/h_0*100, 
													  (data_stacked_1[:,X_cols["weight"]]+data_stacked_1[:,X_cols["weight"]+1]/np.sqrt(n_batches))/h_0*100, color="#800000", alpha=0.5)
		if X_cols.get("weight_l", 0) > 0:
			axes[0].plot(data_stacked_1[:,0], (data_stacked_1[:,X_cols["weight_l"]]/h_0+wls)*100, color="#1f77b4", label=framework_name_1+" (L)", marker='None', zorder=7) # late phase
			if n_batches > 1:
				axes[0].fill_between(data_stacked_1[:,0], ((data_stacked_1[:,X_cols["weight_l"]]-data_stacked_1[:,X_cols["weight_l"]+1]/np.sqrt(n_batches))/h_0+wls)*100, 
														  ((data_stacked_1[:,X_cols["weight_l"]]+data_stacked_1[:,X_cols["weight_l"]+1]/np.sqrt(n_batches))/h_0+wls)*100, color="#1f77b4", alpha=0.5)
		if data_stacked_2 is not None:
			axes[0].plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["weight"]]/h_0*100, color="#330000", linestyle='dashed', label=framework_name_2, marker='None', zorder=10)
			if n_batches > 1:
				axes[0].fill_between(data_stacked_2[:,0], (data_stacked_2[:,X_cols_2["weight"]]-data_stacked_2[:,X_cols_2["weight"]+1]/np.sqrt(n_batches))/h_0*100, 
														  (data_stacked_2[:,X_cols_2["weight"]]+data_stacked_2[:,X_cols_2["weight"]+1]/np.sqrt(n_batches))/h_0*100, color="#330000", alpha=0.5)
			if X_cols_2.get("weight_l", 0) > 0:
				axes[0].plot(data_stacked_2[:,0], (data_stacked_2[:,X_cols_2["weight_l"]]/h_0+wls)*100, color="#082133", linestyle='dashed', label=framework_name_2+" (L)", marker='None', zorder=7) # late phase
				if n_batches > 1:
					axes[0].fill_between(data_stacked_2[:,0], ((data_stacked_2[:,X_cols_2["weight_l"]]-data_stacked_2[:,X_cols_2["weight_l"]+1]/np.sqrt(n_batches))/h_0+wls)*100, 
															  ((data_stacked_2[:,X_cols_2["weight_l"]]+data_stacked_2[:,X_cols_2["weight_l"]+1]/np.sqrt(n_batches))/h_0+wls)*100, color="#082133", alpha=0.5)
			
		# create axis ranges and legend for axes[0]
		if t_max:
			axes[0].set_xlim(0, t_max)
		#axes[0].set_ylim(0,1)
		axes[0].legend()

	#----------------------------------
	# Membrane voltage
	if X_cols.get("voltage", 0) > 0:

		# membrane current on twin axis
		if X_cols.get("current", 0) > 0:
			ax1twin = axes[1].twinx()

		# set axis labels for axes[1]
		axes[1].set_xlabel(f"Time ({time_unit_str})")
		if X_name != "":
			axes[1].set_ylabel(f"{X_name} of membrane potential (mV)")
			if X_cols.get("current", 0) > 0:
				ax1twin.set_ylabel(f"{X_name} of membrane current (nA)")
		else:
			axes[1].set_ylabel("Membrane potential (mV)")
			if X_cols.get("current", 0) > 0:
				ax1twin.set_ylabel("Membrane current (nA)")
		
		# plot data for axes[1]
		if X_cols.get("current", 0) > 0:
			axes[1].plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["voltage"]], color="#ff0000", label=f'Membrane pot. ({framework_name_1})', marker='None', zorder=10)
			ax1twin.plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["current"]], color="#ffee00", label=f'Membrane curr. ({framework_name_1})', marker='None', zorder=10)
		else:
			axes[1].plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["voltage"]], color="#ff0000", label=framework_name_1, marker='None', zorder=10)
		if n_batches > 1:
			axes[1].fill_between(data_stacked_1[:,0], data_stacked_1[:,X_cols["voltage"]]-data_stacked_1[:,X_cols["voltage"]+1]/np.sqrt(n_batches), 
													  data_stacked_1[:,X_cols["voltage"]]+data_stacked_1[:,X_cols["voltage"]+1]/np.sqrt(n_batches), color="#ff0000", alpha=0.5)
		if data_stacked_2 is not None:
			if X_cols_2.get("current", 0) > 0:
				axes[1].plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["voltage"]], color="#330000", linestyle='dashed', label=f'Membrane pot. ({framework_name_2})', marker='None', zorder=10)
				ax1twin.plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["current"]], color="#332f00", linestyle='dashed', label=f'Membrane curr. ({framework_name_2})', marker='None', zorder=10)
			else:
				axes[1].plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["voltage"]], color="#330000", linestyle='dashed', label=framework_name_2, marker='None', zorder=10)
			if n_batches > 1:
				axes[1].fill_between(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["voltage"]]-data_stacked_2[:,X_cols_2["voltage"]+1]/np.sqrt(n_batches), 
														  data_stacked_2[:,X_cols_2["voltage"]]+data_stacked_2[:,X_cols_2["voltage"]+1]/np.sqrt(n_batches), color="#330000", alpha=0.5)

		# plot spikes of specified neuron
		if spike_neuron is not None:
			marker_type = '.' # ',' 
			marker_size = 2
			mask_1 = (spikes_stacked_1[1] == spike_neuron)
			axes[1].plot(spikes_stacked_1[0][mask_1], -90*np.ones(np.sum(mask_1)), marker_type, color="#ff0000", markersize=marker_size)
			if spikes_stacked_2 is not None:
				mask_2 = (spikes_stacked_2[1] == spike_neuron)
				axes[1].plot(spikes_stacked_2[0][mask_2], -91*np.ones(np.sum(mask_2)), marker_type, color="#330000", markersize=marker_size)
	
		# create axis ranges and legend for `axes[0]` (and possibly for `ax1twin``)
		if t_max:
			axes[1].set_xlim(0, t_max)
		#axes[1].set_ylim(0,0.006)
		if X_cols.get("current", 0) > 0:
			ax1twin.set_ylim(-2, 4)
			handles, labels = axes[1].get_legend_handles_labels()
			handles_twin, labels_twin = ax1twin.get_legend_handles_labels()
			axes[1].legend(handles + handles_twin, labels + labels_twin, loc="center left")
		else:
			axes[1].legend()

	# ---------------------------------
	# show calcium dynamics if considering timescales of milliseconds or seconds, otherwise show protein dynamics
	if time_unit_str in ["ms", "s"]:

	#----------------------------------
	# Calcium dynamics
		if X_cols.get("calcium", 0) > 0:
			# set axis labels for `axes[2]`
			axes[2].set_xlabel(f"Time ({time_unit_str})")
			if X_name != "":
				axes[2].set_ylabel(f"{X_name} of calcium concentration")
			else:
				axes[2].set_ylabel("Calcium concentration")
			
			# plot data for `axes[2]`
			axes[2].plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["calcium"]], color="#c8c896", label=framework_name_1, marker='None', zorder=8)
			if n_batches > 1:
				axes[2].fill_between(data_stacked_1[:,0], data_stacked_1[:,X_cols["calcium"]]-data_stacked_1[:,X_cols["calcium"]+1]/np.sqrt(n_batches), 
														data_stacked_1[:,X_cols["calcium"]]+data_stacked_1[:,X_cols["calcium"]+1]/np.sqrt(n_batches), color="#c8c896", alpha=0.5)
			if data_stacked_2 is not None:
				axes[2].plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["calcium"]], color="#333326", linestyle='dashed', label=framework_name_2, marker='None', zorder=8)
				if n_batches > 1:
					axes[2].fill_between(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["calcium"]]-data_stacked_2[:,X_cols_2["calcium"]+1]/np.sqrt(n_batches), 
															data_stacked_2[:,X_cols_2["calcium"]]+data_stacked_2[:,X_cols_2["calcium"]+1]/np.sqrt(n_batches), color="#333326", alpha=0.5)
			axes[2].axhline(y=config["synapses"]["early_phase"]["theta_p"], label='LTP threshold', linestyle='dashed', color="#969664", zorder=7)
			axes[2].axhline(y=config["synapses"]["early_phase"]["theta_d"], label='LTD threshold', linestyle='dashed', color="#969696", zorder=6)

			# create axis ranges and legend for `axes[2]`
			if t_max:
				axes[2].set_xlim(0, t_max)
			#axes[2].set_ylim(0, 30)
			axes[2].legend()
	else:
	
	#----------------------------------
	# Protein dynamics
		if X_cols.get("protein", 0) > 0:
			# set axis labels for `axes[2]`
			axes[2].set_xlabel(f"Time ({time_unit_str})")
			if X_name != "":
				axes[2].set_ylabel(f"{X_name} of protein concentration (µM)")
			else:
				axes[2].set_ylabel("Protein concentration (µM)")
			
			# plot data for `axes[2]`
			axes[2].plot(data_stacked_1[:,0], data_stacked_1[:,X_cols["protein"]], color="#008000", label=framework_name_1, marker='None', zorder=8)
			if n_batches > 1:
				axes[2].fill_between(data_stacked_1[:,0], data_stacked_1[:,X_cols["protein"]]-data_stacked_1[:,X_cols["protein"]+1]/np.sqrt(n_batches), 
														data_stacked_1[:,X_cols["protein"]]+data_stacked_1[:,X_cols["protein"]+1]/np.sqrt(n_batches), color="#008000", alpha=0.5)
			if data_stacked_2 is not None:
				axes[2].plot(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["protein"]], color="#003300", linestyle='dashed', label=framework_name_2, marker='None', zorder=8)
				if n_batches > 1:
					axes[2].fill_between(data_stacked_2[:,0], data_stacked_2[:,X_cols_2["protein"]]-data_stacked_2[:,X_cols_2["protein"]+1]/np.sqrt(n_batches), 
															data_stacked_2[:,X_cols_2["protein"]]+data_stacked_2[:,X_cols_2["protein"]+1]/np.sqrt(n_batches), color="#003300", alpha=0.5)

			# create axis ranges and legend for `axes[2]`
			if t_max:
				axes[2].set_xlim(0, t_max)
			axes[2].set_ylim(-0.05, 1.05)
			axes[2].legend()

	# save figure, e.g., as vector graphics
	fig.savefig(store_path, dpi=200)


#####################################
# plotFullComparison
# Plots comparison of the mean and s.e.m. (across 'n_batches') of mean and standard deviation of trials for synaptic weight, membrane voltage, and calcium/protein concentration,
# across simulators (Arbor, Brian, Stand-alone)
# - protocol: name of the protocol to be used (according '*.json' configuration, file must exist)
# - n_batches: number of batches
# - trials_per_batch: number of trials per batch
# - data_cols [optional]: specifies the columns where the wanted data is located
# - store_dir [optional]: folder in which to store resulting graphics files
# - file_format [optional]: format to use for storing graphics files
# - plot_sampling [optional]: sampling (1: best resolution, >1: samples are left out) -- can be useful to reduce file sizes
def plotFullComparison(protocol, n_batches, trials_per_batch,
                       data_cols = {"voltage": 1, "current": -1, "weight": 3, "calcium": 7, "protein": 9},
                       store_dir = ".", file_format = "svg", plot_sampling = 1):

	# load configuration
	if n_batches > 0:
		data_dir_standalone = os.path.abspath(f"stand-alone_data_{protocol}_{n_batches}x{trials_per_batch}")
		data_dir_arbor = os.path.abspath(f"arbor_data_{protocol}_{n_batches}x{trials_per_batch}")
		data_dir_brian = os.path.abspath(f"brian-heun_data_{protocol}_{n_batches}x{trials_per_batch}")
	else:
		data_dir_standalone = os.path.abspath(f"stand-alone_data_{protocol}")
		data_dir_arbor = os.path.abspath(f"arbor_data_{protocol}")
		data_dir_brian = os.path.abspath(f"brian-heun_data_{protocol}")
	config = json.load(open(os.path.join(data_dir_standalone, f"config_{protocol}.json"), "r")) # configuration must be the same for all simulators

	# load data
	mean_stacked_standalone = np.loadtxt(os.path.join(data_dir_standalone, f"meta_mean_averaged_traces.txt"))
	mean_stacked_arbor = np.loadtxt(os.path.join(data_dir_arbor, f"meta_mean_averaged_traces.txt"))
	mean_stacked_brian = np.loadtxt(os.path.join(data_dir_brian, f"meta_mean_averaged_traces.txt"))
	if trials_per_batch > 1:
		stdev_stacked_standalone = np.loadtxt(os.path.join(data_dir_standalone, f"meta_stdev_averaged_traces.txt"))
		stdev_stacked_arbor = np.loadtxt(os.path.join(data_dir_arbor, f"meta_stdev_averaged_traces.txt"))
		stdev_stacked_brian = np.loadtxt(os.path.join(data_dir_brian, f"meta_stdev_averaged_traces.txt"))

	# create storing directory
	os.makedirs(store_dir, exist_ok=True)

	# consider mean across trials (within batches)
	plotMeanAndSEMofX("Mean", n_batches, config, "Stand-alone", mean_stacked_standalone, None, "Brian", mean_stacked_brian, None,
	                  X_cols = data_cols, X_cols_2 = data_cols, 
					  store_path = os.path.join(store_dir, f"compare_mean_stand-alone_brian-heun_{protocol}.{file_format}"),
		              sampling = plot_sampling, time_unit_1 = 1000, time_unit_2 = 1)
	plotMeanAndSEMofX("Mean", n_batches, config, "Arbor", mean_stacked_arbor, None, "Brian", mean_stacked_brian, None,
	                  X_cols = data_cols, X_cols_2 = data_cols, 
					  store_path = os.path.join(store_dir, f"compare_mean_arbor_brian-heun_{protocol}.{file_format}"),
		              sampling = plot_sampling, time_unit_1 = 1, time_unit_2 = 1)
	plotMeanAndSEMofX("Mean", n_batches, config, "Arbor", mean_stacked_arbor, None, "Stand-alone", mean_stacked_standalone, None,
	                  X_cols = data_cols, X_cols_2 = data_cols, 
					  store_path = os.path.join(store_dir, f"compare_mean_arbor_stand-alone_{protocol}.{file_format}"),
		              sampling = plot_sampling, time_unit_1 = 1, time_unit_2 = 1000)

	# consider standard deviation across trials (within batches)
	if trials_per_batch > 1:
		plotMeanAndSEMofX("St.dev.", n_batches, config, "Stand-alone", stdev_stacked_standalone, None, "Brian", stdev_stacked_brian, None,
			              X_cols = data_cols, X_cols_2 = data_cols, 
						  store_path = os.path.join(store_dir, f"compare_stdev_stand-alone_brian-heun_{protocol}.{file_format}"),
		                  sampling = plot_sampling, wls = 0, time_unit_1 = 1000, time_unit_2 = 1)
		plotMeanAndSEMofX("St.dev.", n_batches, config, "Arbor", stdev_stacked_arbor, None, "Brian", stdev_stacked_brian, None,
			              X_cols = data_cols, X_cols_2 = data_cols, 
						  store_path = os.path.join(store_dir, f"compare_stdev_arbor_brian-heun_{protocol}.{file_format}"),
		                  sampling = plot_sampling, wls = 0, time_unit_1 = 1, time_unit_2 = 1)
		plotMeanAndSEMofX("St.dev.", n_batches, config, "Arbor", stdev_stacked_arbor, None, "Stand-alone", stdev_stacked_standalone, None,
			              X_cols = data_cols, X_cols_2 = data_cols, 
						  store_path = os.path.join(store_dir, f"compare_stdev_arbor_stand-alone_{protocol}.{file_format}"),
		                  sampling = plot_sampling, wls = 0, time_unit_1 = 1, time_unit_2 = 1000)


#####################################
# plotSingle
# Plots the mean and s.e.m. (across 'n_batches') of mean and standard deviation of trials for synaptic weight, membrane voltage, and calcium/protein concentration
# for a single simulator
# - simulator: name of the simulator to be used
# - protocol: name of the protocol to be used (according '*.json' configuration, file must exist)
# - n_batches: number of batches
# - trials_per_batch: number of trials per batch
# - data_cols [optional]: specifies the columns where the wanted data is located
# - store_dir [optional]: folder in which to store resulting graphics files
# - file_format [optional]: format to use for storing graphics files
# - plot_sampling [optional]: sampling (1: best resolution, >1: samples are left out) -- can be useful to reduce file sizes
# - time_unit [optional]: unit of time (1: millisecond, 1000: second, ...)
def plotSingle(simulator, protocol, n_batches, trials_per_batch,
               data_cols = {"voltage": 1, "current": -1, "weight": 3, "weight_l": 5, "calcium": 7, "protein": 9},
               store_dir = ".", file_format = "svg", plot_sampling = 1, time_unit = 1):

	# load configuration
	if n_batches > 0:
		data_dir = os.path.abspath(f"{simulator}_data_{protocol}_{n_batches}x{trials_per_batch}")
	else:
		data_dir = os.path.abspath(f"{simulator}_data_{protocol}")
	config = json.load(open(os.path.join(data_dir, f"config_{protocol}.json"), "r"))

	# load data
	mean_stacked = np.loadtxt(os.path.join(data_dir, f"meta_mean_averaged_traces.txt"))
	if trials_per_batch > 1:
		stdev_stacked = np.loadtxt(os.path.join(data_dir, f"meta_stdev_averaged_traces.txt"))

	# create storing directory
	os.makedirs(store_dir, exist_ok=True)

	# consider mean across trials (within batches)
	plotMeanAndSEMofX("Mean", n_batches, config, simulator, mean_stacked, None,
	                  X_cols = data_cols, store_path = os.path.join(store_dir, f"mean_{simulator}_{protocol}.{file_format}"),
		              sampling = plot_sampling, time_unit_1 = time_unit)

	# consider standard deviation across trials (within batches)
	if trials_per_batch > 1:
		plotMeanAndSEMofX("St.dev.", n_batches, config, simulator, stdev_stacked, None,
			              X_cols = data_cols, store_path = os.path.join(store_dir, f"stdev_{simulator}_{protocol}.{file_format}"),
		                  sampling = plot_sampling, wls = 0, time_unit_1 = time_unit)


#####################################
# plotDetComparison
# Plots a comparison of synaptic weight, membrane voltage, and calcium/protein concentration for deterministic data
# from two simulators
# - simulator_1: name of the first simulator
# - simulator_2: name of the second simulator
# - protocol: name of the protocol to be used (according data folder with 'config.json' file must exist)
# - label [optional]: additional label for output file
# - data_cols_1 [optional]: specifies the columns where the wanted data is located
# - data_cols_2 [optional]: specifies the columns where the wanted data is located
# - spike_neuron [optional]: index of the neuron the spikes of which shall be plotted
# - t_max [optional]: maximum time for plotting the data points (in ms, `None` means all data will be plotted)
# - store_dir [optional]: folder in which to store resulting graphics files
# - file_format [optional]: format to use for storing graphics files
# - plot_sampling [optional]: sampling (1: best resolution, >1: samples are left out) -- can be useful to reduce file sizes
def plotDetComparison(simulator_1, simulator_2, protocol, label = "",
               data_cols_1 = {"voltage": 1, "current": -1, "weight": 3, "calcium": 5, "protein": 6},
               data_cols_2 = {"voltage": 1, "current": -1, "weight": 3, "calcium": 5, "protein": 6},
               spike_neuron = None, t_max = None, store_dir = ".", file_format = "svg",
               plot_sampling = 1):

	# load configuration (from first simulator data directory)
	data_dirs = {simulator_1 : os.path.abspath(f"{simulator_1.lower()}_data_{protocol}"),
	             simulator_2 : os.path.abspath(f"{simulator_2.lower()}_data_{protocol}")}
	config = json.load(open(os.path.join(data_dirs[simulator_1], f"config.json"), "r"))

	# find and load data
	data_stacked = {}
	spikes_stacked = {}
	time_unit = {}
	for simulator in [simulator_1, simulator_2]:
		rawpaths = Path(data_dirs[simulator])
		timestamp = None
		for path in sorted(rawpaths.iterdir()):
			# get tail of file path
			tpath = os.path.split(str(path))[1] 
			# look for files whose name ends with '_traces.txt' or '_data.txt'
			if not path.is_dir() and hasTimestamp(path):
				for data_file_suffix in ["_traces.txt", "_data.txt"]:
					if data_file_suffix in tpath:
						timestamp = tpath.split(data_file_suffix)[0]
		if timestamp is None:
			raise FileNotFoundError(f"No data file found in '{data_dirs[simulator]}'.")
		if simulator.lower() in ["arbor", "brian"]:
			data_stacked[simulator] = np.loadtxt(os.path.join(data_dirs[simulator], timestamp + "_traces.txt"))
			spikes_stacked[simulator] = np.loadtxt(os.path.join(data_dirs[simulator], timestamp + "_spikes.txt")).transpose()
			# leave time unit as it is (ms)
			time_unit[simulator] = 1
		elif simulator.lower() == "stand-alone":
			data_stacked[simulator] = np.loadtxt(os.path.join(data_dirs[simulator], timestamp + "_data.txt"))
			spikes_stacked[simulator] = np.loadtxt(os.path.join(data_dirs[simulator], timestamp + "_spike_raster.txt")).transpose()
			# convert time unit to ms
			time_unit[simulator] = 1000
		else:
			raise ValueError(f"Unsupported simulator: '{simulator}'.")

	# count spikes for the specified neuron
	if spike_neuron is not None:
		if t_max:
			t_max_select = t_max
		else:
			t_max_select = data_stacked[simulator_1][:,0][-1]
		print(f"Spikes of neuron {spike_neuron} in {t_max_select} ms:")
		print(f"  {simulator_1}:", 
		      len(spikes_stacked[simulator_1][0][np.logical_and(spikes_stacked[simulator_1][0]*time_unit[simulator_1] < t_max_select, 
		                                                        spikes_stacked[simulator_1][1] == spike_neuron)]))
		print(f"  {simulator_2}:",
		      len(spikes_stacked[simulator_2][0][np.logical_and(spikes_stacked[simulator_2][0]*time_unit[simulator_2] < t_max_select, 
		                                                        spikes_stacked[simulator_2][1] == spike_neuron)]))

	# create storing directory
	os.makedirs(store_dir, exist_ok=True)

	# consider deterministic ("mean") values
	plotMeanAndSEMofX("", 1, config,
				      simulator_1, data_stacked[simulator_1], spikes_stacked[simulator_1], 
					  simulator_2, data_stacked[simulator_2], spikes_stacked[simulator_2], 
	                  X_cols = data_cols_1, X_cols_2 = data_cols_2, spike_neuron = spike_neuron, t_max = t_max,
					  store_path = os.path.join(store_dir, f"compare_det_{simulator_1.lower()}_{simulator_2.lower()}_{protocol}{label}.{file_format}"),
		              sampling = plot_sampling, time_unit_1 = time_unit[simulator_1], time_unit_2 = time_unit[simulator_2])
	

#####################################
# plotFullDetComparisonSmallnet3
# Plots a comparison of synaptic weight, membrane voltage, and calcium/protein concentration in all neurons and synapses of the 'smallnet3' setting,
# for the simulators Arbor and Stand-Alone
# - protocol [optional]: name of the protocol to be used (according data folder with 'config.json' file must exist)
# - t_max [optional]: maximum time for plotting the data points (in ms, `None` means all data will be plotted)
# - file_format [optional]: format to use for storing graphics files
def plotFullDetComparisonSmallnet3(protocol = "smallnet3_det", t_max = None, file_format = "svg"):
	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_0",
	                  {"voltage": 1, "current": -1, "weight": -1, "calcium": -1, "protein": 7}, 
					  {"voltage": 1, "current": -1, "weight": -1, "calcium": -1, "protein": 3},
					  spike_neuron = 0, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 0, no synapse

	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_1_synapse_0to1",
	                  {"voltage": 8, "current": -1, "weight": 10, "calcium": 12, "protein": 14}, 
					  {"voltage": 4, "current": -1, "weight": 16, "calcium": 18, "protein": 6},
					  spike_neuron = 1, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 1, synapse 0->1

	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_2_synapse_0to2",
	                  {"voltage": 15, "current": -1, "weight": 17, "calcium": 19, "protein": 21}, 
					  {"voltage": 7, "current": -1, "weight": 19, "calcium": 21, "protein": 9},
					  spike_neuron = 2, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 2, synapse 0->2
	
	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_2_synapse_3to2",
	                  {"voltage": 22, "current": -1, "weight": 24, "calcium": 26, "protein": 28}, 
					  {"voltage": 7, "current": -1, "weight": 25, "calcium": 27, "protein": 9},
					  spike_neuron = 2, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 2, synapse 3->2

	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_3_synapse_0to3",
	                  {"voltage": 36, "current": -1, "weight": 38, "calcium": 40, "protein": 42}, 
					  {"voltage": 10, "current": -1, "weight": 22, "calcium": 24, "protein": 12},
					  spike_neuron = 3, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 3, synapse 0->3
	
	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_3_synapse_2to3",
					  {"voltage": 29, "current": -1, "weight": 31, "calcium": 33, "protein": 35}, 
					  {"voltage": 10, "current": -1, "weight": 25, "calcium": 27, "protein": 12},
					  spike_neuron = 3, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 3, synapse 2->3
	
	plotDetComparison("Arbor", "Stand-alone", protocol, "_neuron_4",
	                  {"voltage": 43, "current": -1, "weight": -1, "calcium": -1, "protein": -1}, 
					  {"voltage": 13, "current": -1, "weight": -1, "calcium": -1, "protein": -1},
					  spike_neuron = 4, t_max = t_max,
					  store_dir = f"./t_max_{t_max}", file_format = file_format) # neuron 4 (inh.), no synapse

#####################################
if __name__ == '__main__':
	#plotSingle("stand-alone", "basic_early", 10, 100, file_format = "svg", time_unit=1000)
	plotSingle("stand-alone", "basic_late", 10, 10, file_format = "png", time_unit=1000)
	plotFullComparison("basic_early", 10, 100, file_format = "svg")
	plotFullComparison("basic_late", 10, 10, file_format = "png", \
					    data_cols = {"voltage": 1, "current": -1, "weight": 3, "weight_l": 5, "calcium": 7, "protein": 9})
	#plotFullComparison("basic_late", 10, 10, file_format = "svg", \
	#                   data_cols = {"voltage": -1, "current": -1, "weight": 3, "weight_l": 5, "calcium": 7, "protein": 9}, \
	#                   plot_sampling = 20)
	#plotFullDetComparisonSmallnet3("smallnet3_det_8h-recall", t_max = 20000)
	#plotFullDetComparisonSmallnet3("smallnet3_det_8h-recall")
