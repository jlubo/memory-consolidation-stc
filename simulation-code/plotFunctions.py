##############################################################################################
### Code to create plots of the neurons and synaptic weights in the network via Matplotlib ###
##############################################################################################

### Copyright 2017-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0), except the function "shiftedColorMap"

### Uses the function shiftedColorMap -
### Copyright 2014 stackoverflow.com user Paul H
### licensed under CC-BY-SA 4.0 (https://creativecommons.org/licenses/by-sa/4.0/)

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib import cm
import os.path

plot_folder = 'network_plots/' # folder for plot output
log_shift = 1.03 # shift before taking logarithm
epsilon = 1e-11 # very small number that is counted as zero
rel_epsilon = 1e-4 # fraction of maximum plasticity value that has to be exceeded to acknowledge occurrence of plasticity
np.set_printoptions(precision=8, threshold=1e10, linewidth=200) # precision for number output

# shiftedColorMap
# source: https://stackoverflow.com/questions/7404116/defining-the-midpoint-of-a-colormap-in-matplotlib/20528097#20528097
def shiftedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, name='shiftedcmap'):
    '''
    Function to offset the "center" of a colormap. Useful for
    data with a negative min and positive max and you want the
    middle of the colormap's dynamic range to be at zero.
    
    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower offset). Should be between
          0.0 and `midpoint`.
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0. In
          general, this should be  1 - vmax / (vmax + abs(vmin))
          For example if your data range from -15.0 to +5.0 and
          you want the center of the colormap at 0.0, `midpoint`
          should be set to  1 - 5/(5 + 15)) or 0.75
      stop : Offset from highest point in the colormap's range.
          Defaults to 1.0 (no upper offset). Should be between
          `midpoint` and 1.0.
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }
      
    # regular index to compute the colors
    reg_index = np.linspace(start, stop, 257)

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129, endpoint=True)
    ])
    
    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))
        
    newcmap = cm.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap


# logScale functions
# Enables logarithmic presentation of data containing zero values
# x: the data
# h_0: initial early-phase/total weight
# z_max: maximum late-phase weight
# w_max: maximum total weight
def logScale_h(x, h_0):
	#return x # to disable logarithmic scale
	xprime = x - h_0
	h_max = 2*h_0

	neg_vals = (xprime < -(h_max*rel_epsilon))
	pos_vals = (xprime > (h_max*rel_epsilon))
	noconn = (np.abs(x) < (h_max*rel_epsilon)) # no connection exists

	comp1 = neg_vals * (-h_0 * np.log(-xprime+log_shift)/np.log(log_shift+h_0) + h_0)
	comp2 = pos_vals * (h_0 * np.log(xprime+log_shift)/np.log(log_shift+h_0) + h_0)
	comp3 = np.logical_or(np.logical_not(np.logical_or(neg_vals, pos_vals)), noconn) * h_0 # neither LTP nor LTD, or no connection exists -> at h_0

	return comp1 + comp2 + comp3 # log_shift+h_0 because h_0 is the maximum value the shifted xprime can reach

def logScale_z(x, z_max):
	return x # to disable logarithmic scale

	comp1 = (x < -(z_max*rel_epsilon)) * -np.log(-x+log_shift)/np.log(log_shift+z_max) # log_shift+1 because 1 is the maximum value x could reach
	comp2 = (x > (z_max*rel_epsilon)) * np.log(x+log_shift)/np.log(log_shift+z_max)

	return comp1 + comp2

def logScale_w(x, h_0, z_max):
	xprime = x - h_0
	w_max = (2 + z_max)*h_0

	#return x + (np.abs(x) < (w_max*rel_epsilon/100)) * h_0 # to disable logarithmic scale

	neg_vals = (xprime < -(w_max*rel_epsilon)) # significant LTD
	pos_vals = (xprime > (w_max*rel_epsilon)) # significant LTP
	noconn = (np.abs(x) < (w_max*rel_epsilon)) # no connection exists

	comp1 = neg_vals * (-h_0 * np.log(-xprime*neg_vals+log_shift)/np.log(log_shift+h_0) + h_0)
	comp2 = pos_vals * (h_0 * np.log(xprime*pos_vals+log_shift)/np.log(log_shift+h_0) + h_0)
	comp3 = np.logical_or(np.logical_not(np.logical_or(neg_vals, pos_vals)), noconn) * h_0 # neither LTP nor LTD, or no connection exists -> at h_0

	return comp1 + comp2 + comp3 # log_shift+h_0 because h_0 is the maximum value the shifted xprime can reach

# readWeightMatrixData
# Reads complete weight matrix data from a file
# filename: name of the file to read the data from
# Nl_exc: number of excitatory neurons in one row/column
# return: the adjacency matrix, the early-phase weight matrix, the late-phase weight matrix, the firing rate vector
def readWeightMatrixData(filename, Nl_exc):

	# read weight matrices and firing rates from file
	with open(plot_folder + filename) as f:
		rawdata = f.read()

	rawdata = rawdata.split('\n\n')
	rawmatrix_h = rawdata[0].split('\n')
	rawmatrix_z = rawdata[1].split('\n')
	rawmatrix_v = rawdata[2].split('\n')

	rows = len(rawmatrix_v)

	if (rows != len(rawmatrix_v[0].split('\t\t'))) or (rows != Nl_exc):
		print('Data file error in "' + filename + '"')
		f.close()
		exit()

	v = np.zeros((Nl_exc,Nl_exc))
	h = np.zeros((Nl_exc**2,Nl_exc**2))
	z = np.zeros((Nl_exc**2,Nl_exc**2))

	for i in range(Nl_exc**2):
		if i < Nl_exc:
			value0 = rawmatrix_v[i].split('\t\t')
		value1 = rawmatrix_h[i].split('\t\t')
		value2 = rawmatrix_z[i].split('\t\t')

		for j in range(Nl_exc**2):
			if i < Nl_exc and j < Nl_exc:
				v[i][j] = float(value0[j])
			h[i][j] = float(value1[j])
			z[i][j] = float(value2[j])

	f.close()
	connections = (h > epsilon)

	return connections, h, z, v

# plotWeights
# Plots synaptic weights from a complete weight matrix
# filename: name of the file to read the data from
# h_0: initial early-phase weight
# z_min: minimum late-phase weight
# z_max: maximum late-phase weight
# Nl_exc: number of excitatory neurons in one row/column
# title [optional]: main title of the figure
def plotWeights(filename, h_0, z_min, z_max, Nl_exc, title = "Weight matrices"):

	# normalization factors
	h_max = 2*h_0
	w_max = h_max + h_0 * z_max

	# colormaps for h and z
	cmh0center = shiftedColorMap(cm.seismic, start=0, midpoint=0.5, stop=1.0, name='h0center')
	cmasymmetric = shiftedColorMap(cm.seismic, start=0, midpoint=abs(z_min), stop=z_max+abs(z_min), name='asymmetric')
	cmasymmetric_w = shiftedColorMap(cm.seismic, start=0, midpoint=h_0/w_max, stop=1, name='asymmetric_w')

	# read weight matrix data
	connections, h, z, v = readWeightMatrixData(filename, Nl_exc)
	connections = np.flip(connections, 0)
	h = np.flip(h, 0)
	z = np.flip(z, 0)

	# normalize firing rates
	v_max = np.amax(v)
	if v_max == 0.0:
        	v_max = 1.0
	v_data = v.reshape(int(Nl_exc),int(Nl_exc)) / v_max

	# plotting
	plt.figure(figsize=(20,12))
	grid = gs.GridSpec(2, 2)
	plt.rc('font', size=13)

	ax1 = plt.subplot(grid[0, 0])
	p1 = ax1.imshow(v_data, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Oranges, origin='lower')
	ax1.set_title('Firing rate v')
	ax1.xaxis.set_ticks([])
	ax1.yaxis.set_ticks([])
	plt.colorbar(p1, ticks=[0, 0.5, 1])

	ax2 = plt.subplot(grid[0, 1])
	p2 = ax2.imshow(logScale_w((h + h_0*z), h_0, z_max) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmasymmetric_w, origin='lower') # non-existing connections are also shown in white
	ax2.set_title('Total weight w')
	ax2.xaxis.set_ticks([])
	ax2.yaxis.set_ticks([])
	plt.colorbar(p2, ticks=[0, 1])

	ax3 = plt.subplot(grid[1, 0])
	p3 = ax3.imshow(logScale_h(h, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower') # non-existing connections are also shown in white
	plt.colorbar(p3, ticks=[0, 1])
	ax3.set_title('Early-phase weight h')
	ax3.xaxis.set_ticks([])
	ax3.yaxis.set_ticks([])

	ax4 = plt.subplot(grid[1, 1])
	p4 = ax4.imshow(z, vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
	plt.colorbar(p4, ticks=[z_min, 0, z_max])
	ax4.set_title('Late-phase weight z')
	ax4.xaxis.set_ticks([])
	ax4.yaxis.set_ticks([])

	# set general picture frame
	plt.tight_layout()
	plt.suptitle(title)  # main title
	plt.subplots_adjust(top=0.9)  # vertical space for main title
	plt.subplots_adjust(right=0.9)  # horizontal space for main title
	plt.draw()  # draw plot but do not show
	plt.savefig(plot_folder + filename.replace('.txt', '.png'), bbox_inches='tight', dpi=300)
	plt.close()


# plotWeightDiffs
# Plots the difference between the synaptic weights from a weight matrix and those of another weight matrix
# filename1: name of the file to read the first matrix from
# filename2: name of the file to read the second matrix from
# h_0: initial early-phase weight
# z_min: minimum late-phase weight
# z_max: maximum late-phase weight
# Nl_exc: number of excitatory neurons in one row/column
# title [optional]: main title of the figure
def plotWeightDiffs(filename1, filename2, h_0, z_min, z_max, Nl_exc, title = "Weight matrices"):

	# Colormaps for Calcium simulation
	# colormaps for h and z
	cmh0center = shiftedColorMap(cm.seismic, start=0, midpoint=0.5, stop=1.0, name='h0center')
	cmasymmetric = shiftedColorMap(cm.seismic, start=0, midpoint=abs(z_min), stop=z_max+abs(z_min), name='asymmetric')

	# normalization factors
	h_max = 2*h_0
	w_max = h_max + h_0 * z_max

	# read weight matrix data
	connections1, h1, z1, v1 = readWeightMatrixData(filename1, Nl_exc)
	connections2, h2, z2, v2 = readWeightMatrixData(filename2, Nl_exc)

	if (connections1 != connections2).any():
		print("Not the same connectivity, plot cannot be created!")
		return

	connections = np.flip(connections1, 0)
	h = np.flip(connections1*(h2-h1), 0)
	z = np.flip(connections1*(z2-z1), 0)
	v = v2-v1

	# normalize firing rates
	v_max = np.amax(v)
	if v_max == 0.0:
        	v_max = 1.0
	v_data = v.reshape(int(Nl_exc),int(Nl_exc)) / v_max

	# plotting
	plt.figure(figsize=(20,12))
	grid = gs.GridSpec(2, 2)
	plt.rc('font', size=13)

	ax1 = plt.subplot(grid[0, 0])
	p1 = ax1.imshow(v_data, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Oranges, origin='lower')
	ax1.set_title('Firing rate v')
	ax1.xaxis.set_ticks([])
	ax1.yaxis.set_ticks([])
	plt.colorbar(p1, ticks=[0, 0.5, 1])

	ax2 = plt.subplot(grid[0, 1])
	p2 = ax2.imshow((h + h_0*z) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower') # non-existing connections are also shown in white
	ax2.set_title('Total weight w')
	ax2.xaxis.set_ticks([])
	ax2.yaxis.set_ticks([])
	plt.colorbar(p2, ticks=[0, 0.5, 1])

	ax3 = plt.subplot(grid[1, 0])
	p3 = ax3.imshow(logScale_h(h + h_0, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower') # non-existing connections are also shown in white
	plt.colorbar(p3, ticks=[0, 1])
	ax3.set_title('Early-phase weight h')
	ax3.xaxis.set_ticks([])
	ax3.yaxis.set_ticks([])

	ax4 = plt.subplot(grid[1, 1])
	p4 = ax4.imshow(logScale_z(z, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
	plt.colorbar(p4, ticks=[z_min, 0, z_max])
	ax4.set_title('Late-phase weight z')
	ax4.xaxis.set_ticks([])
	ax4.yaxis.set_ticks([])

	# set general picture frame
	plt.tight_layout()
	plt.suptitle(title)  # main title
	plt.subplots_adjust(top=0.9)  # vertical space for main title
	plt.subplots_adjust(right=0.9)  # horizontal space for main title
	plt.draw()  # draw plot but do not show
	plt.savefig(plot_folder + filename2.replace('.txt', '_diff.png'), bbox_inches='tight', dpi=300)
	plt.close()

# plotAveragedWeights
# Plots averaged synaptic weights, either from a complete weight matrix or from already averaged data that is saved neuronwise
# filename: name of the file to read the data from
# h_0: initial early-phase weight
# z_min: minimum late-phase weight
# z_max: maximum late-phase weight
# Nl_exc: number of excitatory neurons in one row/column
# already_averaged: specifies if a data file shall be created
# title [optional]: main title of the figure
def plotAveragedWeights(filename, h_0, z_min, z_max, Nl_exc, already_averaged, title = "Averaged incoming and outgoing weights"):

	# colormaps for h and z
	cmh0center = shiftedColorMap(cm.seismic, start=0, midpoint=0.5, stop=1.0, name='h0center')
	cmasymmetric = shiftedColorMap(cm.seismic, start=0, midpoint=abs(z_min), stop=z_max+abs(z_min), name='asymmetric')

	# normalization factors
	h_max = 2*h_0
	w_max = h_max + h_0 * z_max

	# using alread averaged data
	if already_averaged:

		# read from file and normalize
		with open(plot_folder + filename) as f:
			rawdata = f.read()

		rawdata = rawdata.split('\n')
		nn = len(rawdata)-1

		if nn != 2*Nl_exc*Nl_exc:
			print('Data file error in "' + filename + '"')
			print(nn)
			f.close()
			exit()

		v = np.zeros((Nl_exc,Nl_exc))
		h_inc = np.zeros((Nl_exc,Nl_exc))
		h_out = np.zeros((Nl_exc,Nl_exc))
		z_inc = np.zeros((Nl_exc,Nl_exc))
		z_out = np.zeros((Nl_exc,Nl_exc))
		w_inc = np.zeros((Nl_exc,Nl_exc))
		w_out = np.zeros((Nl_exc,Nl_exc))

		for n in range(nn):
			n2 = n % (Nl_exc*Nl_exc)
			i = (n2 - (n2 % Nl_exc)) // Nl_exc # row number
			j = n2 % Nl_exc # column number

			if n < Nl_exc*Nl_exc:
				values = rawdata[n].split()
				h_inc[i][j] = logScale_h(float(values[0]), h_0) / h_max
				h_out[i][j] = logScale_h(float(values[1]), h_0) / h_max
				z_inc[i][j] = logScale_z(float(values[2]), z_max)
				z_out[i][j] = logScale_z(float(values[3]), z_max)
				w_inc[i][j] = float(float(values[0]) + float(values[2])) / w_max
				w_out[i][j] = float(float(values[1]) + float(values[3])) / w_max

			else:
				v[i][j] = float(rawdata[n])

		f.close()

		# find firing rate maximum
		v_max = np.amax(v)

		# do the plotting
		plt.figure(figsize=(6,12))
		grid = gs.GridSpec(4, 2)
		plt.rc('font', size=13)

		ax1 = plt.subplot(grid[0, :])
		p1 = ax1.imshow(v, vmin=0, vmax=v_max, interpolation='nearest', cmap=cm.Oranges, origin='lower')
		ax1.set_title('v')
		ax1.xaxis.set_ticks([])
		ax1.yaxis.set_ticks([])
		plt.colorbar(p1, ticks=[0, v_max/2, v_max])

		ax3 = plt.subplot(grid[1, 0])
		p3 = ax3.imshow(w_inc, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Blues, origin='lower')
		cb = plt.colorbar(p3, ticks=[0, 1])
		cb.set_ticklabels(["0", "3 h_0"])
		ax3.set_title('⟨ Incoming w ⟩')
		ax3.xaxis.set_ticks([])
		ax3.yaxis.set_ticks([])

		ax4 = plt.subplot(grid[1, 1])
		p4 = ax4.imshow(w_out, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
		cb = plt.colorbar(p4, ticks=[0, 1])
		cb.set_ticklabels(["0","3 h_0"])
		ax4.set_title('⟨ Outgoing w ⟩')
		ax4.xaxis.set_ticks([])
		ax4.yaxis.set_ticks([])

		ax5 = plt.subplot(grid[2, 0])
		p5 = ax5.imshow(h_inc, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
		cb = plt.colorbar(p5, ticks=[0, 0.5, 1])
		cb.set_ticklabels(["0", "h_0", "2 h_0"])
		ax5.set_title('⟨ Incoming h ⟩')
		ax5.xaxis.set_ticks([])
		ax5.yaxis.set_ticks([])

		ax6 = plt.subplot(grid[2, 1])
		p6 = ax6.imshow(h_out, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
		cb = plt.colorbar(p6, ticks=[0, 0.5, 1])
		cb.set_ticklabels(["0", "h_0", "2 h_0"])
		ax6.set_title('⟨ Outgoing h ⟩')
		ax6.xaxis.set_ticks([])
		ax6.yaxis.set_ticks([])

		ax7 = plt.subplot(grid[3, 0])
		p7 = ax7.imshow(z_inc, vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
		plt.colorbar(p7, ticks=[z_min, 0, z_max])
		ax7.set_title('⟨ Incoming z ⟩')
		ax7.xaxis.set_ticks([])
		ax7.yaxis.set_ticks([])

		ax8 = plt.subplot(grid[3, 1])
		p8 = ax8.imshow(z_out, vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
		plt.colorbar(p8, ticks=[z_min, 0, z_max])
		ax8.set_title('⟨ Outgoing z ⟩')
		ax8.xaxis.set_ticks([])
		ax8.yaxis.set_ticks([])

		# set general picture frame
		plt.tight_layout()
		plt.suptitle(title)  # main title
		plt.subplots_adjust(top=0.9)  # vertical space for main title
		plt.subplots_adjust(right=0.9)  # horizontal space for main title
		plt.text(-1.0, 0.5, 'Firing rate', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax1.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Total weight', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax3.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Early phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax5.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Late phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax7.transAxes) # add label for each measure
		plt.draw()  # draw plot but do not show

		plt.savefig(plot_folder + filename.replace('.txt', '.png'), bbox_inches='tight', dpi=300)
		plt.close()

	# using a complete weight matrix and averaging over it
	else:

		# read weight matrix data
		connections, h, z, v = readWeightMatrixData(filename, Nl_exc)

		# change filename
		filename_av = filename.replace('_net_', '_net_av_')

		# find firing rate maximum and reshape array
		v_max = np.amax(v)
		v_data = v.reshape(int(Nl_exc),int(Nl_exc))

		# average incoming (axis=0) synaptic weights per neuron
		con_count = np.sum(connections, axis=0)
		con_count[con_count == 0] = 1
		h_incoming = np.array(np.sum(h, axis=0) / con_count).reshape(int(Nl_exc),int(Nl_exc))
		z_incoming = np.array(np.sum(z, axis=0) / con_count).reshape(int(Nl_exc),int(Nl_exc))

		# average outgoing (axis=1) synaptic weights per neuron
		con_count = np.sum(connections, axis=1)
		con_count[con_count == 0] = 1
		h_outgoing = np.array(np.sum(h, axis=1) / con_count).reshape(int(Nl_exc),int(Nl_exc))
		z_outgoing = np.array(np.sum(z, axis=1) / con_count).reshape(int(Nl_exc),int(Nl_exc))

		# plotting
		plt.figure(figsize=(6,12))
		grid = gs.GridSpec(4, 2)
		plt.rc('font', size=13)

		ax1 = plt.subplot(grid[0, :])
		p1 = ax1.imshow(v_data, vmin=0, vmax=v_max, interpolation='nearest', cmap=cm.Oranges, origin='lower')
		ax1.set_title('v')
		ax1.xaxis.set_ticks([])
		ax1.yaxis.set_ticks([])
		plt.colorbar(p1, ticks=[0, v_max/2, v_max])

		ax3 = plt.subplot(grid[1, 0])
		p3 = ax3.imshow((h_incoming + h_0*z_incoming) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Blues, origin='lower')
		cb = plt.colorbar(p3, ticks=[0, 1])
		cb.set_ticklabels(["0","3 h_0"])
		ax3.set_title('⟨ Incoming w ⟩')
		ax3.xaxis.set_ticks([])
		ax3.yaxis.set_ticks([])

		ax4 = plt.subplot(grid[1, 1])
		p4 = ax4.imshow((h_outgoing + h_0*z_outgoing) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
		cb = plt.colorbar(p4, ticks=[0, 1])
		cb.set_ticklabels(["0","3 h_0"])
		ax4.set_title('⟨ Outgoing w ⟩')
		ax4.xaxis.set_ticks([])
		ax4.yaxis.set_ticks([])

		ax5 = plt.subplot(grid[2, 0])
		p5 = ax5.imshow(logScale_h(h_incoming, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
		cb = plt.colorbar(p5, ticks=[0, 0.5, 1])
		cb.set_ticklabels(["0", "h_0", "2 h_0"])
		ax5.set_title('⟨ Incoming h ⟩')
		ax5.xaxis.set_ticks([])
		ax5.yaxis.set_ticks([])

		ax6 = plt.subplot(grid[2, 1])
		p6 = ax6.imshow(logScale_h(h_outgoing, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
		cb = plt.colorbar(p6, ticks=[0, 0.5, 1])
		cb.set_ticklabels(["0", "h_0", "2 h_0"])
		ax6.set_title('⟨ Outgoing h ⟩')
		ax6.xaxis.set_ticks([])
		ax6.yaxis.set_ticks([])

		ax7 = plt.subplot(grid[3, 0])
		p7 = ax7.imshow(logScale_z(z_incoming, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
		plt.colorbar(p7, ticks=[z_min, 0, z_max])
		ax7.set_title('⟨ Incoming z ⟩')
		ax7.xaxis.set_ticks([])
		ax7.yaxis.set_ticks([])

		ax8 = plt.subplot(grid[3, 1])
		p8 = ax8.imshow(logScale_z(z_outgoing, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
		plt.colorbar(p8, ticks=[z_min, 0, z_max])
		ax8.set_title('⟨ Outgoing z ⟩')
		ax8.xaxis.set_ticks([])
		ax8.yaxis.set_ticks([])

		# set general picture frame
		plt.tight_layout()
		plt.suptitle(title)  # main title
		plt.subplots_adjust(top=0.9)  # vertical space for main title
		plt.subplots_adjust(right=0.9)  # horizontal space for main title
		plt.text(-1.0, 0.5, 'Firing rate', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax1.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Total weight', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax3.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Early phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax5.transAxes) # add label for each measure
		plt.text(-0.3, 0.5, 'Late phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax7.transAxes) # add label for each measure
		plt.draw()  # draw plot but do not show
		plt.savefig(plot_folder + filename_av.replace('.txt', '.png'), bbox_inches='tight', dpi=300)
		plt.close()

		# output of averaged data to file
		f = open(plot_folder + filename_av, 'w')
		f.write("v_data =\r\n" + str(v_data) + "\r\n\r\n\r\n")
		f.write("h_incoming =\r\n" + str(h_incoming) + "\r\n\r\n\r\n")
		f.write("h_outgoing =\r\n" + str(h_outgoing) + "\r\n\r\n\r\n")
		f.write("z_incoming =\r\n" + str(z_incoming) + "\r\n\r\n\r\n")
		f.write("z_outgoing =\r\n" + str(z_outgoing) + "\r\n\r\n\r\n")
		f.close()

# earlyPhaseWeightsFromCore
# Returns all the early-phase synaptic weights incoming to neuron i from core neurons
# i: neuron index
# adj: the adjacency matrix
# h: the early-phase weights
# core: the neurons belonging to the core
def earlyPhaseWeightsFromCore(i, adj, h, core):

	hi = h[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(hi)), core))
	inc_h = hi[adj_connections]

	return inc_h

# earlyPhaseWeightsFromControl
# Returns all the early-phase synaptic weights incoming to neuron i from control neurons
# i: neuron index
# adj: the adjacency matrix
# h: the early-phase weights
# core: the neurons belonging to the core
def earlyPhaseWeightsFromControl(i, adj, h, core):

	hi = h[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.logical_not(np.in1d(np.arange(len(hi)), core)))
	inc_h = hi[adj_connections]

	return inc_h

# earlyPhaseWeightsToCore
# Returns all the early-phase synaptic weights outgoing to core neurons from neuron i
# i: neuron index
# adj: the adjacency matrix
# h: the early-phase weights
# core: the neurons belonging to the core
def earlyPhaseWeightsToCore(i, adj, h, core):

	hi = h[i, :]
	adj_connections = np.logical_and(adj[i, :] == 1, np.in1d(np.arange(len(hi)), core))
	out_h = hi[adj_connections]

	return out_h

# latePhaseWeightsFromCore
# Returns all the late-phase synaptic weights incoming to neuron i from core neurons
# i: neuron index
# adj: the adjacency matrix
# z: the late-phase weights
# core: the neurons belonging to the core
def latePhaseWeightsFromCore(i, adj,z, core):

	zi = z[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.in1d(np.arange(len(zi)), core))
	inc_z = zi[adj_connections]

	return inc_z

# latePhaseWeightsFromControl
# Returns all the late-phase synaptic weights incoming to neuron i from control neurons
# i: neuron index
# adj: the adjacency matrix
# z: the late-phase weights
# core: the neurons belonging to the core
def latePhaseWeightsFromControl(i, adj,z, core):

	zi = z[:, i]
	adj_connections = np.logical_and(adj[:, i] == 1, np.logical_not(np.in1d(np.arange(len(zi)), core)))
	inc_z = zi[adj_connections]

	return inc_z

# latePhaseWeightsToCore
# Prints and returns all the late-phase synaptic weights outgoing to core neurons from neuron i
# i: neuron index
# adj: the adjacency matrix
# z: the late-phase weights
# core: the neurons belonging to the core
def latePhaseWeightsToCore(i, adj, z, core):

	zi = z[i, :]
	adj_connections = np.logical_and(adj[i, :] == 1, np.in1d(np.arange(len(zi)), core))
	out_z = zi[adj_connections]

	return out_z

# plotAveragedSubPopWeights
# Plots the averaged weight incoming from core neurons and from control neurons, respectively, from a complete weight matrix
# filename: name of the file to read the data from
# h_0: initial early-phase weight
# z_min: minimum late-phase weight
# z_max: maximum late-phase weight
# Nl_exc: number of excitatory neurons in one row/column
# s_CA: size of the cell assembly core
# title [optional]: main title of the figure
def plotAveragedSubPopWeights(filename, h_0, z_min, z_max, Nl_exc, s_CA, title = "Averaged incoming weights from subpopulations"):

	# colormaps for h and z
	cmh0center = shiftedColorMap(cm.seismic, start=0, midpoint=0.5, stop=1.0, name='h0center')
	cmasymmetric = shiftedColorMap(cm.seismic, start=z_min, midpoint=0, stop=z_max, name='asymmetric')

	# normalization factors
	h_max = 2*h_0
	w_max = h_max + h_0 * z_max

	# read weight matrix data
	connections, h, z, v = readWeightMatrixData(filename, Nl_exc)
	print("=========================\nData from: " + filename)

	# change filename
	filename_av = filename.replace('_net_', '_net_subpop_av_')

	# normalize firing rates
	v_max = np.amax(v)
	if v_max == 0.0:
        	v_max = 1.0
	v_data = v.reshape(int(Nl_exc),int(Nl_exc)) / v_max

	# average synaptic weight from core neurons and from control neurons for all neurons
	core = np.arange(s_CA)

	h_core_inc = np.zeros((Nl_exc*Nl_exc))
	z_core_inc = np.zeros((Nl_exc*Nl_exc))
	#h_core_out = np.zeros((Nl_exc*Nl_exc))
	#z_core_out = np.zeros((Nl_exc*Nl_exc))
	h_control_inc = np.zeros((Nl_exc*Nl_exc))
	z_control_inc = np.zeros((Nl_exc*Nl_exc))

	for i in range(Nl_exc*Nl_exc):
		weights = earlyPhaseWeightsFromCore(i, connections, h, core)
		if weights != []:
			h_core_inc[i] = np.sum(weights) / len(weights)

		weights = latePhaseWeightsFromCore(i, connections, z, core)
		if weights != []:
			z_core_inc[i] = np.sum(weights) / len(weights)

		#weights = earlyPhaseWeightsToCore(i, connections, h, core)
		#h_core_out[i] = np.sum(weights) / len(weights)

		#weights = latePhaseWeightsToCore(i, connections, z, core)
		#z_core_out[i] = np.sum(weights) / len(weights)

		weights = earlyPhaseWeightsFromControl(i, connections, h, core)
		if weights != []:
			h_control_inc[i] = np.sum(weights) / len(weights)

		weights = latePhaseWeightsFromControl(i, connections, z, core)
		if weights != []:
			z_control_inc[i] = np.sum(weights) / len(weights)


	# compute mean synaptic weight from core to core, from core to control, from control to core and from control to control
	core_mask = np.in1d(np.arange(Nl_exc*Nl_exc), core) # boolean mask of exc. population, entries are True for core neurons
	control_mask = np.logical_not(core_mask) # boolean mask of exc. population, entries are True for control neurons
	core_size = len(core)
	control_size = Nl_exc*Nl_exc - core_size

	mean_h_core_core = np.sum(h_core_inc[core_mask]) / core_size
	mean_h_core_control = np.sum(h_core_inc[control_mask]) / control_size
	mean_h_control_core = np.sum(h_control_inc[core_mask]) / core_size
	mean_h_control_control = np.sum(h_control_inc[control_mask]) / control_size

	mean_z_core_core = np.sum(z_core_inc[core_mask]) / core_size
	mean_z_core_control = np.sum(z_core_inc[control_mask]) / control_size
	mean_z_control_core = np.sum(z_control_inc[core_mask]) / core_size
	mean_z_control_control = np.sum(z_control_inc[control_mask]) / control_size

	print("mean_h_core_core = " + str(mean_h_core_core))
	print("mean_h_core_control = " + str(mean_h_core_control))
	print("mean_h_control_core = " + str(mean_h_control_core))
	print("mean_h_control_control = " + str(mean_h_control_control))

	print("mean_z_core_core = " + str(mean_z_core_core))
	print("mean_z_core_control = " + str(mean_z_core_control))
	print("mean_z_control_core = " + str(mean_z_control_core))
	print("mean_z_control_control = " + str(mean_z_control_control))

	# reshape neuron arrays for plotting
	h_core_inc = h_core_inc.reshape(int(Nl_exc),int(Nl_exc))
	z_core_inc = z_core_inc.reshape(int(Nl_exc),int(Nl_exc))
	#h_core_out = h_core_out.reshape(int(Nl_exc),int(Nl_exc))
	#z_core_out = z_core_out.reshape(int(Nl_exc),int(Nl_exc))
	h_control_inc = h_control_inc.reshape(int(Nl_exc),int(Nl_exc))
	z_control_inc = z_control_inc.reshape(int(Nl_exc),int(Nl_exc))

	# plotting
	plt.figure(figsize=(6,12))
	grid = gs.GridSpec(4, 2)
	plt.rc('font', size=13)

	ax1 = plt.subplot(grid[0, :])
	p1 = ax1.imshow(v_data, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Oranges, origin='lower')
	ax1.set_title('v')
	ax1.xaxis.set_ticks([])
	ax1.yaxis.set_ticks([])
	plt.colorbar(p1, ticks=[0, 0.5, 1])

	ax3 = plt.subplot(grid[1, 0])
	p3 = ax3.imshow((h_core_inc + h_0*z_core_inc) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Blues, origin='lower')
	cb = plt.colorbar(p3, ticks=[0, 1])
	cb.set_ticklabels(["0","3 h_0"])
	ax3.set_title('⟨ Inc. w from core ⟩')
	ax3.xaxis.set_ticks([])
	ax3.yaxis.set_ticks([])

	ax4 = plt.subplot(grid[1, 1])
	#p4 = ax4.imshow((h_core_out + h_0*z_core_out) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
	#ax4.set_title('⟨ Outgoing w to core ⟩')
	p4 = ax4.imshow((h_control_inc + h_0*z_control_inc) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
	cb = plt.colorbar(p4, ticks=[0, 1])
	cb.set_ticklabels(["0","3 h_0"])
	ax4.set_title('⟨ Inc. w from control ⟩')
	ax4.xaxis.set_ticks([])
	ax4.yaxis.set_ticks([])

	ax5 = plt.subplot(grid[2, 0])
	p5 = ax5.imshow(logScale_h(h_core_inc, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
	ax5.set_title('⟨ Inc. h from core ⟩')
	cb = plt.colorbar(p5, ticks=[0, 0.5, 1])
	cb.set_ticklabels(["0", "h_0", "2 h_0"])
	ax5.xaxis.set_ticks([])
	ax5.yaxis.set_ticks([])

	ax6 = plt.subplot(grid[2, 1])
	#p6 = ax6.imshow(logScale_h(h_core_out, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
	#ax6.set_title('⟨ Outgoing h to core ⟩')
	p6 = ax6.imshow(logScale_h(h_control_inc, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmh0center, origin='lower')
	ax6.set_title('⟨ Inc. h from control ⟩')
	cb = plt.colorbar(p6, ticks=[0, 0.5, 1])
	cb.set_ticklabels(["0", "h_0", "2 h_0"])
	ax6.xaxis.set_ticks([])
	ax6.yaxis.set_ticks([])

	ax7 = plt.subplot(grid[3, 0])
	p7 = ax7.imshow(logScale_z(z_core_inc, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
	ax7.set_title('⟨ Inc. z from core ⟩')
	plt.colorbar(p7, ticks=[z_min, 0, z_max])
	ax7.xaxis.set_ticks([])
	ax7.yaxis.set_ticks([])

	ax8 = plt.subplot(grid[3, 1])
	#p8 = ax8.imshow(logScale_z(z_core_out, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
	#ax8.set_title('⟨ Outgoing z to core⟩')
	p8 = ax8.imshow(logScale_z(z_control_inc, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmasymmetric, origin='lower')
	ax8.set_title('⟨ Inc. z from control ⟩')
	plt.colorbar(p8, ticks=[z_min, 0, z_max])
	ax8.xaxis.set_ticks([])
	ax8.yaxis.set_ticks([])

	# set general picture frame
	plt.tight_layout()
	plt.suptitle(title)  # main title
	plt.subplots_adjust(top=0.9)  # vertical space for main title
	plt.subplots_adjust(right=0.9)  # horizontal space for main title
	plt.text(-1.0, 0.5, 'Firing rate', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax1.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Total weight', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax3.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Early phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax5.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Late phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax7.transAxes) # add label for each measure
	plt.draw()  # draw plot but do not show
	plt.savefig(plot_folder + filename_av.replace('.txt', '.png'), bbox_inches='tight', dpi=300)
	plt.close()

	# output of averaged data to file
	f = open(plot_folder + filename_av, 'w')
	f.write("v_data =\r\n" + str(v_data) + "\r\n\r\n\r\n")
	f.write("h_core_inc =\r\n" + str(h_core_inc) + "\r\n\r\n\r\n")
	#f.write("h_core_out =\r\n" + str(h_core_out) + "\r\n\r\n\r\n")
	f.write("h_control_inc =\r\n" + str(h_control_inc) + "\r\n\r\n\r\n")
	f.write("z_core_inc =\r\n" + str(z_core_inc) + "\r\n\r\n\r\n")
	#f.write("z_core_out =\r\n" + str(z_core_out) + "\r\n\r\n\r\n")
	f.write("z_control_inc =\r\n" + str(z_control_inc) + "\r\n\r\n\r\n")
	f.close()


# plotTotalSubPopWeights
# Plots the total weight incoming from core neurons and from control neurons, respectively, from a complete weight matrix
# filename: name of the file to read the data from
# h_0: initial early-phase weight
# z_min: minimum late-phase weight
# z_max: maximum late-phase weight
# Nl_exc: number of excitatory neurons in one row/column
# s_CA: size of the cell assembly core
# title [optional]: main title of the figure
def plotTotalSubPopWeights(filename, h_0, z_min, z_max, Nl_exc, s_CA, title = "Total incoming weights from subpopulations"):

	# colormaps for h and z
	cmh0center = shiftedColorMap(cm.seismic, start=0, midpoint=0.5, stop=1.0, name='h0center')
	cmz = shiftedColorMap(cm.seismic, start=-1, midpoint=0, stop=1, name='asymmetric')

	# normalization factors
	h_max = 2*h_0
	w_max = h_max + h_0 * z_max

	# read weight matrix data
	connections, h, z, v = readWeightMatrixData(filename, Nl_exc)
	print("=========================\nData from: " + filename)

	# change filename
	filename_av = filename.replace('_net_', '_net_subpop_tot_')

	# normalize firing rates
	v_max = np.amax(v)
	if v_max == 0.0:
        	v_max = 1.0
	v_data = v.reshape(int(Nl_exc),int(Nl_exc)) / v_max

	# average synaptic weight from core neurons and from control neurons for all neurons
	core = np.arange(s_CA)

	h_core_inc = np.zeros((Nl_exc*Nl_exc))
	z_core_inc = np.zeros((Nl_exc*Nl_exc))
	#h_core_out = np.zeros((Nl_exc*Nl_exc))
	#z_core_out = np.zeros((Nl_exc*Nl_exc))
	h_control_inc = np.zeros((Nl_exc*Nl_exc))
	z_control_inc = np.zeros((Nl_exc*Nl_exc))

	for i in range(Nl_exc*Nl_exc):
		weights = earlyPhaseWeightsFromCore(i, connections, h, core)
		h_core_inc[i] = np.sum(weights)

		weights = latePhaseWeightsFromCore(i, connections, z, core)
		z_core_inc[i] = np.sum(weights)

		#weights = earlyPhaseWeightsToCore(i, connections, h, core)
		#h_core_out[i] = np.sum(weights)

		#weights = latePhaseWeightsToCore(i, connections, z, core)
		#z_core_out[i] = np.sum(weights)

		weights = earlyPhaseWeightsFromControl(i, connections, h, core)
		h_control_inc[i] = np.sum(weights)

		weights = latePhaseWeightsFromControl(i, connections, z, core)
		z_control_inc[i] = np.sum(weights)

	# reshape neuron arrays for plotting
	h_core_inc = h_core_inc.reshape(int(Nl_exc),int(Nl_exc))
	z_core_inc = z_core_inc.reshape(int(Nl_exc),int(Nl_exc))
	#h_core_out = h_core_out.reshape(int(Nl_exc),int(Nl_exc))
	#z_core_out = z_core_out.reshape(int(Nl_exc),int(Nl_exc))
	h_control_inc = h_control_inc.reshape(int(Nl_exc),int(Nl_exc))
	z_control_inc = z_control_inc.reshape(int(Nl_exc),int(Nl_exc))

	# plotting
	plt.figure(figsize=(6,12))
	grid = gs.GridSpec(4, 2)
	plt.rc('font', size=13)

	ax1 = plt.subplot(grid[0, :])
	p1 = ax1.imshow(v_data, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Oranges, origin='lower')
	ax1.set_title('v')
	ax1.xaxis.set_ticks([])
	ax1.yaxis.set_ticks([])
	plt.colorbar(p1, ticks=[0, 0.5, 1])

	ax3 = plt.subplot(grid[1, 0])
	p3 = ax3.imshow((h_core_inc + h_0*z_core_inc), vmin=0, vmax=1, interpolation='nearest', cmap=cm.Blues, origin='lower')
	ax3.set_title('Inc. w from core')
	ax3.xaxis.set_ticks([])
	ax3.yaxis.set_ticks([])
	plt.colorbar(p3)

	ax4 = plt.subplot(grid[1, 1])
	#p4 = ax4.imshow((h_core_out + h_0*z_core_out) / w_max, vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
	#ax4.set_title('Outg. w to core')
	p4 = ax4.imshow((h_control_inc + h_0*z_control_inc), vmin=0, vmax=1, interpolation='nearest', cmap=cm.Greens, origin='lower')
	ax4.set_title('Inc. w from control')
	ax4.xaxis.set_ticks([])
	ax4.yaxis.set_ticks([])
	plt.colorbar(p4)

	ax5 = plt.subplot(grid[2, 0])
	p5 = ax5.imshow(h_core_inc, interpolation='nearest', cmap=cm.Greens, origin='lower')
	ax5.set_title('Inc. h from core')
	plt.colorbar(p5)
	ax5.xaxis.set_ticks([])
	ax5.yaxis.set_ticks([])

	ax6 = plt.subplot(grid[2, 1])
	#p6 = ax6.imshow(logScale_h(h_core_out, h_0)/h_max, vmin=0, vmax=1, interpolation='nearest', cmap=cmap=cm.Greens, origin='lower')
	#ax6.set_title('Outg. h to core')
	p6 = ax6.imshow(h_control_inc, interpolation='nearest', cmap=cm.Greens, origin='lower')
	ax6.set_title('Inc. h from control')
	plt.colorbar(p6)
	ax6.xaxis.set_ticks([])
	ax6.yaxis.set_ticks([])

	ax7 = plt.subplot(grid[3, 0])
	p7 = ax7.imshow(z_core_inc, interpolation='nearest', cmap=cm.Greens, origin='lower')
	ax7.set_title('Inc. z from core')
	plt.colorbar(p7)
	ax7.xaxis.set_ticks([])
	ax7.yaxis.set_ticks([])

	ax8 = plt.subplot(grid[3, 1])
	#p8 = ax8.imshow(logScale_z(z_core_out, z_max), vmin=z_min, vmax=z_max, interpolation='nearest', cmap=cmap=cm.Greens, origin='lower')
	#ax8.set_title('⟨ Outgoing z to core⟩')
	p8 = ax8.imshow(z_control_inc, interpolation='nearest', cmap=cm.Greens, origin='lower')
	ax8.set_title('Inc. z from control')
	plt.colorbar(p8)
	ax8.xaxis.set_ticks([])
	ax8.yaxis.set_ticks([])

	# set general picture frame
	plt.tight_layout()
	plt.suptitle(title)  # main title
	plt.subplots_adjust(top=0.9)  # vertical space for main title
	plt.subplots_adjust(right=0.9)  # horizontal space for main title
	plt.text(-1.0, 0.5, 'Firing rate', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax1.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Total weight', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax3.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Early phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax5.transAxes) # add label for each measure
	plt.text(-0.3, 0.5, 'Late phase', horizontalalignment='center', verticalalignment='center', rotation='vertical', fontsize=14, transform=ax7.transAxes) # add label for each measure
	plt.draw()  # draw plot but do not show
	plt.savefig(plot_folder + filename_av.replace('.txt', '.png'), bbox_inches='tight', dpi=300)
	plt.close()

	# output of averaged data to file
	f = open(plot_folder + filename_av, 'w')
	f.write("v_data =\r\n" + str(v_data) + "\r\n\r\n\r\n")
	f.write("h_core_inc =\r\n" + str(h_core_inc) + "\r\n\r\n\r\n")
	#f.write("h_core_out =\r\n" + str(h_core_out) + "\r\n\r\n\r\n")
	f.write("h_control_inc =\r\n" + str(h_control_inc) + "\r\n\r\n\r\n")
	f.write("z_core_inc =\r\n" + str(z_core_inc) + "\r\n\r\n\r\n")
	#f.write("z_core_out =\r\n" + str(z_core_out) + "\r\n\r\n\r\n")
	f.write("z_control_inc =\r\n" + str(z_control_inc) + "\r\n\r\n\r\n")
	f.close()
	
# plotMinOverview
# Plots the results of a simulation of synaptic weights changing based on calcium dynamics
# data_file: data file containing the values of the membrane potential, weights, calcium amount, etc. over time
# col_neur: the first column containing data of the targeted neuron (the membrane potential, the next two columns contain membrane current and, if applicable, protein amount)
# col_syn: the first column containing data of the targeted synapse (the early-phase weight, the next two columns contain late-phase weight and calcium amount)
# h_0: initial weight
# theta_tag: tagging threshold
# theta_pro: protein synthesis threshold
# theta_p: potentiation threshold for Ca dynamics
# theta_d: depression threshold for Ca dynamics
# store_path [optional]: path to the resulting graphics file
def plotMinOverview(data_file, col_neur, col_syn, h_0, theta_tag, theta_pro, theta_p, theta_d, store_path = './traces.svg'):

	data_stacked = np.loadtxt(data_file)
	
	xlim_0 = None
	xlim_1 = None
	xlim_auto = True
	
	fig, axes = plt.subplots(nrows=3, ncols=1, sharex=False, figsize=(10, 10))

	# set axis labels for axes[0]
	axes[0].set_xlabel("Time (ms)")
	axes[0].set_ylabel("Synaptic weight (%)")
	axes[0].set_xlim(xlim_0, xlim_1, auto = xlim_auto)

	# convert x-axis values to ms
	data_stacked[:,0] *= 1000
	
	# plot data for axes[0]
	if col_syn >= 1: # if synapse data exist
		axes[0].plot(data_stacked[:,0], data_stacked[:,col_syn]/h_0*100, color="#800000", label='h', marker='None', zorder=10)
		axes[0].plot(data_stacked[:,0], (data_stacked[:,col_syn+1]+1)*100, color="#1f77b4", label='z', marker='None', zorder=9)
		axes[0].axhline(y=(theta_pro/h_0+1)*100, label='Protein thresh.', linestyle='-.', color="#dddddd", zorder=5)
		axes[0].axhline(y=(theta_tag/h_0+1)*100, label='Tag thresh.', linestyle='dashed', color="#dddddd", zorder=4)
		# total weight: color="#ff7f0e"
	
		# create legend for axes[0]
		axes[0].legend() #loc=(0.75,0.65)) #"center right")

	# set axis labels for axes[1] (and twin axis ax1twin)
	axes[1].set_xlabel("Time (ms)")
	axes[1].set_ylabel("Membrane potential (mV)")
	axes[1].set_xlim(xlim_0, xlim_1, auto = xlim_auto)
	ax1twin = axes[1].twinx()
	ax1twin.set_ylabel("Current (nA)")
	
	# plot data for axes[1] (and twin axis ax1twin)
	ax1g1 = axes[1].plot(data_stacked[:,0], data_stacked[:,col_neur], color="#ff0000", label='Membrane pot.', marker='None', zorder=10)
	ax1g2 = ax1twin.plot(data_stacked[:,0], data_stacked[:,col_neur+1], color="#ffee00", label='Membrane curr.', marker='None', zorder=10)
	
	# create common legend for axes[1] and ax1twin
	#fig.legend(loc=(0.72,0.45)) #"center right")
	#axes[1].legend(loc=(0.75,0.35)) #"center right")
	handles, labels = axes[1].get_legend_handles_labels()
	handles_twin, labels_twin = ax1twin.get_legend_handles_labels()
	axes[1].legend(handles + handles_twin, labels + labels_twin)
	
	# set axis labels for axes[2]
	axes[2].set_xlabel("Time (ms)")
	axes[2].set_ylabel("Calcium or protein amount")
	axes[2].set_xlim(xlim_0, xlim_1, auto = xlim_auto)
	
	# plot data for axes[2]
	if col_syn >= 1: # if synapse data exist
		axes[2].plot(data_stacked[:,0], data_stacked[:,col_syn+2], color="#c8c896", label='Ca', marker='None', zorder=8)
		axes[2].plot(data_stacked[:,0], data_stacked[:,col_neur+2], color="#008000", label='p', marker='None', zorder=7) # protein amount is also only shown if synapse data exist (which means that an E->E synapse is considered)
		axes[2].axhline(y=theta_p, label='LTP thresh.', linestyle='dashed', color="#969664", zorder=7)
		axes[2].axhline(y=theta_d, label='LTD thresh.', linestyle='dashed', color="#969696", zorder=6)
	
		# create legend for axes[2]
		axes[2].legend() #loc=(0.75,0.65)) #"center right")
	
	# save figure as vector graphics
	fig.savefig(store_path)
