#########################################################################################
### Routine to determine the distribution of spikes per time bin in the whole network ###
###                and in different assemblies from multiple data files               ###
#########################################################################################

### Copyright 2019-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import pandas as pd
from scipy import optimize
import matplotlib.pyplot as plt
import pylab
import os
from pathlib import Path

# search in current directory for a "*_CA_time_series.txt" file (created by assemblyAvalancheStatistics.py)
rawpaths = Path(".")
df = None
for x in sorted(rawpaths.iterdir()):

	full_path = str(x)
	tpath = os.path.split(full_path)[1] # take tail
	if "_CA_time_series.txt" in tpath:
		print("Reading", tpath)
		df_new = pd.read_table(tpath, header=None, sep="\t\t", engine='python')
		if df is None:
			df = df_new
		else:
			df = df.append(df_new)

if df is None:
	print("No data found. Exiting...")
	exit()

num_pars = 2 # number of fit parameters (2: standard power-law fit, 3: power-law fit with shift)
fit_range_start = 10 # first datapoint to be used for fit
fit_range_end = 200 # last datapoint to be used to fit

# creating histograms
valcA = df[df.columns[2]].value_counts()
print("index(A) =", valcA.index)
print("values(A) =", valcA.values)
print("quantiles for A:\n", df[df.columns[2]].quantile([.25, .5, .99]))

valcB = df[df.columns[3]].value_counts()
print("index(B) =", valcB.index)
print("values(B) =", valcB.values)
print("quantiles for B:\n", df[df.columns[3]].quantile([.25, .5, .99]))

valcC = df[df.columns[4]].value_counts()
print("index(C) =", valcC.index)
print("values(C) =", valcC.values)
print("quantiles for C:\n", df[df.columns[4]].quantile([.25, .5, .99]))

df_all = df[df.columns[5]]
valcAll = df_all.value_counts()
print("index(all) =", valcAll.index)
print("values(all) =", valcAll.values)
print("quantiles for all:\n", df_all.quantile([.25, .5, .99]))

# plotting and fitting preparations
fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
fig.set_size_inches(10,10)
ax1 = plt.subplot2grid((6, 4), (0, 0), colspan=1, rowspan=6)
ax2 = plt.subplot2grid((6, 4), (0, 1), colspan=1, rowspan=6, sharey=ax1)
ax3 = plt.subplot2grid((6, 4), (0, 2), colspan=1, rowspan=6, sharey=ax1)
ax4 = plt.subplot2grid((6, 4), (0, 3), colspan=1, rowspan=6, sharey=ax1)
fig.suptitle('Distribution of spikes per assembly in one bin')
if num_pars == 2:
	fitfunc = lambda par, x: par[0] / (x**par[1])
	par0 = [12000, 1.55] # initial guess
else:
	fitfunc = lambda par, x: par[0] / (x**par[1] + par[2])
	par0 = [7000, 1.55, 1] # initial guess
errfunc = lambda par, x, y: fitfunc(par, x) - y # distance to the target function

# fitting histogram for A
if num_pars == 2:
	x_A = valcA.index[np.logical_and(valcA.index>=fit_range_start, valcA.index<=fit_range_end)]
	y_A = valcA.values[np.logical_and(valcA.index>=fit_range_start, valcA.index<=fit_range_end)]
else:
	x_A = valcA.index
	y_A = valcA.values
fit_result_A = optimize.least_squares(errfunc, par0[:], args=(x_A, y_A))
print("Fit parameters for A:", fit_result_A.x, ", success:", fit_result_A.success, \
      ", optimality:", fit_result_A.optimality, ", residual:", np.sum(np.square(fit_result_A.fun))/fit_result_A.fun.shape[0])

# fitting histogram for B
if num_pars == 2:
	x_B = valcB.index[np.logical_and(valcB.index>=fit_range_start, valcB.index<=fit_range_end)]
	y_B = valcB.values[np.logical_and(valcB.index>=fit_range_start, valcB.index<=fit_range_end)]
else:
	x_B = valcB.index
	y_B = valcB.values
fit_result_B = optimize.least_squares(errfunc, par0[:], args=(x_B, y_B))
print("Fit parameters for B:", fit_result_B.x, ", success:", fit_result_B.success, \
      ", optimality:", fit_result_B.optimality, ", residual:", np.sum(np.square(fit_result_B.fun))/fit_result_B.fun.shape[0])

# fitting histogram for C
if num_pars == 2:
	x_C = valcC.index[np.logical_and(valcC.index>=fit_range_start, valcC.index<=fit_range_end)]
	y_C = valcC.values[np.logical_and(valcC.index>=fit_range_start, valcC.index<=fit_range_end)]
else:
	x_C = valcC.index
	y_C = valcC.values
fit_result_C = optimize.least_squares(errfunc, par0[:], args=(x_C, y_C))
print("Fit parameters for C:", fit_result_C.x, ", success:", fit_result_C.success, \
      ", optimality:", fit_result_C.optimality, ", residual:", np.sum(np.square(fit_result_C.fun))/fit_result_C.fun.shape[0])

# fitting histogram for all
if num_pars == 2:
	x_all = valcAll.index[np.logical_and(valcAll.index>=fit_range_start, valcAll.index<=fit_range_end)]
	y_all = valcAll.values[np.logical_and(valcAll.index>=fit_range_start, valcAll.index<=fit_range_end)]
else:
	x_all = valcAll.index
	y_all = valcAll.values
fit_result_all = optimize.least_squares(errfunc, par0[:], args=(x_all, y_all))
print("Fit parameters for C:", fit_result_all.x, ", success:", fit_result_all.success, \
      ", optimality:", fit_result_all.optimality, ", residual:", np.sum(np.square(fit_result_all.fun))/fit_result_all.fun.shape[0])

# plotting
range_x = np.linspace(1, np.amax([valcA.index.max(), valcB.index.max(), valcC.index.max()]), 2000)
range_x_all = np.linspace(1, valcAll.index.max(), 2000)

ax1.text(1, 0.35, r"$\gamma$ = " + str(round(fit_result_A.x[1],2)))
ax2.text(1, 0.35, r"$\gamma$ = " + str(round(fit_result_B.x[1],2)))
ax3.text(1, 0.35, r"$\gamma$ = " + str(round(fit_result_C.x[1],2)))
ax4.text(1, 0.35, r"$\gamma$ = " + str(round(fit_result_all.x[1],2)))

ax1.set_xlabel('Number of spikes')
ax1.set_ylabel('Frequency')
ax2.set_xlabel('Number of spikes')
ax3.set_xlabel('Number of spikes')
ax4.set_xlabel('Number of spikes')
ax1.set_ylim(0.1, 1e5)

ax1.loglog([df[df.columns[2]].quantile([.99]), df[df.columns[2]].quantile([.99])], [0.1, 1e5], "-", color="#cccccc", dashes=[6, 2])
plot1d,plot1f = ax1.loglog(valcA.index, valcA.values, "ro", range_x, fitfunc(fit_result_A.x, range_x), "-", color="#004586", label='A') # plot of the data and the fit
ax2.loglog([df[df.columns[3]].quantile([.99]), df[df.columns[3]].quantile([.99])], [0.1, 1e5], "-", color="#cccccc", dashes=[6, 2])
plot2d,plot2f = ax2.loglog(valcB.index, valcB.values, "ro", range_x, fitfunc(fit_result_B.x, range_x), "-", color="#ff420e", label='B') # plot of the data and the fit
ax3.loglog([df[df.columns[4]].quantile([.99]), df[df.columns[4]].quantile([.99])], [0.1, 1e5], "-", color="#cccccc", dashes=[6, 2])
plot3d,plot3f = ax3.loglog(valcC.index, valcC.values, "ro", range_x, fitfunc(fit_result_C.x, range_x), "-", color="#ffd320", label='C') # plot of the data and the fit
ax4.loglog([df[df.columns[5]].median(), df[df.columns[5]].median()], [0.1, 1e5], "-", color="#333333", dashes=[6, 2])
plot4d,plot4f = ax4.loglog(valcAll.index, valcAll.values, "ro", range_x_all, fitfunc(fit_result_all.x, range_x_all), "-", color="#111111", label='Overall') # plot of the data and the fit

#plt.savefig("numberOfSpikesInBins.svg")
plt.show()
