########################################################################
###  Routine to compute the frequency spectrum of spike raster data  ###
###                   from multiple data files                       ###
########################################################################

### Copyright 2019-2022 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: jannik.lubo[at]gmx.de

import numpy as np
import pandas as pd
import scipy
import scipy.fftpack
from scipy import pi
import matplotlib.pyplot as plt
import pylab
import os
from pathlib import Path

dt = 0.0002 # duration of one time bin

# search in current directory for a "*_spike_raster.txt" file
rawpaths = Path(".")
df = None
for x in sorted(rawpaths.iterdir()):

	full_path = str(x)
	tpath = os.path.split(full_path)[1] # take tail
	if "_spike_raster.txt" in tpath:
		print("Reading", tpath)
		df_new = pd.read_table(tpath, header=None, sep="\t\t", engine='python')
		if df is None:
			df = df_new
		else:
			df = df.append(df_new)

if df is None:
	print("No data found. Exiting...")
	exit()

# count the number of spikes per time bin (!!! not per 10 ms bin !!!)
spike_counts = df[df.columns[0]].value_counts().sort_index()
time_all = spike_counts.index.to_numpy()
spikes_whole_all = spike_counts.to_numpy()
print(time_all)
print(spikes_whole_all)

# Fast Fourier Transform
FFT = abs(scipy.fft.fft(spikes_whole_all))
freqs = scipy.fftpack.fftfreq(spikes_whole_all.size, dt)

pylab.subplot(211)
pylab.xlabel('Time (s)')
pylab.ylabel('Spikes')
pylab.plot(time_all, spikes_whole_all, '-')
pylab.subplot(212)
pylab.xlabel('Frequency (Hz)')
pylab.ylabel('Amplitude')
plt.semilogy(freqs, FFT, '-', color="darkgreen")
pylab.show()

#plt.savefig("frequencyAnalysis.svg")
