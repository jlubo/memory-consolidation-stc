##############################################################################################
### Script to average data from the same columns in data files stored in different folders ###
##############################################################################################

### Copyright 2017-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
from pathlib import Path

# averageFileColumns
# Averages data from three columns across data files located in directories which names contain the string <protocol>
# outname: name of the file to write the averaged data to
# skip: number of characters in the filename of the data file before the suffix
# suffix: suffix of the data file to be read
# protocol: string that the path to the file has to contain
# col: number of the first column in the data file to be read
# col2: number of the second column in the data file to be read
# col3: number of the third column in the data file to be read
def averageFileColumns(outname, skip, suffix, protocol, col, col2, col3):
	print("Averaging columns " + str(col) + "," + str(col2) + "," + str(col3) + " in *" + suffix + "for protocol " + protocol + "...")

	# gather folder names
	rawpaths = Path('.')
	paths = np.array([str(x) for x in rawpaths.iterdir() if x.is_dir() and protocol in str(x)])

	# read data and average
	for fn in range(paths.size):

		filename = paths[fn] + "/" + paths[fn][0:skip] + suffix # create the filename using the directory name
		with open(filename) as f:
			rawdata = f.read()

		rawdata = rawdata.split('\n')
		del rawdata[0] # leave out comment line

		if fn == 0: # read number of rows and create data arrays
			num_rows = len(rawdata)-1
			time = np.zeros(num_rows)
			data = np.zeros(num_rows)
			data2 = np.zeros(num_rows)
			data3 = np.zeros(num_rows)
			data_var = np.zeros(num_rows)
			data2_var = np.zeros(num_rows)
			data3_var = np.zeros(num_rows)
		elif num_rows != len(rawdata)-1:
			print("Error: wrong number of rows")

		for i in range(num_rows):
		    values = np.double(rawdata[i].split('\t\t'))
		    time[i] += values[0]
		    data[i] += values[col-1]
		    data2[i] += values[col2-1]
		    data3[i] += values[col3-1]

		f.close()

	time = time / paths.size
	data = data / paths.size
	data2 = data2 / paths.size
	data3 = data3 / paths.size

	# read data and compute variance
	for fn in range(paths.size):

		filename = paths[fn] + "/" + paths[fn][0:skip] + suffix # create the filename using the directory name
		with open(filename) as f:
			rawdata = f.read()

		rawdata = rawdata.split('\n')
		del rawdata[0] # leave out comment line

		for i in range(num_rows):
		    values = np.double(rawdata[i].split('\t\t'))
		    data_var[i] += np.power(values[col-1]-data[i], 2)
		    data2_var[i] += np.power(values[col2-1]-data2[i], 2)
		    data3_var[i] += np.power(values[col3-1]-data3[i], 2)

		f.close()

	data_var = np.sqrt(data_var) / (paths.size - 1)
	data2_var = np.sqrt(data2_var) / (paths.size - 1)
	data3_var = np.sqrt(data3_var) / (paths.size - 1)

	# write averaged data
	fout = open(outname + '.txt', 'w')
	for i in range(num_rows):
		fout.write(str(time[i]) + "\t\t" + str(data[i]) + "\t\t" + str(data_var[i]) + "\t\t" + \
		                                   str(data2[i]) + "\t\t" + str(data2_var[i]) + "\t\t" + \
		                                   str(data3[i]) + "\t\t" + str(data3_var[i]) + "\r\n")
	fout.close()

averageFileColumns('averaged_STET', 17, '_data.txt', 'STET', 8, 9, 10)
averageFileColumns('averaged_WTET', 17, '_data.txt', 'WTET', 8, 9, 10)
averageFileColumns('averaged_SLFS', 17, '_data.txt', 'SLFS', 8, 9, 10)
averageFileColumns('averaged_WLFS', 17, '_data.txt', 'WLFS', 8, 9, 10)
