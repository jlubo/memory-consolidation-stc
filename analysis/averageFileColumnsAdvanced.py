##############################################################################################
### Script to average data from the same columns in data files stored in different folders ###
##############################################################################################

### Copyright 2017-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
import os
from pathlib import Path

# averageFileColumns
# Averages specified data columns across data files located in directories which names contain a specific string
# and computes the standard deviation
# outname: name of the file to write the averaged data to
# rootpath: path in which to look for data folders
# protocol: string that the data folders have to contain
# suffix: suffix in the filename of data files to be read
# columns: list of numbers of the columns in the data file to be read and averaged (e.g., [1, 3] for first and third column)
# first_column_par [optional]: indicates if first column is to be treated as parameter (e.g., time) - it is then added regardless of 'columns'
# comment_line [optional]: if True, leaves out the first line
def averageFileColumns(outname, rootpath, protocol, suffix, columns, first_column_par=True, comment_line=False):
	print("Averaging columns " + str(columns) + " from files matching '*" + suffix + "' in folders of the protocol '" + protocol + "'...")
	sample_number = 0
	col_sep = '\t\t' # character(s) separating the columns

	# find the folders with the protocol in their name
	rawpaths = Path(rootpath)
	paths = np.array([str(x) for x in rawpaths.iterdir() if x.is_dir() and protocol in os.path.split(str(x))[1]])

	if paths.size == 0:
		raise FileNotFoundError("No folders found that contain the string '" + protocol + "' in their name.")
	print("According folders found:\n", paths)

	# read data and average
	# loop over directories
	for i in range(paths.size):

		# find the files with the suffix in their name
		subrawpaths = Path(paths[i])
		subpaths = np.array([str(x) for x in subrawpaths.iterdir() if suffix in str(x) and str(x).find(suffix) >= len(str(x))-len(suffix)])

		if subpaths.size == 0:
			raise FileNotFoundError("No files found matching '*" + suffix + "' in '" + paths[i] + "'.")

		print("According files found in '" + paths[i] + "':\n", subpaths)
		sample_number += subpaths.size

		# loop over files in each directory
		for j in range(subpaths.size):

			with open(subpaths[j]) as f:
				rawdata = f.read()

			rawdata = rawdata.split('\n')
			if comment_line:
				del rawdata[0] # leave out comment line
			if rawdata[-1] == "":
				del rawdata[-1] # delete empty line

			if i == 0 and j == 0: # first file found: read number of rows and create data arrays
				num_rows = len(rawdata)
				num_cols = len(columns)
				time = np.zeros(num_rows)
				data = np.zeros((num_rows, num_cols))
				data_var = np.zeros((num_rows, num_cols))
			elif num_rows != len(rawdata):
				raise IndexError("In '" + subpaths[j] + "': wrong number of rows: " + str(len(rawdata)-1) + " (" + str(num_rows) + " expected).")

			for k in range(num_rows):
				values = rawdata[k].split(col_sep)
				try:
					time[k] += np.double(values[0]) # read first/parameter column
				except:
					print("Error computing mean: in line " + str(k+1) + ", column 1\n\tin '" + subpaths[j] + "'.")
				for l in range(num_cols):
					try:
						data[k][l] += np.double(values[columns[l]-1]) # read data columns
					except:
						print("Error computing mean: in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'.")

			f.close()

	time = time / sample_number
	data = data / sample_number

	# read data and compute variance
	# loop over directories
	for i in range(paths.size):

		# loop over files in each directory
		for j in range(subpaths.size):

			with open(subpaths[j]) as f:
				rawdata = f.read()

			rawdata = rawdata.split('\n')

			if comment_line:
				del rawdata[0] # leave out comment line
			if rawdata[-1] == "":
				del rawdata[-1] # delete empty line

			for k in range(num_rows):
				values = rawdata[k].split(col_sep)
				for l in range(num_cols):
					try:
						data_var[k][l] += np.power(np.double(values[columns[l]-1])-data[k][l], 2) # read data columns
					except:
						print("Error computing variance: in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'.")

			f.close()

	data_stdev = np.sqrt(data_var / (sample_number - 1))

	# write averaged data
	fout = open(outname + '.txt', 'w')
	for k in range(num_rows):
		if first_column_par:
			fout.write(str(time[k]) + "\t\t")
		for l in range(num_cols):
			fout.write(str(data[k][l]) + "\t\t" + str(data_stdev[k][l]))
			if l < num_cols-1: # as long as last column is not yet reached
				 fout.write("\t\t")
			else: # after the last column
				 fout.write("\n")
	fout.close()

## average weight over time traces for standard plasticity-inducing protocols STET, WTET; SLFS; WLFS:
suffix = '_data.txt'
averageFileColumns('averaged_STET', '.', 'STET', suffix, [8, 9, 10], comment_line=True)
averageFileColumns('averaged_WTET', '.', 'WTET', suffix, [8, 9, 10], comment_line=True)
averageFileColumns('averaged_SLFS', '.', 'SLFS', suffix, [8, 9, 10], comment_line=True)
averageFileColumns('averaged_WLFS', '.', 'WLFS', suffix, [8, 9, 10], comment_line=True)

