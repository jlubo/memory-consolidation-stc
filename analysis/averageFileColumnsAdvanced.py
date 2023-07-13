##############################################################################################
### Script to average data from the same columns in data files stored in different folders ###
###  -- v1.2 --                                                                            ###
##############################################################################################

### Copyright 2017-2023 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import os
from pathlib import Path

# averageFileColumns
# Averages specified data columns from multiple files that may be located in different folders, and computes the standard deviation
# - outname: name of the file to write the averaged data to
# - rootpath: path in which to look for data folders
# - protocol: string that data folders have to contain in their name
# - req_contains: list of strings that data files have to contain anywhere in their name
# - req_suffix: string that the suffix of the name of a data file has to match
# - columns: list of numbers of columns to be read from the data files and to be averaged over (e.g., [2, 3] for second and third column)
# - columns_sum [optional]: list of numbers of columns to be read from the data files and to be summed and averaged over (e.g., [2, 3] for the sum of the second and third column)
# - first_column_anyway [optional]: indicates if first column is to be treated as parameter (e.g., time) - it is then added regardless of 'columns'
# - skip_first_line [optional]: if True, leaves out the first line (needed if it contains column headings, for example)
# - col_sep [optional]: character(s) separating the columns
def averageFileColumns(outname, rootpath, protocol, req_contains, req_suffix, columns, \
                       columns_sum = [], first_column_anyway=True, skip_first_line=False, col_sep = '\t\t'):
       
	description = "Task: averaging columns " + str(columns)
	if columns_sum:
		description += ", and sum over columns " + str(columns_sum)
	description += ", from files matching ["
	for x in req_contains:
		description += "'*" + x + "*', "
	description += "'*" + req_suffix + "'] in folders of the protocol '" + protocol + "'"
	print(description)
	
	sample_number = 0

	# find the folders with the protocol in their name
	rawpaths = Path(rootpath)
	paths = np.sort(np.array([str(x) for x in rawpaths.iterdir() if x.is_dir() and protocol in os.path.split(str(x))[1]]))

	if paths.size == 0:
		raise FileNotFoundError("No folders found that contain the string '" + protocol + "' in their name.")
	print("Found", len(paths), "matching folder" + ("s" if len(paths) > 1 else "") + ":\n", paths)

	# read data and average
	# loop over directories
	for i in range(paths.size):

		# find the files that have the strings from the 'req_contains' list in their name
		subrawpaths = Path(paths[i])
		subpaths_list = []
		for x in subrawpaths.iterdir():
			filename = str(x)
			if (req_suffix in filename and filename.find(req_suffix) >= len(filename)-len(req_suffix)) or req_suffix == "": # check if suffix matches
				use = True
				for req_str in req_contains: # loop over all requested strings
					if not req_str in filename:  # check if filename does not contain string
						use = False
				if use:
					subpaths_list.append(filename)
		subpaths = np.sort(np.array(subpaths_list))
		
		if subpaths.size == 0:
			raise FileNotFoundError("No files found matching '*" + req_suffix + "' in '" + paths[i] + "'.")

		print("Found", len(subpaths), "matching file" + ("s" if len(subpaths) > 1 else "") + " in '" + paths[i] + "':\n", subpaths)
		sample_number += subpaths.size

		# loop over files in each directory
		for j in range(subpaths.size):

			with open(subpaths[j]) as f:
				rawdata = f.read()

			rawdata = rawdata.split('\n')
			if skip_first_line:
				del rawdata[0] # leave out the first line
			if rawdata[-1] == "":
				del rawdata[-1] # delete empty line

			if i == 0 and j == 0: # first file found: read number of rows and create data arrays
				num_rows = len(rawdata)
				num_cols = len(columns)
				time = np.zeros(num_rows)
				if columns_sum:
					num_cols_eff = num_cols + 1
				else:
					num_cols_eff = num_cols
				data = np.zeros((num_rows, num_cols_eff))
				data_var = np.zeros((num_rows, num_cols_eff))
			elif num_rows != len(rawdata):
				raise IndexError("In '" + subpaths[j] + "': wrong number of rows: " + str(len(rawdata)-1) + " (" + str(num_rows) + " expected).")

			# loop over lines
			for k in range(num_rows):
				
				# time/parameter to average
				values = rawdata[k].split(col_sep)
				try:
					time[k] += float(values[0]) # read first/parameter column
				except:
					raise ValueError("Error computing mean: in line " + str(k+1) + ", column 1\n\tin '" + subpaths[j] + "'. Is 'col_sep' set correctly?")
				
				# loop over data columns to average
				for l in range(num_cols):
					try:
						data[k][l] += float(values[columns[l]-1]) # read data columns
					except:
						raise ValueError("Error computing mean: in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'. Is 'col_sep' set correctly?")
				
				# loop over data columns to sum and average
				if columns_sum:
					try:
						tmp_summed_data = 0
						for l in range(len(columns_sum)):
							tmp_summed_data += float(values[columns_sum[l]-1]) # read data columns and add them to yield a new random variable
						data[k][num_cols] += tmp_summed_data
					except:
						raise ValueError("Error computing mean: in line " + str(k+1) + ", sum over columns " + str(columns_sum) + "\n\tin '" + subpaths[j] + "'. Is 'col_sep' set correctly?")

			f.close()

	time /= sample_number
	data /= sample_number

	# read data and compute variance
	# loop over directories
	for i in range(paths.size):

		# loop over files in each directory
		for j in range(subpaths.size):

			with open(subpaths[j]) as f:
				rawdata = f.read()

			rawdata = rawdata.split('\n')

			if skip_first_line:
				del rawdata[0] # leave out first line
			if rawdata[-1] == "":
				del rawdata[-1] # delete empty line

			# loop over lines
			for k in range(num_rows):
				values = rawdata[k].split(col_sep)

				# loop over data columns to average
				for l in range(num_cols):
					try:
						data_var[k][l] += np.power(float(values[columns[l]-1]) - data[k][l], 2) # read data columns and use the mean to compute the variance
					except:
						raise ValueError("Error computing variance: in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'.")

				# loop over data columns to sum and average
				if columns_sum:
					try:
						tmp_summed_data = 0
						for l in range(len(columns_sum)):
							tmp_summed_data += float(values[columns_sum[l]-1]) # read data columns and add them to yield a new random variable
						data_var[k][num_cols] += np.power(tmp_summed_data - data[k][num_cols], 2) # use the mean to compute the variance
					except:
						raise ValueError("Error computing variance: in line " + str(k+1) + ", sum over columns " + str(columns_sum) + "\n\tin '" + subpaths[j] + "'.")

			f.close()

	data_var /= (sample_number - 1)
	data_stdev = np.sqrt(data_var)

	# write averaged data
	fout = open(outname, 'w')
	for k in range(num_rows):
		if first_column_anyway:
			fout.write(str(time[k]) + col_sep)
		for l in range(num_cols_eff):
			fout.write(str(data[k][l]) + col_sep + str(data_stdev[k][l]))
			if l < num_cols_eff-1: # as long as last column is not yet reached
				 fout.write(col_sep)
			else: # after the last column
				 fout.write("\n")
	fout.close()
	
#####################################
### examples:
#if __name__ == '__main__':
''' 
## average weight over time traces for standard plasticity-inducing protocols STET, WTET; SLFS; WLFS (Arbor implementation):
suffix = '_traces.txt'
averageFileColumns('averaged_STET.txt', '.', 'data_STET', [], suffix, [2, 3, 4, 5], [2, 3], skip_first_line=False, col_sep=' ')
averageFileColumns('averaged_WTET.txt', '.', 'data_WTET', [], suffix, [2, 3, 4, 5], [2, 3], skip_first_line=False, col_sep=' ')
averageFileColumns('averaged_SLFS.txt', '.', 'data_SLFS', [], suffix, [2, 3, 4, 5], [2, 3], skip_first_line=False, col_sep=' ')
averageFileColumns('averaged_WLFS.txt', '.', 'data_WLFS', [], suffix, [2, 3, 4, 5], [2, 3], skip_first_line=False, col_sep=' ')

## average weight over time traces for standard plasticity-inducing protocols STET, WTET; SLFS; WLFS (stand-alone implementation):
suffix = '_data.txt'
averageFileColumns('averaged_STET.txt', '.', 'STET', [], suffix, [8, 9, 10], skip_first_line=True)
averageFileColumns('averaged_WTET.txt', '.', 'WTET', [], suffix, [8, 9, 10], skip_first_line=True)
averageFileColumns('averaged_SLFS.txt', '.', 'SLFS', [], suffix, [8, 9, 10], skip_first_line=True)
averageFileColumns('averaged_WLFS.txt', '.', 'WLFS', [], suffix, [8, 9, 10], skip_first_line=True)

## average weight distributions in nm=* folders:
suffix = "_totweight_dist_28810.0.txt"
averageFileColumns('averaged_0.02_28810.0.txt', '.', 'nm=0.02', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.08_28810.0.txt', '.', 'nm=0.08', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.14_28810.0.txt', '.', 'nm=0.14', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.20_28810.0.txt', '.', 'nm=0.20', [], suffix, [2, 3, 4, 5])
suffix = "_totweight_dist_20.0.txt"
averageFileColumns('averaged_0.02_20.0.txt', '.', 'nm=0.02', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.08_20.0.txt', '.', 'nm=0.08', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.14_20.0.txt', '.', 'nm=0.14', [], suffix, [2, 3, 4, 5])
averageFileColumns('averaged_0.20_20.0.txt', '.', 'nm=0.20', [], suffix, [2, 3, 4, 5])

## average activity distributions in nm=* folders:
suffix = "_act_dist_28810.1.txt"
averageFileColumns('averaged_fr_0.02_28810.1.txt', '.', 'nm=0.02', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.08_28810.1.txt', '.', 'nm=0.08', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.14_28810.1.txt', '.', 'nm=0.14', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.20_28810.1.txt', '.', 'nm=0.20', [], suffix, [2, 3])
suffix = "_act_dist_20.0.txt"
averageFileColumns('averaged_fr_0.02_20.0.txt', '.', 'nm=0.02', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.08_20.0.txt', '.', 'nm=0.08', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.14_20.0.txt', '.', 'nm=0.14', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.20_20.0.txt', '.', 'nm=0.20', [], suffix, [2, 3])
suffix = "_act_dist_20.1.txt"
averageFileColumns('averaged_fr_0.02_20.1.txt', '.', 'nm=0.02', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.08_20.1.txt', '.', 'nm=0.08', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.14_20.1.txt', '.', 'nm=0.14', [], suffix, [2, 3])
averageFileColumns('averaged_fr_0.20_20.1.txt', '.', 'nm=0.20', [], suffix, [2, 3])

## average mean total weights from A, B, or C protocol:
contains = ["cores_mean_tot_weights"]
suffix = ".txt"
averageFileColumns('NOOVERLAP_A.txt', 'NOOVERLAP', 'A', contains, suffix, [2, 4, 6])
averageFileColumns('NOOVERLAP_B.txt', 'NOOVERLAP', 'B', contains, suffix, [2, 4, 6])
averageFileColumns('NOOVERLAP_C.txt', 'NOOVERLAP', 'C', contains, suffix, [2, 4, 6])
averageFileColumns('OVERLAP10_A.txt', 'OVERLAP10', 'A', contains, suffix, [2, 4, 6])
averageFileColumns('OVERLAP10_B.txt', 'OVERLAP10', 'B', contains, suffix, [2, 4, 6])
averageFileColumns('OVERLAP10_C.txt', 'OVERLAP10', 'C', contains, suffix, [2, 4, 6])
'''

