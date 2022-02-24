##############################################################################################
### Script to average data from the same columns in data files stored in different folders ###
##############################################################################################

### Copyright 2017-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

import numpy as np
import os
from pathlib import Path
from mergeRawData import *

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
	col_sep = '= ' #'\t\t' character(s) separating the columns

	# find the folders with the protocol in their name
	rawpaths = Path(rootpath)
	paths = np.array([str(x) for x in rawpaths.iterdir() if x.is_dir() and protocol in str(x)])

	if paths.size == 0:
		raise FileNotFoundError("No folders found that contain the string '" + protocol + "' in their name.")
	print("According folders found:\n", paths)

	# read data and average
	# loop over directories
	for i in range(paths.size):

		# find the files with the suffix in their name
		subrawpaths = Path(paths[i])
		subpaths = np.array([str(x) for x in subrawpaths.iterdir() if str(x).find(suffix) >= len(str(x))-len(suffix)])

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
				if len(values) < 2:
					values = [np.nan,np.nan] # to avoid problems reading descriptions
				try:
					time[k] += np.double(values[0]) # read first/parameter column
				except ValueError:
					pass#print("Computing mean: conversion error in line " + str(k+1) + ", column 1\n\tin '" + subpaths[j] + "'.")
				for l in range(num_cols):
					try:
						data[k][l] += np.double(values[columns[l]-1]) # read data columns
					except ValueError:
						pass#print("Computing mean: conversion error in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'.")

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
				#if len(values) < 2:
				#	values = [np.nan,np.nan] # to avoid problems reading descriptions
				for l in range(num_cols):
					try:
						data_var[k][l] += np.power(np.double(values[columns[l]-1])-data[k][l], 2) # read data columns
						
					except:
						pass#print("Computing variance: conversion error in line " + str(k+1) + ", column " + str(columns[l]) + "\n\tin '" + subpaths[j] + "'.")
					#except IndexError:
					#	print("INDEX ERROR")

			f.close()

	data_stdev = np.sqrt(data_var / (sample_number - 1))

	# write averaged data
	fout = open(outname + '.txt', 'w')

	for k in range(num_rows): ## ADAPTED

		if k >=4 and k <=7: # only need those four rows!
			for l in range(num_cols):
				fout.write(str(data[k][l]) + "\t" + str(data_stdev[k][l]))
				if (k+1) % 4 == 0 and l >= num_cols-1:  # after the last column and after 4 rows have been clutched together
					fout.write("\n")
				else: # as long as last column is not yet reached
					fout.write("\t")
	fout.close()

f = open("p_act_summary_temp_0names.txt", "w")
f.write("NOOVERLAP, A primed\n")
f.write("NOOVERLAP, B primed\n")
f.write("NOOVERLAP, C primed\n")
f.write("OVERLAP10, A primed\n")
f.write("OVERLAP10, B primed\n")
f.write("OVERLAP10, C primed\n")
f.write("OVERLAP10 no AC, no ABC, A primed\n")
f.write("OVERLAP10 no AC, no ABC, B primed\n")
f.write("OVERLAP10 no AC, no ABC, C primed\n")
f.write("OVERLAP10 no BC, no ABC, A primed\n")
f.write("OVERLAP10 no BC, no ABC, B primed\n")
f.write("OVERLAP10 no BC, no ABC, C primed\n")
f.close()

# 10 min
averageFileColumns("p_act_averaged_Aprimed", "2. Switching after 10 min/NOOVERLAP/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "2. Switching after 10 min/NOOVERLAP/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "2. Switching after 10 min/NOOVERLAP/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_10min_NOOVERLAP.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "2. Switching after 10 min/OVERLAP10/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "2. Switching after 10 min/OVERLAP10/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "2. Switching after 10 min/OVERLAP10/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_10min_OVERLAP10.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "2. Switching after 10 min/OVERLAP10 no AC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "2. Switching after 10 min/OVERLAP10 no AC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "2. Switching after 10 min/OVERLAP10 no AC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_10min_OVERLAP10_noAC_noABC.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "2. Switching after 10 min/OVERLAP10 no BC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "2. Switching after 10 min/OVERLAP10 no BC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "2. Switching after 10 min/OVERLAP10 no BC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_10min_OVERLAP10_noBC_noABC.txt", remove_raw=True, sep_str='\n')

os.system('cat "p_act_10min_NOOVERLAP.txt" > "p_act_summary_temp_10min.txt"')
os.system('cat "p_act_10min_OVERLAP10.txt" >> "p_act_summary_temp_10min.txt"')
os.system('cat "p_act_10min_OVERLAP10_noAC_noABC.txt" >> "p_act_summary_temp_10min.txt"')
os.system('cat "p_act_10min_OVERLAP10_noBC_noABC.txt" >> "p_act_summary_temp_10min.txt"')
os.system('rm -R -f "p_act_10min"*')
mergeRawData(".", "p_act_summary_temp_", "p_act_summary_10min.txt", remove_raw=False, sep_str='\t') # remove_raw=False to keep _temp_0names file
os.system('rm -f p_act_summary_temp_10min.txt')

# 1 h
averageFileColumns("p_act_averaged_Aprimed", "3. Switching after 1 h/NOOVERLAP/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "3. Switching after 1 h/NOOVERLAP/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "3. Switching after 1 h/NOOVERLAP/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_1h_NOOVERLAP.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "3. Switching after 1 h/OVERLAP10/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "3. Switching after 1 h/OVERLAP10/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "3. Switching after 1 h/OVERLAP10/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_1h_OVERLAP10.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "3. Switching after 1 h/OVERLAP10 no AC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "3. Switching after 1 h/OVERLAP10 no AC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "3. Switching after 1 h/OVERLAP10 no AC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_1h_OVERLAP10_noAC_noABC.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "3. Switching after 1 h/OVERLAP10 no BC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "3. Switching after 1 h/OVERLAP10 no BC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "3. Switching after 1 h/OVERLAP10 no BC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_1h_OVERLAP10_noBC_noABC.txt", remove_raw=True, sep_str='\n')

os.system('cat "p_act_1h_NOOVERLAP.txt" > "p_act_summary_temp_1h.txt"')
os.system('cat "p_act_1h_OVERLAP10.txt" >> "p_act_summary_temp_1h.txt"')
os.system('cat "p_act_1h_OVERLAP10_noAC_noABC.txt" >> "p_act_summary_temp_1h.txt"')
os.system('cat "p_act_1h_OVERLAP10_noBC_noABC.txt" >> "p_act_summary_temp_1h.txt"')
os.system('rm -R -f "p_act_1h"*')
mergeRawData(".", "p_act_summary_temp_", "p_act_summary_1h.txt", remove_raw=False, sep_str='\t')
os.system('rm -f p_act_summary_temp_1h.txt')

# 4 h
averageFileColumns("p_act_averaged_Aprimed", "4. Switching after 4 h/NOOVERLAP/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "4. Switching after 4 h/NOOVERLAP/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "4. Switching after 4 h/NOOVERLAP/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_4h_NOOVERLAP.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "4. Switching after 4 h/OVERLAP10/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "4. Switching after 4 h/OVERLAP10/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "4. Switching after 4 h/OVERLAP10/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_4h_OVERLAP10.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "4. Switching after 4 h/OVERLAP10 no AC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "4. Switching after 4 h/OVERLAP10 no AC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "4. Switching after 4 h/OVERLAP10 no AC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_4h_OVERLAP10_noAC_noABC.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "4. Switching after 4 h/OVERLAP10 no BC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "4. Switching after 4 h/OVERLAP10 no BC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "4. Switching after 4 h/OVERLAP10 no BC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_4h_OVERLAP10_noBC_noABC.txt", remove_raw=True, sep_str='\n')

os.system('cat "p_act_4h_NOOVERLAP.txt" > "p_act_summary_temp_4h.txt"')
os.system('cat "p_act_4h_OVERLAP10.txt" >> "p_act_summary_temp_4h.txt"')
os.system('cat "p_act_4h_OVERLAP10_noAC_noABC.txt" >> "p_act_summary_temp_4h.txt"')
os.system('cat "p_act_4h_OVERLAP10_noBC_noABC.txt" >> "p_act_summary_temp_4h.txt"')
os.system('rm -R -f "p_act_4h"*')
mergeRawData(".", "p_act_summary_temp_", "p_act_summary_4h.txt", remove_raw=False, sep_str='\t')
os.system('rm -f p_act_summary_temp_4h.txt')

# 7 h
averageFileColumns("p_act_averaged_Aprimed", "5. Switching after 7 h/NOOVERLAP/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "5. Switching after 7 h/NOOVERLAP/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "5. Switching after 7 h/NOOVERLAP/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_7h_NOOVERLAP.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "5. Switching after 7 h/OVERLAP10/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "5. Switching after 7 h/OVERLAP10/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "5. Switching after 7 h/OVERLAP10/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_7h_OVERLAP10.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "5. Switching after 7 h/OVERLAP10 no AC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "5. Switching after 7 h/OVERLAP10 no AC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "5. Switching after 7 h/OVERLAP10 no AC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_7h_OVERLAP10_noAC_noABC.txt", remove_raw=True, sep_str='\n')

averageFileColumns("p_act_averaged_Aprimed", "5. Switching after 7 h/OVERLAP10 no BC, no ABC/A", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Bprimed", "5. Switching after 7 h/OVERLAP10 no BC, no ABC/B", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
averageFileColumns("p_act_averaged_Cprimed", "5. Switching after 7 h/OVERLAP10 no BC, no ABC/C", "avalanche_statistics_0.01_10", "_CA_probabilities.txt", [2], first_column_par=False)
mergeRawData(".", "p_act_averaged_", "p_act_7h_OVERLAP10_noBC_noABC.txt", remove_raw=True, sep_str='\n')

os.system('cat "p_act_7h_NOOVERLAP.txt" > "p_act_summary_temp_7h.txt"')
os.system('cat "p_act_7h_OVERLAP10.txt" >> "p_act_summary_temp_7h.txt"')
os.system('cat "p_act_7h_OVERLAP10_noAC_noABC.txt" >> "p_act_summary_temp_7h.txt"')
os.system('cat "p_act_7h_OVERLAP10_noBC_noABC.txt" >> "p_act_summary_temp_7h.txt"')
os.system('rm -R -f "p_act_7h"*')
mergeRawData(".", "p_act_summary_temp_", "p_act_summary_7h.txt", remove_raw=True, sep_str='\t')

os.system('cp "./1. Priming/NOOVERLAP/A/"*/*"_PARAMS.txt" .')

