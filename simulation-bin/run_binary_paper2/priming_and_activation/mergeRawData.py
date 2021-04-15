###################################################################
### Script to merge data from multiple files into a single file ###
###################################################################

### Copyright 2020-2021 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)

from pathlib import Path
import os

######################################
# mergeRawData
# Looks in a specified directory for files with a certain substring in the filename and merges them
# (merging the content of the lines) to a single file
# rootpath: relative path to the output directory
# substr: string that the filename of files to be merged has to contain
# output_file: name of the output file
# remove_raw [optional]: removes the raw data files
# sep_str [optional]: the character or string by which to separate the lines in the output file
def mergeRawData(rootpath, substr, output_file, remove_raw=False, sep_str='\t\t'):

	path = Path(rootpath)
	num_rows = -1
	all_data = []

	for x in sorted(path.iterdir()): # loop through files in the output directory
		x_str = str(x)
		if not x.is_dir() and substr in x_str:

			f = open(x_str)
			single_trial_data = f.read()
			f.close()

			single_trial_data = single_trial_data.split('\n')

			if single_trial_data[-1] == "":
				del single_trial_data[-1] # delete empty line

			if len(single_trial_data) != num_rows:
				if num_rows == -1:
					num_rows = len(single_trial_data)
					all_data = single_trial_data
				else:
					raise Exception("Wrong number of rows encountered in: " + x_str)
			else:
				for i in range(num_rows):
					all_data[i] += sep_str + single_trial_data[i]

			if remove_raw:
				os.remove(x_str)

	fout = open(os.path.join(rootpath, output_file), "w")
	for i in range(num_rows):
		fout.write(all_data[i] + '\n')
	fout.close()
