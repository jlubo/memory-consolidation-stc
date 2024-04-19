##################################################################################################
### Script that averages over trials of Q and MI data, employing 'extractParamsQMIfromSpikes.py' #
##################################################################################################

### Copyright 2023 Jannik Luboeinski
### licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0)
### Contact: mail[at]jlubo.net

import numpy as np
import pandas as pd
import os
import sys
import extractParamsQMIfromSpikes

# extractAverageQMI
# Extracts recall performance measures Q and MI from data, averages over trials
# dir_path: the directory where the data is located
# col_sep [optional]: characters separating columns in the data (and output) files
# N_exc [optional]: the number of excitatory neurons
# N_inh [optional]: the number of inhibitory neurons
# pattern_sizes [optional]: list of pattern sizes to be considered
# output [optional]: decreasing level of console output ("full", "final", or "none")
def extractAverageQMI(dir_path, col_sep = "\x20", N_exc = 1600, N_inh = 400, pattern_sizes = [100, 150, 200, 250],
                      output = "full"):
	
	# ---------------------------------------------------------------------------------------
	# define output files
	params_Q_MI_file = os.path.join(dir_path, "Params_Q_MI.txt")
	averaged_Q_MI_file = os.path.join(dir_path, "data_Q_MI.txt")
	fout_1 = open(params_Q_MI_file, "w")
	fout_2 = open(averaged_Q_MI_file, "w")

	# ---------------------------------------------------------------------------------------
	# extract Q and MI measures from data
	extractParamsQMIfromSpikes.extractRecursion(dir_path, fout_1, col_sep, N_exc, N_inh, output = output)
	fout_1.close()

	# ---------------------------------------------------------------------------------------
	# loop over pattern sizes
	for pattern_size in pattern_sizes:

			# ---------------------------------------------------------------------------------------
			# select rows for specific pattern size and diffusivity value
			df_all_trials = pd.read_table(params_Q_MI_file, header=None, sep=col_sep, engine='python')
			df_selected_trials = df_all_trials[df_all_trials[1] == pattern_size]
			num_trials = len(df_selected_trials.index)

			# ---------------------------------------------------------------------------------------
			# compute mean and std. dev. across trials (for firing rates, Q, MI)
			nu_exc_mean = df_selected_trials.loc[:,6].mean(axis=0)
			nu_exc_sd = df_selected_trials.loc[:,6].std(axis=0)
			nu_as_mean = df_selected_trials.loc[:,8].mean(axis=0)
			nu_as_sd = df_selected_trials.loc[:,8].std(axis=0)
			nu_ans_mean = df_selected_trials.loc[:,10].mean(axis=0)
			nu_ans_sd = df_selected_trials.loc[:,10].std(axis=0)
			nu_ctrl_mean = df_selected_trials.loc[:,12].mean(axis=0)
			nu_ctrl_sd = df_selected_trials.loc[:,12].std(axis=0)
			nu_inh_mean = df_selected_trials.loc[:,14].mean(axis=0)
			nu_inh_sd = df_selected_trials.loc[:,14].std(axis=0)
			Q_mean = df_selected_trials.loc[:,16].mean(axis=0)
			Q_sd = df_selected_trials.loc[:,16].std(axis=0)
			MI_mean = df_selected_trials.loc[:,18].mean(axis=0)
			MI_sd = df_selected_trials.loc[:,18].std(axis=0)
			selfMIL_mean = df_selected_trials.loc[:,19].mean(axis=0)
			selfMIL_sd = df_selected_trials.loc[:,19].std(axis=0)
			if output in ["full", "final"]:
				print(f"Pattern size {pattern_size} ({num_trials} trials):")
				print(f"  nu_exc = {nu_exc_mean} +- {nu_exc_sd}")
				print(f"  nu_as = {nu_as_mean} +- {nu_as_sd}")
				print(f"  nu_ans = {nu_ans_mean} +- {nu_ans_sd}")
				print(f"  nu_ctrl = {nu_ctrl_mean} +- {nu_ctrl_sd}")
				print(f"  nu_inh = {nu_inh_mean} +- {nu_inh_sd}")
				print(f"  Q = {Q_mean} +- {Q_sd}")
				print(f"  MI = {MI_mean} +- {MI_sd}")
				print(f"  selfMIL = {selfMIL_mean} +- {selfMIL_sd}")
			fout_2.write(str(pattern_size) + col_sep +
					     str(Q_mean) + col_sep + str(Q_sd) + col_sep +
					     str(MI_mean) + col_sep + str(MI_sd) + col_sep +
					     str(selfMIL_mean) + col_sep + str(selfMIL_sd) + col_sep +
			             str(num_trials) + "\n")
	fout_2.close()

# ---------------------------------------------------------------
# main
if __name__ == '__main__':
	extractAverageQMI(".")

