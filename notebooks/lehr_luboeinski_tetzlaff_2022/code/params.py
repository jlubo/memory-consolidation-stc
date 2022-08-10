# Copyright 2022 Andrew Lehr
# The MIT LICENSE

import numpy as np
import os.path as path

'''
Paths to required directories and files
'''

# raw data paths
root = path.abspath(path.join(__file__ ,"../.."))
raw_data_dir = root + '/raw_data/'
spike_data_dir = raw_data_dir + 'spike_data/'
outgoing_weight_data_dir = raw_data_dir + 'outgoing_weight_data/8h/'
average_weight_data_dir = raw_data_dir + 'average_weight_data/'
average_rate_data_dir = raw_data_dir + 'average_rate_data/'

# average weights and firing rates (in raw data, but these are preprocessed)
avg_weights_file = average_weight_data_dir + 'Averaged_Weights_28810_all.csv'
avg_rates_file = average_rate_data_dir + 'Averaged_Firing_Rates_28810.csv'

# processed data directories, this is where analysed data gets stored
processed_data_dir = root + '/processed_data/'
event_data_dir = processed_data_dir + 'event_data/'
pca_data_dir = processed_data_dir + 'pca_data/'
temporal_stability_data_dir = processed_data_dir + 'temporal_stability_data/'
regression_data_dir = processed_data_dir + 'regression_data/'

# figures directory
figures = root + '/figures/'

# all processed data directories, used for check that these exist and if not
# they are created
all_directories = [processed_data_dir, 
                   event_data_dir,
                   pca_data_dir,
                   temporal_stability_data_dir,
                   regression_data_dir,
                   figures]

'''
Simulation parameters
'''

# basic parameters
dt = 0.0002 # ms
n_T = 100001
n_Pop = 2000
n_E = 1600
n_I = 400

# two sessions, before and after consolidation
sessions = ['10s', '8h']
n_sessions = len(sessions)

# all events
events = ['learn1','learn2','learn3','recall']

# stimulation time
stim_time = 0.100 # s

# times of events (in seconds)
t_events = {'10s': {'learn1': 10.0,
                    'learn2': 10.5,
                    'learn3': 11.0,
                    'recall': 20.0},
            '8h':  {'learn1': 10.0,
                    'learn2': 10.5,
                    'learn3': 11.0,
                    'recall': 28810.0}
            }

# range of ids for excitatory neurons
id_range = [0, 1599]
n_ids_sel = id_range[1] - id_range[0] + 1

# neuromodulator concentrations
min_nm = 0.00
max_nm = 0.18
step_nm = 0.02 
nmods = np.arange(min_nm, max_nm+step_nm, step_nm).round(2)
n_nmods = len(nmods)

# stimulation frequencies
min_f = 10 # Hz
max_f = 100 # Hz
step_f = 10 # Hz
freqs = np.arange(min_f, max_f+step_f, step_f)
n_freqs = len(freqs)

# number of trials per simulation batch
n_trials_per_batch = 10
n_batches = 5 # files separated into 5 folders

# total number of files
n_total_files = n_nmods * n_freqs * n_sessions


''' 
Figure Parameters
'''

colors = {0: ['black', 'black'],
          1: ['yellow', '#E7A42C'],
          2: ['green', '#4AA42C'],
          3: ['blue', '#2C69AC'],
          }


colors_dark = {0: ['black', 'black'],
               1: ['yellow', '#d1911f'],
               2: ['green', '#35851b'],
               3: ['blue', '#13559c'],
               }


# number of trials
n_total_trials = 50

# students t for 50 trials
students_t_50 = 2.01 
students_t_500 = 1.96   # int(self.n_trials*self.n_batches)*n_outputs

# number of bins, default 100, 1ms bins
n_bins = 100

# regression parameters
n_p_outs = 10
n_outputs = 10
p_outs = np.linspace(0.1,1,n_p_outs)