# Copyright 2022 Andrew Lehr
# The MIT LICENSE

import os
import datetime
import warnings
import numpy as np
import seaborn
from sklearn.linear_model import Ridge
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from collections import defaultdict
import re
import scipy.stats as stats
from statsmodels.formula.api import ols
import statsmodels.api as sm
import pandas as pd
from tabulate import tabulate
import os, shutil, pathlib
import _pickle as cPickle
import time
from utilityFunctions import event_data_obj


class stc_data_obj:
    '''
    This class contains methods to load the raw data, extract revelant
    information, and do required analysis steps to produce figures from 
    Lehr, Luboeinski, & Tetzlaff 2022.
    
    See the accompanying jupyter notebook to reproduce figures.
    '''
    
    # time of creation and data file    
    now = datetime.datetime.now()
    now_str = now.strftime('%Y%m%d-%H%M%S')
    
    def __init__(self, params):   
        ''' 
        initializer
        '''
        self.raw_data_dir = params.raw_data_dir
        self.spike_data_dir = params.spike_data_dir
        self.outgoing_weight_data_dir = params.outgoing_weight_data_dir
        self.processed_data_dir = params.processed_data_dir
        self.regression_data_dir = params.regression_data_dir
        self.average_weight_data_dir = params.average_weight_data_dir
        self.avg_weights_file = params.avg_weights_file
        self.average_rate_data_dir = params.average_rate_data_dir
        self.avg_rates_file = params.avg_rates_file
        self.event_data_dir = params.event_data_dir
        self.pca_data_dir = params.pca_data_dir
        self.temporal_stability_data_dir = params.temporal_stability_data_dir
        self.directories = params.all_directories
        self.figures = params.figures
        self.nmods = params.nmods
        self.freqs = params.freqs
        self.n_nmods = params.n_nmods
        self.n_freqs = params.n_freqs
        self.n_total_files = params.n_total_files
        self.n_trials_per_batch = params.n_trials_per_batch
        self.n_bins = params.n_bins
        self.id_range = params.id_range
        self.n_ids_sel = params.n_ids_sel
        self.dt = params.dt
        self.n_T = params.n_T
        self.n_Pop = params.n_Pop
        self.n_E = params.n_E
        self.n_I = params.n_I
        self.stim_time = params.stim_time                # stim len in seconds
        self.stim_n_ts = int(params.stim_time/params.dt) # no. of timestamps
        self.t_events = params.t_events
        self.sessions = params.sessions
        self.events = params.events
        self.n_batches = params.n_batches
        self.n_total_trials = params.n_total_trials
        self.students_t_50 = params.students_t_50
        self.students_t_500 = params.students_t_500
        self.n_p_outs = params.n_p_outs
        self.n_outputs = params.n_outputs
        self.p_outs = params.p_outs
        self.colors = params.colors
        self.colors_dark = params.colors_dark
        
        # make directories if they do not already exist
        self.ensure_directories_exist()
    
    def ensure_directories_exist(self):
        '''
        check that required directories exist, otherwise make them
        '''
        # make sure raw data directories exist
        if not os.path.isdir(self.raw_data_dir):
            warnings.warn('The raw data is missing! Expected ' +
                          self.raw_data_dir + ' to exist with subdirectories' +
                          ' for 10s and 8h recall conditions.')
        
        # make sure directories for proecessed data and figures exist, else make
        for directory in self.directories:
            if not os.path.isdir(directory):
                os.mkdir(directory)
                print('Created ' + directory)
       
    
    def get_file_names(self, str_list, data_dir):
        '''
        gets the list of file names from the data directory 
        
        Args:
            str_list: list of strings to search for in file names
            data_dir: directory to get files from
            
        Returns:
            sorted list of file names
        '''
        full_list = os.listdir(data_dir)
        final_list = [nm for ps in str_list for nm in full_list if ps in nm]
        return np.sort(final_list)
    
    
    def get_spks(self, full_path):
        '''
        loads spike times and ids given full path
        
        Args:
            full_path: path to the spike data
            
        Returns:
            spike times and ids as lists
        '''
        times = np.loadtxt(full_path, usecols=0)
        ids = np.loadtxt(full_path, usecols=1).astype(int)
        return times, ids
    
    
    def get_event_data(self, data_dir, session, n_batches, nm_i, f_j): 
        '''
        loads spike data for all trials in a given data directory and stores
        times and ids from each event (learning 1,2,3 and recall) for each trial
        
        Args:
            data_dir: path to raw spike data
            session: name of session (10s or 8h)
            n_batches: number of simulation batches (e.g. 5 batches of 10 simulations)
            nm_i: neuromodulator level index, int
            f_j: learning stimulation frequency index, int
            
        Returns:
            event_data: object containing spike times and ids and other info
        '''
        
        # get time at start
        t1 = time.time()
        
        # initialise event data object
        event_data = event_data_obj()
        
        # generate filename 
        filename = (self.event_data_dir + 
                    'session_' + session + 
                    '_nm_' + str(self.nmods[nm_i]) + 
                    '_freq_' + str(self.freqs[f_j]) + '.npy')
        
        # add filename to event data object
        event_data.which = filename
        
        # set file counter to zero, counts number of trials
        trial = 0
        
        # loop through each batch (i.e subdirectory)
        for batch in range(1, n_batches+1):
            
            # path for this nmod concentration and stimulation freq
            path = (data_dir + session + '/' + str(batch) + 
                    '/' 'nm=' + "%.2f" % self.nmods[nm_i] + 
                    ',f=' + str(self.freqs[f_j]) + 'Hz/')  
            # get and sort the spike raster file names from this directory
            files = self.get_file_names(['raster'], path)
            
            # loop through each of the files, get spike times and ids
            for file_i, file in enumerate(files):   
                # load spike times and ids for current file
                ts, ids = self.get_spks(path+file)
                
                # loop through events, i.e. learning 1,2,3 and recall
                for event_i, event in enumerate(self.events):
                    # extract event data
                    ts_sel, ids_sel = self.select_event_spks(session, 
                                                             event, 
                                                             ts, 
                                                             ids)
                    
                    # store selected times and ids 
                    event_data.ts[trial, event_i] = ts_sel
                    event_data.ids[trial, event_i] = ids_sel
                
                # store file and path
                event_data.files[trial] = file
                event_data.path[trial] = path
                
                # increase file (i.e. trial) counter
                trial += 1   
            
        # store number of files
        event_data.n_files = trial
        
        # get time at end to measure compute time
        t2 = time.time()
                    
        # store time required in event data object
        event_data.time_required = t2 - t1
        
        # return the event data object
        return event_data
        
    
    def select_event_spks(self, session, event, ts, ids):
        '''
        selects the spike times and ids belonging to a particular event
        
        Args:
            session: specifies session (10s or 8h)
            event: specifies event (learn1, learn2, learn3, recall)
            ts: time stamps of spikes, list
            ids: ids of neurons that spikes, list
        
        Returns:
            selected spike times and ids based on session/event 
        '''
        # condition 1: only times after t_event
        cond1 = ts > self.t_events[session][event]
        # condition 2: only times before t_event + stim_time
        cond2 = ts < (self.t_events[session][event] + self.stim_time)
        # condition 3: only ids in range
        cond3 = (ids >= self.id_range[0]) & (ids <= self.id_range[1])
        # select if all conditions met
        sel = cond1 & cond2 & cond3
        
        # select spike times and shift them by t_event
        ts_sel = ts[sel]
        # select ids
        ids_sel = ids[sel]
        
        return ts_sel, ids_sel
    
    
    def save_event_data_obj(self, event_data_obj):
        '''
        save event data object
        
        Args:
            event_data_obj: takes an event data object (see class)
        '''
        filename = event_data_obj.which
        with open(filename, "wb") as output_file:
            cPickle.dump(event_data_obj, output_file)
    
    
    def load_event_data_obj(self, session, nmod, freq):
        '''
        load event data object
    
        Args:
            session: specifies session (10s or 8h)
            nmod: specifies neurmodulator concentration
            freq: specifies learning stimulation frequency
        '''
        filename = (self.event_data_dir + 'session_' + str(session) +
                           '_nm_' + str(nmod) + '_freq_' + str(freq) + '.npy')
        # load and return event data
        with open(filename, "rb") as input_file:
            return cPickle.load(input_file)
    
    
    def report_progress(self, session, nm_i, f_j, now, start, last_time):
        '''
        computes and prints the progress and compute times
        
        Args:
            session: current session
            nm_i: current index of neuromodulator concentration 
            f_j: current index of learning stim. frequency
            now: current time
            start: time at beginning of computation
            last_time: time after the last iteration
        '''
        # update progress
        this_time = np.round(now - last_time, 2)
        total_time = np.round(now - start, 2)
        
        print('\r' + 'session: ' + session + ', ' 
              + str(np.round(100*(self.n_nmods*nm_i + f_j)
                                 /(len(self.nmods)*len(self.freqs)), 2))
              + '%, nm=' + str(self.nmods[nm_i]) + ', f=' + str(self.freqs[f_j]) 
              + ', t=' + str(this_time) 
              + ', T=' + str(total_time), end="")
    
    
    def _perform_event_data_extraction(self):
        '''
        loop through sessions, nmods, and freqs and get spike data for each event
        '''
        start = time.time()
        last_time = start
        for session_i, session in enumerate(self.sessions): 
            for nm_i, nm in enumerate(self.nmods):
                for f_j, f in enumerate(self.freqs):
                    # get event data for all trials, e.g. spikes, ids
                    event_data = self.get_event_data(self.spike_data_dir, 
                                                     session, 
                                                     self.n_batches, 
                                                     nm_i, 
                                                     f_j)
                    
                    # save event data object
                    self.save_event_data_obj(event_data)
                    
                    # get current time for progress bar
                    now = time.time()
                    
                    # report progress
                    self.report_progress(session, nm_i, f_j, now, start, last_time)
                    
                    # update time keeper
                    last_time = time.time()
    
    
    def extract_all_event_data(self):
        '''
        checks if event data is already stored. if yes, notifies where the data is.
        if not, loops through sessions, nmods, and freqs and gets spike data for each event.
        in this case helper method _perform_event_data_extraction() is called.
        '''
   
        # get all files in event directory, ignore hidden files
        all_files = [f for f in os.listdir(self.event_data_dir) if not f.startswith('.')]

        # count files
        n_event_files = len(all_files)

        # extract event data if need be
        if n_event_files == 0:
            self._perform_event_data_extraction()
        if n_event_files == self.n_total_files:
            print('Event data seems to already exist in the directory: /processed_data/event_data/.')
        elif (n_event_files != self.n_total_files) and (n_event_files != 0):
            warnings.warn("The number of files in .../event_data/ should match the number of parameter configurations (i.e. n_total_files).")
                    
                 
    def bin_data(self, ts, n_bins):
        '''
        min spike times
        
        Args:
            ts: list of spike times
            n_bins: total number of bins into which to bin
        
        Returns:
            binned spike times
        '''
        # compute bin width, 1000 ms * stimulation time / number of bins
        time_steps_per_bin = (self.stim_time/self.dt) // n_bins
        
        # bin the times and return as integer indices
        ts_binned = ((ts/self.dt)//time_steps_per_bin).astype(int)
        
        return ts_binned
                 
        
    def make_raster(self, ts, ids, size):
        '''
        make raster given spike times and ids
        
        Args:
            ts: spike times
            ids: neuron ids corresponding to spike times
            size: size of the raster (n_neurons, n_times)
        
        Return:
            spike raster of dimension size
        '''
        raster = np.zeros(size)
        raster[ids, ts] = 1
        return raster
    
    
    def make_figure_2a(self):
        '''
        make figure 2a: weight distributions after 8h for CA and outgoing
        '''
        
        # columns of file to read
        columns = [np.arange(1,33,8), np.arange(3,33,8)]
        
        # labels for neuromdulator and synapse groups
        labels = ['None', 'Low', 'Moderate', 'High']
        groups = ['CA', 'Outgoing']
        
        # number of labels and groups
        n_labels = len(labels)
        n_groups = len(groups)
        
        # get number rows
        n_rows = len(np.genfromtxt(self.avg_weights_file, 
                                   delimiter=';', 
                                   skip_header=2))
        
        # initialise data matrix, n_rows x 17 columns 
        weight_dist_matrix = np.zeros((n_rows, 33))

        # get data from file
        for row in range(n_rows):
            weight_dist_matrix[row,:] = np.genfromtxt(self.avg_weights_file, 
                                                      delimiter=';', 
                                                      dtype=float, 
                                                      skip_header=2)[row]
            
        # get bin labels for weights, y-axis
        weights = weight_dist_matrix[:,0]
        
        # rescale bin labels to percentages
        weights = 100*weights/weights[0]

        # initialise data matrices for weight values
        weights_not_pot = np.zeros((n_groups, n_labels))
        weights_not_pot_conf_int = np.zeros((n_groups, n_labels))
        
        # plot for cell assembly and outgoing synapses (i.e. groups)
        for group_i, group in enumerate(groups):
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(1,1,1)
            
            # go through each neuromodulator level
            for nm_i in range(len(labels)):
                weight_dist = weight_dist_matrix[:,columns[group_i][nm_i]]
                std = weight_dist_matrix[:,columns[group_i][nm_i]+1]

                conf_int_95 = self.students_t_50 * std / np.sqrt(self.n_total_trials)

                weights_low = weight_dist - conf_int_95
                weights_high = weight_dist + conf_int_95

                # Potentation under a couple percent considered "not potentiated"
                # STD via Gaussian Error Propagation
                # combine the first two bins (very low or no potentiation)
                
                # weights
                weights_not_pot[group_i, nm_i] = sum(weight_dist[0:2])
                
                # confidence interval 
                conf_int = np.sqrt((self.students_t_50*std[0]/np.sqrt(self.n_total_trials))**2 + 
                                   (self.students_t_50*std[1]/np.sqrt(self.n_total_trials))**2)
                weights_not_pot_conf_int[group_i, nm_i] = conf_int
                
                # PLOT
                ax.plot(weights, weight_dist, self.colors[nm_i][1], label=labels[nm_i])
                ax.fill_between(weights, weights_low, weights_high, 
                                facecolor=self.colors[nm_i][0], alpha=0.1)

            plt.title('Weight distribution after 8h, ' + group)
            plt.xlabel('Mean weight (%)')
            plt.ylabel('Fraction of cells')
            plt.xlim(104, 182)
            plt.ylim(0,0.575)
            if group == 'Outgoing':
                plt.legend(loc = 'upper right')

            fname = 'Fig2a_' + group
            plt.savefig(self.figures+fname, format='svg')
        plt.show()
        
        # weights not potentiated
        plt.figure(figsize=(6,4))
        plt.bar(labels, weights_not_pot[0,:], yerr=weights_not_pot_conf_int[0,:], 
                width=0.7, color=['black', '#E7A42C', '#4AA42C', '#2C69AC'])
        plt.title('Unpotentiated weights, CA')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        fname = 'Fig2a_CA_not_potentiated'
        plt.savefig(self.figures+fname, format='svg')
        plt.show()

        plt.figure(figsize=(6,4))
        plt.bar(labels, weights_not_pot[1,:], yerr=weights_not_pot_conf_int[1,:], 
                width=0.7, color=['black', '#E7A42C', '#4AA42C', '#2C69AC'])
        plt.title('Unpotentiated weights, Outgoing synapses')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        fname = 'Fig2a_Outgoing_not_potentiated'
        plt.savefig(self.figures+fname, format='svg')
        plt.show()
        
    
    def make_figure_2c(self):
        '''
        make figure 2c: example spike rasters
        '''
        
        # parameters specific to this figure
        batch = 1    
        file_ex = 4
        session = self.sessions[1]
        freq = self.freqs[5]                 # 60 Hz
        nms = [self.nmods[3], self.nmods[9]]  # low: 0.06, high: 0.18
        labels = ['low', 'high']
        
        for nm_i, nm in enumerate(nms):
            path = (self.spike_data_dir + session + '/' + str(batch) + '/' +
                    'nm=' + "%.2f" % nm + ',f=' + str(freq) + 'Hz/')

            # load and sort the spike raster files from this directory
            file = self.get_file_names(['raster'], path)[file_ex]

            # load times, list
            ts = np.loadtxt(path+file, usecols=0)

            # load neuron indices corresponding to spike times
            ids = np.loadtxt(path+file, usecols=1).astype(int)

            # separate into assembly, control neurons, and inhibitory
            times_a = ts[ids< 150]
            inds_a = ids[ids< 150]

            times_c = ts[(ids >= 150) & (ids < 1600)]
            inds_c = ids[(ids >= 150) & (ids < 1600)]

            times_i = ts[ids > 1600]
            inds_i = ids[ids > 1600]

            # full raster plot
            plt.figure(figsize=(10,4))
            plt.plot(times_a, inds_a, color=self.colors[2*nm_i+1][1], 
                     marker='.', linestyle = 'None', markersize=10)
            plt.plot(times_c, inds_c, markerfacecolor='lightgrey', 
                     markeredgecolor=self.colors[2*nm_i+1][1], marker='.', 
                     linestyle = 'None', markersize=10)
            plt.xlabel('Time (ms)')
            plt.ylabel('Neuron')
            plt.xticks(np.linspace(28810, 28810.11, 12), 
                       np.linspace(0, 110, 12, dtype=int))
            plt.xlim((28809.999, 28810.113))
            
            # save and show plot
            fname = 'Fig2c_raster_' + labels[nm_i]
            plt.savefig(self.figures+fname, format='svg')
            plt.show()

            # zoomed in raster plot
            plt.figure(figsize=(10,4))
            plt.plot(times_a, inds_a, color=self.colors[2*nm_i+1][1], 
                     marker='.', linestyle = 'None', markersize=10)
            plt.plot(times_c, inds_c, markerfacecolor='lightgrey', 
                     markeredgecolor=self.colors[2*nm_i+1][1], marker='.', 
                     linestyle = 'None', markersize=10)
            plt.xlabel('Time (ms)')
            plt.ylabel('Neuron')
            plt.xticks(np.linspace(28810.070, 28810.075, 6), 
                       np.linspace(70, 75, 6, dtype=int))
            plt.xlim((28810.070, 28810.075))

            # save and show plot
            fname = 'Fig2c_raster_' + labels[nm_i] + '_zoom_70-75ms'
            plt.savefig(self.figures+fname, format='svg')
            plt.show()
        
    
    def make_figure_2d(self):
        '''
        make figure 2d: firing rate distribution, non-core
        '''

        n_rows = len(np.genfromtxt(self.avg_rates_file, delimiter=';', skip_header=2))
        firing_rate_matrix = np.zeros((n_rows, 17))

        for row in range(n_rows):
            firing_rate_matrix[row,:] = np.genfromtxt(self.avg_rates_file, delimiter=';', 
                                                      dtype=float, 
                                                      skip_header=2)[row]
        rates = firing_rate_matrix[:,0]

        columns = [[3, 7, 11, 15], [1, 5, 9, 13]]
        labels = ['None', 'Low', 'Moderate', 'High']
        groups = ['Non-CA', 'CA']
        
        # non-CA
        group_i = 0
        group = groups[group_i]
        
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1, 1, 1)
        for nm_i in range(4):
            rate_dist = firing_rate_matrix[:,columns[group_i][nm_i]]
            std = firing_rate_matrix[:,columns[group_i][nm_i]+1]

            conf_int_95 = self.students_t_50 * std / np.sqrt(self.n_total_trials)

            to_plot = rate_dist>0
            rate_dist = rate_dist[to_plot]
            conf_int_95 = conf_int_95[to_plot]

            rate_low = rate_dist - conf_int_95
            rate_high = rate_dist + conf_int_95

            ax.plot(rates[to_plot], rate_dist, self.colors[nm_i][1], 
                    label=labels[nm_i], marker='.')
            ax.fill_between(rates[to_plot], rate_low, rate_high, 
                            facecolor=self.colors[nm_i][0], alpha=0.1)

        ax.set_yscale('log')
        plt.xlabel('Firing rate (Hz)')
        plt.ylabel('Fraction of cells')
        plt.title('Activity of non-core neurons')
        plt.legend()

        fname = 'Fig2d_' + group
        plt.savefig(self.figures+fname, format='svg')
        plt.show()

           
    def make_figure_3h(self):
        '''
        make fig 3h: sorted firing rate distributions
        ''' 
        
        colors = self.colors
        figures = self.figures
        freq = self.freqs[5]                  # 60 Hz
        nms = [self.nmods[3], self.nmods[9]]  # low: 0.06, high: 0.18
        sessions = self.sessions
        labels = ['low', 'high']
        all_events = self.events
        events_to_plot = ['learn3', 'recall']
        n_events = len(events_to_plot)
        bins = np.arange(0,self.n_E+1)
        
        # loop through both conditions, 10s and 8h recall
        for session_i, session in enumerate(sessions):
            
            # load event data for low and high neuromodulator conditions
            event_data_low = self.load_event_data_obj(session, nms[0], freq)
            event_data_high = self.load_event_data_obj(session, nms[1], freq)
            
            # create arrays to store firing rate distributions for low and high 
            # neuromodulator levels
            low = np.zeros((n_events, self.n_E, int(self.n_total_trials)))
            high = np.zeros((n_events, self.n_E, int(self.n_total_trials)))
            
            # for learn3 and recall events, for each trial 
            # compute firing rate distributions
            for k, event in enumerate(events_to_plot):
                # get index of event 
                event_index = all_events.index(event)
                for trial_i in range(self.n_total_trials):
                    
                    # get neuron ids for low and high nmod
                    ids_l = event_data_low.ids[trial_i, event_index]
                    ids_h = event_data_high.ids[trial_i, event_index]
                    
                    # compute firing rate distributions, low and high nmod
                    low[k, :, trial_i] = (np.histogram(ids_l, bins=bins)[0] / 
                                          self.stim_time)
                    high[k, :, trial_i] = (np.histogram(ids_h, bins=bins)[0] / 
                                           self.stim_time)
        
            # SORT BY FIRING RATE
            # arrays to store sorted firing rate distributions 
            low_sorted = np.zeros((n_events, self.n_E, int(self.n_total_trials)))
            high_sorted = np.zeros((n_events, self.n_E, int(self.n_total_trials)))

            # sort by firing rate independently for stimulated assembly 
            # neurons (a), non-stimulated assembly neurons (na), and control 
            # neurons (c)
            for trial_i in range(int(self.n_total_trials)):
                # get sorted indices for low nmod
                inds_low_a = np.argsort(low[0, 0:75, trial_i])[::-1]
                inds_low_na = 75 + np.argsort(low[0, 75:150, trial_i])[::-1]
                inds_low_c = 150 + np.argsort(low[0, 150:1600,trial_i])[::-1]
                inds_low = np.concatenate([inds_low_a, 
                                           inds_low_na, 
                                           inds_low_c]).astype(int)

                # get sorted indices for high nmod
                inds_high_a = np.argsort(high[0, 0:75, trial_i])[::-1]
                inds_high_na = 75 + np.argsort(high[0, 75:150, trial_i])[::-1]
                inds_high_c = 150 + np.argsort(high[0, 150:1600,trial_i])[::-1]
                inds_high = np.concatenate([inds_high_a, 
                                            inds_high_na, 
                                            inds_high_c]).astype(int)
                
                # for each event sort by the new index order
                for event_i in range(n_events):
                    low_sorted[event_i, :, trial_i] = low[event_i, 
                                                          inds_low, 
                                                          trial_i]
                    high_sorted[event_i, :, trial_i] = high[event_i, 
                                                            inds_high, 
                                                            trial_i]

            # MOVING AVERAGE FOR RECALL
            window = 20
            low_sorted_mean = np.zeros(self.n_E)
            high_sorted_mean = np.zeros(self.n_E)
            # go through each neuron, if too close to boundary, 
            # only take next N - i values
            for i in range(self.n_E):
                if np.abs(75-i)<window & ((75-i)>0):
                    start = i
                    stop = 75
                elif (np.abs(150-i)<window) & ((150-i)>0):
                    start = i
                    stop = 150
                else:
                    start = i
                    stop = i+window

                # compute average for moving average
                low_sorted_mean[i] = np.mean(low_sorted[1,start:stop,:])
                high_sorted_mean[i] = np.mean(high_sorted[1,start:stop,:])


            # list of neuron ids
            N = np.arange(self.n_E)

            # generate the figure
            fig = plt.figure(figsize=(10,4))
            ax = fig.add_subplot(1, 1, 1)
            
            # Plot data points for high and low neuromodulator, at recall
            ax.plot(N, np.mean(high_sorted[1,:,:], axis=1).T, '.', 
                    color=colors[3][1], markersize=6, alpha=0.3)
            ax.plot(N, np.mean(low_sorted[1,:,:], axis=1).T, '.', 
                    color=colors[1][1], markersize=6, alpha=0.2) 
            
            # Plot average firing rate distribution from learning as reference
            ax.plot(N, np.mean(high_sorted[0,:,:], axis=1).T, 
                    label='Learning, high NM', linewidth=3, color=colors[3][1])
            ax.plot(N, np.mean(low_sorted[0,:,:], axis=1).T, 
                    label='Learning, low NM', linewidth=2, color=colors[1][1])

            # Plot high neuromodulator values, moving average
            ax.plot(N[0:75], high_sorted_mean[0:75], '--', label='Recall, high NM', 
                    linewidth=3, color=colors[3][1], alpha=0.9)
            ax.plot(N[75:150], high_sorted_mean[75:150], '--', 
                    linewidth=3, color=colors[3][1], alpha=0.9)
            ax.plot(N[150:], high_sorted_mean[150:], '--', 
                    linewidth=3, color=colors[3][1], alpha=0.9)

            # Plot low neuromodulator values, moving average
            ax.plot(N[0:75], low_sorted_mean[0:75], '--', label='Recall, low NM', 
                    linewidth=2, color=colors[1][1], alpha=0.9)
            ax.plot(N[75:150], low_sorted_mean[75:150], '--', 
                    linewidth=2, color=colors[1][1], alpha=0.9)
            ax.plot(N[150:], low_sorted_mean[150:], '--', 
                    linewidth=2, color=colors[1][1], alpha=0.9)

            # plot specs
            plt.title('During ' + str(session) + '-recall')
            plt.xlim(0,1600)
            plt.xlabel('Neuron (sorted)')
            plt.ylabel('Firing rate')
            ax.set_yscale('log')
            plt.legend(handlelength=3)
            
            # save and show figure
            fname = 'Fig3h_' + session
            plt.savefig(figures+fname, format='svg')
            plt.show()
        
    
    def compute_pca(self, event, n_bins_for_spikes):
        '''
        compute PCA for each binned raster and check how many dimensions 
        needed for var_thr percent variance
        
        Args:
            event: for which event to compute pca
            n_bins_for_spikes: total number of bins to bin spike times
        
        Return:
            the average dimensionality over trials for each session, 
            neuromodulator concentration, and learning stimulation freq.
        '''
        
        # for timing, progress reporting
        start = time.time()
        last_time = start
        
        # data directory
        pca_data_dir = self.pca_data_dir
        
        # parameters
        event_i = self.events.index(event)
        size = (self.n_E, n_bins_for_spikes)
        var_thr = 0.7
        
        # set random seed
        np.random.seed(0)
        
        # array to store dimensionality
        dim = np.zeros((len(self.sessions), self.n_nmods, self.n_freqs, 
                        int(self.n_total_trials)))
        
        # go through each session, nmod, freq
        for session_i, session in enumerate(self.sessions):
            for nm_i, nm in enumerate(self.nmods):
                for f_j, f in enumerate(self.freqs):
                    
                    # load event data
                    event_data = self.load_event_data_obj(session, nm, f)
                    
                    # for each trial, compute pca and dimensionality
                    for trial in range(int(self.n_total_trials)):
                        
                        # time stamps and ids to be binned
                        ts = (event_data.ts[trial, event_i] - 
                              self.t_events[session][event])
                        ids = event_data.ids[trial, event_i]
                        
                        # bin data
                        ts_binned = self.bin_data(ts, n_bins_for_spikes)
                        
                        # make raster
                        raster = self.make_raster(ts_binned, ids, size)
                        
                        # do PCA
                        pca = PCA(n_components=n_bins_for_spikes)
                        pca.fit(raster.T)

                        # check dimensionality
                        k = 0
                        cumul_var = np.sum(pca.explained_variance_ratio_[:k])
                        while (k <= n_bins_for_spikes) & (cumul_var <= var_thr):
                            k += 1
                            cumul_var = np.sum(pca.explained_variance_ratio_[:k])
                            
                        dim[session_i, nm_i, f_j, trial] = k
                
                    # get current time for progress bar
                    now = time.time()

                    # report progress
                    self.report_progress(session, nm_i, f_j, now, start, last_time)

                    # update time keeper
                    last_time = time.time()

                
        # take average and std over trials
        dim_mean = np.mean(dim, axis=3)
        dim_ste = np.std(dim, axis=3)/np.sqrt(int(self.n_total_trials))

        # save data
        np.save(pca_data_dir + 'pca_dim_mean_' + str(n_bins_for_spikes) + '.npy' , dim_mean)
        np.save(pca_data_dir + 'pca_dim_ste_' + str(n_bins_for_spikes) + '.npy' , dim_ste)
        
        # return the average dimensionality
        return dim_mean
    
    def make_figure_5(self):
        '''
        makes figure 5: dimensionality via pca
        '''
        # parameters
        figures = self.figures
        processed_data = self.processed_data_dir
        pca_data_dir = self.pca_data_dir
        event = 'recall'
        size = (self.n_E, self.n_bins)
        
        file = 'pca_dim_mean_' + str(self.n_bins) + '.npy' 
        files = [f for f in os.listdir(pca_data_dir) if not f.startswith('.')]
        
        # if data already stored, load it, else compute pcas and store it
        if file in files:
            dim_mean = np.load(pca_data_dir + file)
        else:
            # compute the pca for all event data
            dim_mean = self.compute_pca(event, self.n_bins)

        # plot means
        for session_i, session in enumerate(self.sessions):
            plt.figure()
            plt.imshow(dim_mean[session_i,:,:].T, origin='lower', cmap='viridis', 
                       vmin=np.min(dim_mean), vmax=np.max(dim_mean))
            plt.xlabel('Neuromodulatior Concentration', fontsize=18)
            plt.ylabel('Stimulation Frequency', fontsize=18)
            plt.xticks(np.arange(10), self.nmods, fontsize=11)
            plt.yticks(np.arange(10), self.freqs, fontsize=11)
            plt.colorbar()
            fname = 'Fig5_' + session
            plt.savefig(figures+fname, format='svg')
            plt.show()
            
        plt.figure()
        plt.imshow(dim_mean[1,:,:].T - dim_mean[0,:,:].T, origin='lower', 
                   cmap='seismic', vmin=-10, vmax=10)
        plt.xlabel('Neuromodulatior Concentration', fontsize=18)
        plt.ylabel('Stimulation Frequency', fontsize=18)
        plt.xticks(np.arange(10), self.nmods, fontsize=11)
        plt.yticks(np.arange(10), self.freqs, fontsize=11)
        plt.colorbar()
        
        fname = 'Fig5_diff'
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
    def compute_stability(self, nms, freq, filenames):
        '''
        compute spike time stability for figure 6
        
        Args:
            nms: list of neuromodulator concentrations
            freq: learning stimulation frequency
            filenames: list of file names for storing data
        
        Return:
            spike time stability for real data and control (shuffled)
        '''
        
        # vars
        session = '8h'
        spk_thr = 0
        size = (self.n_E, self.n_bins)
        
        # initialise arrays
        stability = np.zeros((len(nms), self.n_total_trials, self.n_E))
        stability_ctrl = np.zeros((len(nms), self.n_total_trials, self.n_E))
        n_spks_learning = np.zeros((len(nms), self.n_total_trials, self.n_E))
        
        # set random seed 
        np.random.seed(0)
        
        # loop through nm concentrations
        for nm_i, nm in enumerate(nms):
                
            # load event data
            event_data = self.load_event_data_obj(session, nm, freq)

            for trial in range(self.n_total_trials):
                for n in range(self.n_E):
                    # note that this method of computing stability doesn't 
                    # work if more than one spike occurs per bin. this would 
                    # need to be taken into account. here bins are 1ms so due
                    # to refractory period, cannot have more than one spike.

                    # get time stamps for neuron n, third learning stimulus and recall
                    ts_learn3 = event_data.ts[trial, 2][event_data.ids[trial,2] == n] - self.t_events['8h']['learn3']
                    ts_recall = event_data.ts[trial, 3][event_data.ids[trial,3] == n] - self.t_events['8h']['recall']

                    # bin data
                    learn3 = self.bin_data(ts_learn3, self.n_bins)
                    recall = self.bin_data(ts_recall, self.n_bins)

                    # number of spikes at learning
                    n_spks_learning[nm_i, trial, n] = len(learn3)

                    # number of spikes that match between learning and recall
                    n_same = len(np.intersect1d(learn3, recall))

                    # CONTROL: shuffled data
                    # make shuffled indices for control condition
                    new_inds_learn3 = np.random.choice(self.n_bins, self.n_bins, replace=False)
                    new_inds_recall = np.random.choice(self.n_bins, self.n_bins, replace=False)

                    # shuffle spike order
                    learn3_shuffled = np.where(np.isin(new_inds_learn3, learn3))[0]
                    recall_shuffled = np.where(np.isin(new_inds_recall, recall))[0]

                    # number of (shuffled!) spikes that match between learning and recall
                    n_same_shuffled = len(np.intersect1d(learn3_shuffled, recall_shuffled))

                    # compute stability for real and control
                    if n_spks_learning[nm_i, trial, n] > spk_thr:
                        stability[nm_i, trial, n] = n_same / n_spks_learning[nm_i, trial, n]
                        stability_ctrl[nm_i, trial, n] = n_same_shuffled / n_spks_learning[nm_i, trial, n]
                    # if not enough spikes, put nan values
                    else:
                        stability[nm_i, trial, n] = float('nan')
                        stability_ctrl[nm_i, trial, n] = float('nan')
                        
            # print progress
            print('Progress: ' + str(nm_i+1) + '/' + str(len(nms)))
                    
        # save data
        np.save(self.temporal_stability_data_dir + filenames[0], stability)
        np.save(self.temporal_stability_data_dir + filenames[1], stability_ctrl)
        np.save(self.temporal_stability_data_dir + filenames[2], n_spks_learning)
                    
        return stability, stability_ctrl
    
    def _make_figure_6b(self):
        '''
        make figure 6b: spike time stability as a function of real vs. shuffled
        '''
        # directories 
        temporal_stability_data_dir = self.temporal_stability_data_dir
        figures = self.figures
        
        # define colors for plot locally, modify black to a dark grey
        colors = {0: ['black', '#434445'],
                  1: ['yellow', '#E7A42C'],
                  2: ['green', '#4AA42C'],
                  3: ['blue', '#2C69AC'],
                  }
        
        # neuromodulator concentrations and learning stim frequency
        nms = [0.0, 0.06, 0.12, 0.18]
        freq = 60

        # file names for stability data
        file_stab = 'stability.npy' 
        file_ctrl = 'stability_ctrl.npy' 
        file_nspks = 'n_spks_learn3.npy'
        filenames = [file_stab, file_ctrl, file_nspks]
        
        # get files in directory
        files = [f for f in os.listdir(temporal_stability_data_dir) if not f.startswith('.')]
        
        # if data already stored, load it, else compute pcas and store it
        if (file_stab in files) and (file_ctrl in files) and (file_nspks in files):
            stability = np.load(temporal_stability_data_dir + file_stab)
            stability_ctrl = np.load(temporal_stability_data_dir + file_ctrl)
            print('Data was loaded from ' + temporal_stability_data_dir)
        else:
            # compute the temporal stability
            stability, stability_ctrl = self.compute_stability(nms, freq, filenames)
            print('Temporal stability was computed and data stored in ' + temporal_stability_data_dir)
        
        # MAKE PLOT
        # select subset of neurons, 0-1600 for all
        n1 = 0
        n2 = 1600
        
        # storage for stability data
        stab = np.zeros((2*len(nms), 50))
        count = 0
        for i in range(len(nms)):
            stab[count, :] = np.nanmean(stability[i, :, n1:n2], axis=1)
            count += 1
            stab[count, :] = np.nanmean(stability_ctrl[i, :, n1:n2], axis=1)
            count += 1
            
            
        fname = 'Fig6b_stability_violinplot'
        fig, ax = plt.subplots(figsize=(10,5))
        # make the violin plot without median, quartiles and whiskers
        bp1 = seaborn.violinplot(data=stab.T, 
                                 palette=[colors[i][1] for i in range(len(nms))
                                          for j in range(2)], 
                                 scale='width', 
                                 cut=2,
                                 inner=None,
                                 saturation=0.8)
        
        bp1.set_xlabel('Neuromodulator level', fontsize=16)
        bp1.set_ylabel('Stability', fontsize=16)
        bp1.set_xticklabels(['real', 'shuffled'] * len(nms), fontsize=14)
        
        # make the boxplot inset and add median as a white point
        for i in range(8):
            ax.boxplot(stab[i,:], whis=[0,100], positions=np.array([i]),
                       showcaps=False, widths=0.05, patch_artist=True,
                       boxprops=dict(color='black', facecolor='black'),
                       whiskerprops=dict(color='black', linewidth=2),
                       medianprops=dict(linewidth=0), zorder=5)
            ax.scatter(i, np.median(stab[i,:]), c='white', s=1, zorder=10)
        
        plt.savefig(figures + fname, format='svg')
        plt.show()
        
    def extract_outgoing_weights(self, nms, freq):
        '''
        extract the outgoing late phase weights from the core to each other neuron
        
        Args:
            nms: list of neuromodulator concentrations
            freq: learning stimulation frequency
        
        Return:
            array containing the outgoing late phase weights from each core neuron
            for each nm, trial, and neuron
        '''
        data_dir = self.outgoing_weight_data_dir
        
        mean_LPW_core_all_to_n = 3  # index for late phase mean weights
        data = np.zeros((16000, 12))
        outgoing_mean = np.zeros(16000)
        weight = np.zeros((len(nms),self.n_total_trials, self.n_E))

        for nm_i, nm in enumerate(nms):
            # for nm=0 all late phase weights are zero, they are not stored, just populated here
            if nm == 0:
                weight[nm_i, :, :] = np.zeros((self.n_total_trials, self.n_E))
            else:
                for batch in range(self.n_batches):
                    weight_file = str(batch+1)+'/nm='+str(nm)+',f='+str(freq)+'Hz/Params_OutgoingWeights.txt'

                    f = open(data_dir+weight_file, 'r')
                    raw = f.read()
                    line_split = re.split('\n', raw)
                    for i in range(16000):
                        row_split = re.split('\t\t', line_split[i])
                        data[i, :] = row_split[int(i % 1600 == 0):]  # this line gets rid of time stamp at every 1600th entry
                        outgoing_mean[i] = data[i, mean_LPW_core_all_to_n]

                    weight[nm_i, batch*self.n_trials_per_batch:(batch*self.n_trials_per_batch + self.n_trials_per_batch), :] = outgoing_mean.reshape(10,1600)
        return weight
        
    def _make_figure_6c(self):
        '''
        make figure 6c: scatter plots stability vs late phase outgoing weight
        for each level of neuromodulator
        '''
        
        figures = self.figures
        
        colors = self.colors
        colors_dark = self.colors_dark
        
        nms = [0.0, 0.06, 0.12, 0.18]
        freq = 60
        
        # get weights
        weight = self.extract_outgoing_weights(nms, freq)
        
        stability = np.load(self.temporal_stability_data_dir + 'stability.npy')
        num_spks = np.load(self.temporal_stability_data_dir + 'n_spks_learn3.npy')
        
        n1 = 150
        n2 = 1600
        tr0 = 0
        tr1 = 50

        n_thr_min = 0 

        line_lens = [[0,      0     ],
                     [0.0001, 0.0038],
                     [0.0005, 0.0070],
                     [0.0005, 0.0110]]
        
        fnames = ['Fig6c', 'Fig6d', 'Fig6e', 'Fig6f']

        for nm_i, nm in enumerate(nms):
            
            w = weight[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))
            n = num_spks[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))
            s = stability[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))

            not_nan_inds = ~np.isnan(w) & ~np.isnan(s)
            n_cutoff_inds = (n > n_thr_min) #& (n < n_thr_max)

            w = w[not_nan_inds & n_cutoff_inds] # & pos_s_inds]
            n = n[not_nan_inds & n_cutoff_inds]
            s = s[not_nan_inds & n_cutoff_inds] # & pos_s_inds]

            plt.figure(figsize=(5,5))
            if nm == 0:
                plt.scatter(w, s, s=10, c=colors[nm_i][1], alpha=0.4)
            else:
                plt.scatter(w, s, s=2, c=colors[nm_i][1], alpha=0.4)

            if nm != 0:
                z = np.polyfit(w, s, 1)
                p = np.poly1d(z)

                w_sel = w[(w >= line_lens[nm_i][0]) & (w <= line_lens[nm_i][1])]
                plt.plot(w_sel, p(w_sel), color=colors_dark[nm_i][1], linewidth=4)
                print('corr: ', stats.spearmanr(w,s))
                print('n: ', len(w))
            
            plt.xlabel('Mean Late Phase Weight from Core Neurons to Non-Core Neuron n')
            plt.ylabel('Stability of spike times')
            plt.xlim(-0.00001,0.014)
            plt.ylim(-0.015, 1.05)

            plt.savefig(figures + fnames[nm_i] + '_stability_vs_outgoing_weights_nm_' + str(nm), format='svg')
            plt.show()
        
        
    def make_figure_6(self):
        '''
        makes figure 6: see self.make_figure_6b() and self.make_figure_6c()
        '''
        self._make_figure_6b()
        self._make_figure_6c()
        
        
    def compute_regression(self, neuron_start, neuron_stop, name):
        '''
        Computes the ridge regression using the binned spike data on random walk
        target functions for different percent output connectivities across 
        neuromodulator concentrations and stimulation frequencies
        
        Used to produce Figure 7 and related supplementary figures
        
        Args:
            neuron_start: defines the neuron group to use for regression
            neuron_stop: defines the neuron group to use for regression
            name: a suffix for the files that get saved (e.g. core assembly, 'ca')
            
        Returns:
            Filenames for where data is stored
        '''
        
        # directories
        regression_data_dir = self.regression_data_dir
        
        # file names for regression data
        file_r_sq = 'r_sq_r_learn_' + name + '.npy'
        file_y_hat = 'y_hat_store_' + name + '.npy'
        file_prmse = 'prmse_r_learn_' + name + '.npy'
        file_y = 'y_store_' + name + '.npy'
        filenames = [file_r_sq, file_y_hat, file_prmse, file_y]
        
        # get which files exist in directory
        files = [f for f in os.listdir(regression_data_dir) if not f.startswith('.')]
        
        # if data already stored we can load it, else compute pcas and store it
        if all(file in files for file in filenames):
            print('Data exists and will be loaded from ' + regression_data_dir)
            return file_r_sq, file_y, file_y_hat
        else:
            # need to compute regression
            #self.compute_regression(name)
            print('Regression will be computed and data stored in ' + regression_data_dir)
        
        # initialize rige regression 
        ridgereg = Ridge(alpha=0.1) 

        # function to compute percent root mean square error
        def prmse(estimation, target):
            return 100*np.sqrt(np.mean((estimation - target)**2)/np.sum(target**2))

        # array of all neurons that will be used for regression
        neurons = np.arange(neuron_start, neuron_stop)

        # number of sessions (i.e. two, 10s and 8h)
        n_sessions = len(self.sessions)

        # number of target functions and percent output connectivity
        n_outputs = self.n_outputs
        n_p_outs = self.n_p_outs
        p_outs = self.p_outs

        # size of raster
        size = (self.n_E, self.n_bins)

        # arrays to store R square, error, predicted output and actual target function
        r_sq_r_learn = np.zeros((n_sessions, self.n_nmods, self.n_freqs, 
                                 int(self.n_total_trials), n_outputs, n_p_outs))
        prmse_r_learn = np.zeros((n_sessions, self.n_nmods, self.n_freqs, 
                                  int(self.n_total_trials), n_outputs, n_p_outs))
        y_hat_store = np.zeros((n_sessions, self.n_nmods, self.n_freqs, 
                                int(self.n_total_trials), n_outputs, n_p_outs, 
                                self.n_bins))
        y_store = np.zeros((n_outputs, self.n_bins))

        # for keeping track of compute time / progress bar
        start = time.time()
        
        # set seed for reproducibility
        np.random.seed(0) 
        for output in range(n_outputs):
            # generate a random walk
            y = np.zeros(self.n_bins)
            y[0] = 1
            std = 1
            for i in range(self.n_bins-1):
                y[i+1] = y[i] + std*np.random.normal(0)

            # store target function
            y_store[output, :] = y

            # tile target function, three learning stimuli + recall = 4
            Y = np.tile(y, 4)

            # loop over all percent output connectivities
            for p_i, p_out in enumerate(p_outs): 
                
                # compute number of neurons given percent output
                n_neurons = int(p_out*len(neurons))
                
                # select random neurons
                ids_sel = np.random.choice(neurons, n_neurons, replace=False)

                for session_i, session in enumerate(self.sessions):
                    for nm_i, nm in enumerate(self.nmods):
                        for f_j, f in enumerate(self.freqs):

                            # load event data
                            event_data = self.load_event_data_obj(session, nm, f)
                            raster = np.zeros((4, self.n_E, self.n_bins))

                            # for each trial, compute pca and dimensionality
                            for trial in range(int(self.n_total_trials)):

                                for i, e in enumerate(self.events):

                                    # time stamps and ids to be binned
                                    ts = (event_data.ts[trial, i] - 
                                          self.t_events[session][e])
                                    ids = event_data.ids[trial, i]

                                    # bin data
                                    ts_binned = self.bin_data(ts, self.n_bins)

                                    # make raster
                                    raster[i, :, :] = self.make_raster(ts_binned, 
                                                                       ids, size)

                                # training data
                                data =  np.concatenate((raster[0, ids_sel, :].T,
                                                        raster[1, ids_sel, :].T,
                                                        raster[2, ids_sel, :].T,
                                                        raster[3, ids_sel, :].T))
                                # test data
                                test = raster[3, ids_sel, :].T   

                                # compute ridge regression
                                model_r = ridgereg.fit(data, Y)
                                y_hat_r = model_r.predict(test)

                                # store data
                                r_sq_r_learn[session_i, nm_i, f_j, trial, output, p_i] = model_r.score(test, Y[:self.n_bins])
                                y_hat_store[session_i, nm_i, f_j, trial, output, p_i, :] = y_hat_r
                                prmse_r_learn[session_i, nm_i, f_j, trial, output, p_i] = prmse(Y[:self.n_bins], y_hat_r)
                
                # get current time for progress bar
                now = time.time()
                
                # print progress
                print('\r' + 'Target function: ' + str(output + 1) + '/' + str(n_outputs) + '; ' +
                      'Percent connectivity: ' + str(p_i + 1) + '/' + str(n_p_outs) + '; ' + 
                      'Elapsed time: ' + str(np.round(now - start, 2)), end='')

        np.save(regression_data_dir + 'r_sq_r_learn_' + name, r_sq_r_learn)
        np.save(regression_data_dir + 'y_hat_store_' + name, y_hat_store)
        np.save(regression_data_dir + 'prmse_r_learn_' + name, prmse_r_learn)
        np.save(regression_data_dir + 'y_store_' + name, y_store)
        
        return file_r_sq, file_y, file_y_hat
            
    def _get_data_for_regression_rasters(self, file_r_sq, name):
        '''
        Makes the raster of learning stim frequency vs. neuromodulator concentration
        for 10s and 8h recall, as well as the difference between them.
        
        Used for Fig 7abc and related supplemenary figure for core assembly only
        and non-core (i.e. all other neurons) only
        
        Args:
            file_r_sq: the file name where the R squared data is stored
            name: suffix for files and figures that will get stored (e.g. 'all', 'ca', or 'nc')
        
        Returns:
            Regression results for 10s, 8h and the difference
        '''
        
        # directories
        regression_data_dir = self.regression_data_dir
        
        # load data
        r_sq_r_learn = np.load(regression_data_dir + file_r_sq) 
        
        # make raster -----------------------------------------------------------------
        r_sq_raster = np.zeros((2, self.n_p_outs, self.n_nmods, self.n_freqs))
        r_sq_raster_ste = np.zeros((2, self.n_p_outs, self.n_nmods, self.n_freqs))
        for session_i, session in enumerate(self.sessions):
            for nm_i, nm in enumerate(self.nmods):
                for f_j, f in enumerate(self.freqs):
                    for p_k, p_out in enumerate(self.p_outs):
                        r_sq_raster[session_i, p_k, nm_i, f_j] = np.mean(r_sq_r_learn[session_i, nm_i, f_j, :, :, p_k], axis=(0,1))
                        r_sq_raster_ste[session_i, p_k, nm_i, f_j] = np.std(r_sq_r_learn[session_i, nm_i, f_j, :, :, p_k], axis=(0,1))/(self.n_nmods*self.n_freqs)

        r_sq_10s = np.mean(r_sq_raster[0,:,:,:], axis=0)
        r_sq_8h = np.mean(r_sq_raster[1,:,:,:], axis=0)
        r_sq_diff = r_sq_8h - r_sq_10s

        np.save(regression_data_dir + 'r_sq_10s_' + name, r_sq_10s)
        np.save(regression_data_dir + 'r_sq_8h_' + name, r_sq_8h)
        np.save(regression_data_dir + 'r_sq_diff_' + name, r_sq_diff)
        
        return r_sq_10s, r_sq_8h, r_sq_diff
        
        
    def _make_regression_rasters_plot(self, file_r_sq, fig, name):
        '''
        makes the rasters of R squared values for the regression results 
        for figure 7 and related supplementary figures
        
        Args:
            file_r_sq: path to data file
            fig: name of figure (e.g. Fig7, SuppFig5)
            name: suffix for processed data files to be saved
        '''
        
        # directories
        figures = self.figures
        
        # first make the rasters and store the data
        r_sq_10s, r_sq_8h, r_sq_diff = self._get_data_for_regression_rasters(file_r_sq, name)
        
        # plot and store figures
        fname = fig + 'a'
        plt.figure()
        plt.imshow(r_sq_10s.T, origin='lower')
        plt.xlabel('Neuromodulatior Concentration', fontsize=18)
        plt.ylabel('Stimulation Frequency', fontsize=18)
        plt.colorbar()
        plt.savefig(figures+fname, format='svg')
        plt.show()

        fname = fig + 'b'
        plt.figure()
        plt.imshow(r_sq_8h.T, origin='lower')
        plt.xlabel('Neuromodulatior Concentration', fontsize=18)
        plt.ylabel('Stimulation Frequency', fontsize=18)
        plt.colorbar()
        plt.savefig(figures+fname, format='svg')
        plt.show()

        fname = fig + 'c'
        plt.figure()
        plt.imshow(r_sq_diff.T, origin='lower', cmap='seismic', vmin=-.3, vmax=.3)
        plt.xlabel('Neuromodulatior Concentration', fontsize=18)
        plt.ylabel('Stimulation Frequency', fontsize=18)
        plt.colorbar()
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
        
    def _get_data_for_r_squared_vs_perc_conn_plot(self, file_r_sq):
        '''
        gets r squared data from file and computes means and ste
        
        Args:
            file_r_sq: path to data file
        
        Returns:
            mean and ste of r squared data
        '''
        
        # compute means and standard devations for example: FREQ = 60Hz 
        session_ex = 1 # 8h
        freq_ex = 5    # 60 Hz
        
        r_mean = np.zeros((self.n_nmods, self.n_p_outs))
        r_ste = np.zeros((self.n_nmods, self.n_p_outs))

        # load data
        r_sq_r_learn = np.load(self.regression_data_dir + file_r_sq)

        for nm_i, nm in enumerate(self.nmods):
            for p_i, p_out in enumerate(self.p_outs): 
                r_mean[nm_i, p_i] = np.mean(r_sq_r_learn[session_ex,nm_i,freq_ex,:,:,p_i])
                r_ste[nm_i, p_i] = np.std(r_sq_r_learn[session_ex,nm_i,freq_ex,:,:,p_i])/np.sqrt(int(self.n_total_trials)*self.n_outputs)

        return r_mean, r_ste
        
    def _make_r_sq_vs_perc_conn_plot(self, file_r_sq, fig):
        '''
        plots r squared vs percent connectivity for data from given file
        
        Args:
            file_r_sq: path to data file
            fig: name of figure (e.g. Fig7, SuppFig5)
        
        Returns:
            mean and ste of r squared data
        '''
        
        figures = self.figures
        
        low_nm = 3    # 0.06
        high_nm = 9   # 0.16
        
        # compute mean R squared
        r_mean, r_ste = self._get_data_for_r_squared_vs_perc_conn_plot(file_r_sq)
        
        # t_value, confidence interval 95
        r_conf_95 = self.students_t_500 * r_ste

        # plot for low and high neuromodulator
        plt.figure(figsize=(6,6))
        for nm_i, nm in enumerate(self.nmods):
            if nm_i == low_nm:
                plt.plot(self.p_outs, r_mean[nm_i, :], self.colors[1][1], label='Low Concentration')
                plt.fill_between(self.p_outs, r_mean[nm_i, :] - r_conf_95[nm_i, :], 
                                 r_mean[nm_i, :] + r_conf_95[nm_i, :], facecolor=self.colors[1][1], alpha=0.2)
            elif nm_i == high_nm:
                plt.plot(self.p_outs, r_mean[nm_i, :], self.colors[3][1], label='High Concentration')
                plt.fill_between(self.p_outs, r_mean[nm_i, :] - r_conf_95[nm_i, :], 
                                 r_mean[nm_i, :] + r_conf_95[nm_i, :], facecolor=self.colors[3][1], alpha=0.2)
                
        plt.xlim(0.1,1)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.xlabel('Proportion output connectivity', fontsize=14)
        plt.ylabel('$R^2$', fontsize=14)
        plt.legend(loc='lower right', fontsize=14)
        fname = fig + 'd'
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
    
    def _get_data_for_random_walk_example(self, file_y_hat):
        '''
        gets the y_hat values for the random walk example
        
        Args:
            file_y_hat: path to data file
        
        Returns:
            average and ste of predicted output y_hat
        '''
        
        # compute means and standard devations for example: FREQ = 60Hz 
        session_ex = 1 # 8h
        freq_ex = 5    # 60 Hz
        output_ex = 2

        y_hat_store = np.load(self.regression_data_dir + file_y_hat)

        y_hat_mean = np.zeros((self.n_nmods, self.n_p_outs, self.n_bins))
        y_hat_ste = np.zeros((self.n_nmods, self.n_p_outs, self.n_bins))
        for nm_i, nm in enumerate(self.nmods):
            for p_i, p_out in enumerate(self.p_outs): 
                y_hat_mean[nm_i, p_i, :] = np.mean(y_hat_store[session_ex, nm_i, freq_ex, :, output_ex, p_i, :], axis=0)
                y_hat_ste[nm_i, p_i, :] = np.std(y_hat_store[session_ex, nm_i, freq_ex, :, output_ex, p_i, :], axis=0)/np.sqrt(int(self.n_total_trials))

        return y_hat_mean, y_hat_ste
    
    
    def _make_example_random_walk_plot(self, file_y, file_y_hat, fig):
        '''
        make random walk example plot
        
        Args:
            file_y: path to actual target functions
            file_y_hat: path to data file
            fig: name of figure
        '''
        
        nm_high = 9  # nmod concentration 0.18
        nm_low = 3 # nmod concentration 0.06
        p_i = 5
        output_ex = 2
        
        colors = self.colors
        figures = self.figures
        
        y_hat_mean, y_hat_ste = self._get_data_for_random_walk_example(file_y_hat)
        y_store = np.load(self.regression_data_dir + file_y)

        # t_value, confidence interval 95
        y_hat_conf_95 = self.students_t_50 * y_hat_ste
        
        times_ms = np.arange(100)
        plt.figure()
        plt.subplot(2,1,1)
        plt.plot(times_ms, y_hat_mean[nm_low, p_i, :], colors[1][1], label='Low Concentration')
        plt.fill_between(times_ms, y_hat_mean[nm_low, p_i, :] - y_hat_conf_95[nm_low, p_i, :], 
                         y_hat_mean[nm_low, p_i, :] + y_hat_conf_95[nm_low, p_i, :], facecolor=colors[1][1], alpha=0.1)
        plt.plot(times_ms, y_hat_mean[nm_high, p_i, :], colors[3][1], label='High Concentration')
        plt.fill_between(times_ms, y_hat_mean[nm_high, p_i, :] - y_hat_conf_95[nm_high, p_i, :], 
                         y_hat_mean[nm_high, p_i, :] + y_hat_conf_95[nm_high, p_i, :], facecolor=colors[3][1], alpha=0.1)

        plt.plot(times_ms, y_store[output_ex,:], color='darkgrey')
        plt.xlabel('')
        plt.xlim(0,100)
        plt.xticks([])
        plt.yticks([0, 5, 10])
        plt.ylabel('Output')

        plt.subplot(2,1,2)

        # absolute error
        error_low = y_hat_mean[nm_low, p_i, :] - y_store[output_ex,:]
        error_high = y_hat_mean[nm_high, p_i, :] - y_store[output_ex,:]

        plt.fill_between(times_ms, np.zeros(100), error_low, facecolor=colors[1][1], alpha=0.3)
        plt.fill_between(times_ms, np.zeros(100), error_high, facecolor=colors[3][1], alpha=0.3)
        plt.xlabel('Time (ms)')
        plt.ylabel('Absolute Error')
        plt.xlim(0,100)
        fname = fig + 'e'
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
        
    def make_figure_7(self):
        '''
        makes figure 7: computes regression or loads data from file, then makes all subplots
        '''
        # directories
        regression_data_dir = self.regression_data_dir
        
        # name of figure
        fig = 'Fig7'
        
        # suffix for storing data and figures, all neurons
        name = 'all'
        
        # all neurons
        neuron_start = 0
        neuron_stop = 1600
        
        # compute regression, or load data if already exists
        # returns filenames where data is stored
        file_r_sq, file_y, file_y_hat = self.compute_regression(neuron_start,
                                                                neuron_stop,
                                                                name)

        # make subfigures
        self._make_regression_rasters_plot(file_r_sq, fig, name)
        self._make_r_sq_vs_perc_conn_plot(file_r_sq, fig)
        self._make_example_random_walk_plot(file_y, file_y_hat, fig)
        
        
    def make_supplementary_figure_1a(self):
        '''
        make supplementary fig 1a: weight distribution, incoming and control
        '''
        # locally define params
        n_total_trials = self.n_total_trials
        avg_weights_file = self.avg_weights_file
        students_t_50 = self.students_t_50
        colors = self.colors
        figures = self.figures
        
        # first CA, then NON-CA
        columns = [np.arange(1,33,8), np.arange(3,33,8), np.arange(5,33,8), np.arange(7,33,8)]
        labels = ['None', 'Low', 'Moderate', 'High']
        groups = ['CA', 'Outgoing', 'Incoming', 'Control']
        
        # get number rows
        n_rows = len(np.genfromtxt(avg_weights_file, 
                                   delimiter=';', 
                                   skip_header=2))
        
        # initialise data matrix, n_rows x 17 columns 
        weight_dist_matrix = np.zeros((n_rows, 33))

        # get data from file
        for row in range(n_rows):
            weight_dist_matrix[row,:] = np.genfromtxt(avg_weights_file, 
                                                      delimiter=';', 
                                                      dtype=float, 
                                                      skip_header=2)[row]
        
        # get bin labels for weights, y-axis
        weights = weight_dist_matrix[:,0]
        
        # rescale bin labels to percentages
        weights = 100*weights/weights[0]

        # initialise data matrices for weight values
        weights_not_pot = np.zeros((len(groups), 4))
        weights_not_pot_conf_int = np.zeros((len(groups), 4))
        
        # plot for cell assembly and outgoing synapses (i.e. groups)
        for group_i, group in enumerate(groups):
            if group == 'Incoming' or group == 'Control':
                fig = plt.figure(figsize=(10,4))
                ax = fig.add_subplot(1,1,1)

                # go through each neuromodulator level
                for nm_i in range(len(labels)):
                    weight_dist = weight_dist_matrix[:,columns[group_i][nm_i]]
                    std = weight_dist_matrix[:,columns[group_i][nm_i]+1]

                    conf_int_95 = students_t_50 * std / np.sqrt(n_total_trials)

                    weights_low = weight_dist - conf_int_95
                    weights_high = weight_dist + conf_int_95

                    # Potentation under a couple percent considered "not potentiated"
                    # STD via Gaussian Error Propagation
                    # combine the first two bins (very low or no potentiation)

                    # weights
                    weights_not_pot[group_i, nm_i] = sum(weight_dist[0:2])
                    
                    # confidence interval 
                    conf_int = np.sqrt((students_t_50*std[0]/np.sqrt(n_total_trials))**2 + 
                                       (students_t_50*std[1]/np.sqrt(n_total_trials))**2)
                    weights_not_pot_conf_int[group_i, nm_i] = conf_int

                    # PLOT
                    ax.plot(weights, weight_dist, colors[nm_i][1], label=labels[nm_i])
                    ax.fill_between(weights, weights_low, weights_high, 
                                    facecolor=colors[nm_i][0], alpha=0.1)

                ax.set_yscale('log')
                plt.title('Weight distribution after 8h, ' + group)
                plt.xlabel('Mean weight (%)')
                plt.ylabel('Fraction of cells')
                plt.xlim(104, 182)
                plt.ylim(10e-8,0.56)
                if group == 'Outgoing':
                    plt.legend(loc = 'upper right')

                fname = 'SuppFig1a_' + group
                plt.savefig(figures+fname, format='svg')
                plt.show()
                
                # weights not potentiated 
                plt.figure(figsize=(7,5))
                plt.bar(labels, weights_not_pot[group_i,:], 
                        yerr=weights_not_pot_conf_int[group_i,:], width=0.7, 
                        color=['black', '#E7A42C', '#4AA42C', '#2C69AC'])
                plt.xticks(fontsize=20)
                plt.yticks(fontsize=20)
                fname = 'SuppFig1a_' + group + '_not_potentiated'
                plt.savefig(figures+fname, format='svg')
                plt.show()
        
        
    def make_supplementary_figure_1b(self):
        '''
        make supplementary fig 1b: firing rate distribution, CA
        '''
        n_total_trials = self.n_total_trials
        students_t_50 = self.students_t_50
        file = self.avg_rates_file
        figures = self.figures

        n_rows = len(np.genfromtxt(file, delimiter=';', skip_header=2))
        firing_rate_matrix = np.zeros((n_rows, 17))

        for row in range(n_rows):
            firing_rate_matrix[row,:] = np.genfromtxt(file, delimiter=';', 
                                                      dtype=float, 
                                                      skip_header=2)[row]

        rates = firing_rate_matrix[:,0]

        columns = [[3, 7, 11, 15], [1, 5, 9, 13]]
        labels = ['None', 'Low', 'Moderate', 'High']
        groups = ['Non-CA', 'CA']
        
        # non-CA
        group_i = 1
        group = groups[group_i]
        
        fig = plt.figure(figsize=(8,4))
        ax = fig.add_subplot(1, 1, 1)
        for nm_i in range(4):
            rate_dist = firing_rate_matrix[:,columns[group_i][nm_i]]
            std = firing_rate_matrix[:,columns[group_i][nm_i]+1]

            conf_int_95 = students_t_50 * std / np.sqrt(n_total_trials)

            to_plot = rate_dist>0
            rate_dist = rate_dist[to_plot]
            conf_int_95 = conf_int_95[to_plot]

            rate_low = rate_dist - conf_int_95
            rate_high = rate_dist + conf_int_95

            ax.plot(rates[to_plot], rate_dist, self.colors[nm_i][1], 
                    label=labels[nm_i], marker='.')
            ax.fill_between(rates[to_plot], rate_low, rate_high, 
                            facecolor=self.colors[nm_i][0], alpha=0.1)

        ax.set_yscale('log')
        plt.xlabel('Firing rate (Hz)')
        plt.ylabel('Fraction of cells')
        plt.title('Activity of core neurons')
        plt.legend()

        fname = 'SuppFig1b_' + group
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
        
    def make_supplementary_figure_1(self):
        '''
        makes supplementary figure 1a and 1b
        '''
        
        self.make_supplementary_figure_1a()
        self.make_supplementary_figure_1b()
        
    def make_supplementary_figure_3(self):
        '''
        makes supplementary figure 3
        '''
        
        # paths
        pca_data_dir = self.pca_data_dir
        figures = self.figures
        colors = self.colors
        
        # parameters
        session = '8h'
        nms = [0.06, 0.18]
        f = 60
        students_t_50 = self.students_t_50
        trial_ex = 0
        event = 'recall'
        event_i = self.events.index(event)
        size = (self.n_E, self.n_bins)
        var_thr = 0.7

        # set random seed
        np.random.seed(0)

        cumulative_var = np.zeros((2, self.n_total_trials, self.n_bins))
        raster_store = {}

        for nm_i, nm in enumerate(nms):

            # load event data
            event_data = self.load_event_data_obj(session, nm, f)

            # for each trial, compute pca and dimensionality
            for trial in range(int(self.n_total_trials)):

                # time stamps and ids to be binned
                ts = (event_data.ts[trial, event_i] - 
                      self.t_events[session][event])
                ids = event_data.ids[trial, event_i]

                # bin data
                ts_binned = self.bin_data(ts, self.n_bins)

                # make raster
                raster = self.make_raster(ts_binned, ids, size)

                # store rasters
                raster_store[nm_i, trial] = raster

                # do PCA
                pca = PCA(n_components=self.n_bins)
                pca.fit(raster.T)

                for k in range(self.n_bins):
                    cumulative_var[nm_i, trial, k] = np.sum(pca.explained_variance_ratio_[:k])
          
        # Supplementary Figure 3a
        # project into the space, LOW NM
        pca.fit(raster_store[0, trial_ex].T)
        data_proj = np.dot(pca.components_, raster_store[0, 0])

        std = np.std(data_proj)
        avg = np.mean(data_proj)

        plt.figure()
        plt.imshow(data_proj, aspect='auto', origin='lower', vmin=avg-.86, vmax=avg+.86)
        plt.xlabel('Time (ms)', fontsize=16)
        plt.ylabel('Principal component', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.colorbar()
        fname = 'SuppFig3a_pca_projection_f_60Hz_low_nm_06_origin_lower-avg+-std'
        plt.savefig(figures+fname, format='svg')
        plt.show()

        # Supplementary Figure 3b
        # project into the space, HIGH NM
        pca.fit(raster_store[1, trial_ex].T)
        data_proj = np.dot(pca.components_, raster_store[1, 0])

        std = np.std(data_proj)
        avg = np.mean(data_proj)

        plt.figure()
        plt.imshow(data_proj, aspect='auto', origin='lower', vmin=avg-std, vmax=avg+std)
        plt.xlabel('Time (ms)', fontsize=16)
        plt.ylabel('Principal component', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.colorbar()
        fname = 'SuppFig3b_pca_projection_f_60Hz_high_nm_18_origin_lower-avg+-std'
        plt.savefig(figures+fname, format='svg')
        plt.show()
        
        # Supplementary Figure 3c
        cumulative_var_mean = np.mean(cumulative_var, axis=1)
        cumulative_var_ste = np.std(cumulative_var, axis=1)/np.sqrt(self.n_total_trials)
        cumulative_var_95 = students_t_50 * cumulative_var_ste
        
        # array of dimensions for plotting
        dims = np.arange(self.n_bins)

        plt.figure()
        plt.plot(dims, var_thr*np.ones(self.n_bins), color='darkgrey', linestyle='--')
        plt.plot(dims, cumulative_var_mean[0,:], colors[1][1], label='Low Concentration')
        plt.fill_between(dims, cumulative_var_mean[0, :] - cumulative_var_95[0, :], 
                         cumulative_var_mean[0, :] + cumulative_var_95[0, :], 
                         facecolor=colors[1][1], alpha=0.1)
        plt.plot(dims, cumulative_var_mean[1,:], colors[3][1], label='High Concentration')
        plt.fill_between(dims, cumulative_var_mean[1, :] - cumulative_var_95[1, :], 
                         cumulative_var_mean[1, :] + cumulative_var_95[1, :], 
                         facecolor=colors[3][1], alpha=0.1)
        plt.xlim(0,self.n_bins)
        plt.xlabel('Number of PCs', fontsize=16)
        plt.ylabel('Proportion explained variance', fontsize=16)
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        plt.legend(fontsize=14)
        fname = 'SuppFig3c_cumulative_variance'
        plt.savefig(figures+fname, format='svg')
        plt.show()  
        
    
    def make_supplementary_figure_4(self):
        '''
        make supplementary fig 4: pca dimensionality for different bin sizes
        '''
        # parameters
        figures = self.figures
        processed_data = self.processed_data_dir
        pca_data_dir = self.pca_data_dir
        event = 'recall'
        
        # for 2ms and 4ms bins
        for n_bins_i in [50, 25]:
            size = (self.n_E, n_bins_i)
            
            file = 'pca_dim_mean_' + str(n_bins_i) + '.npy' 
            files = [f for f in os.listdir(pca_data_dir) if not f.startswith('.')]

            # if data already stored, load it, else compute pcas and store it
            if file in files:
                print('Data exists and will be loaded from ' + pca_data_dir)
                dim_mean = np.load(pca_data_dir + file)
            else:
                # compute the pca for all event data
                dim_mean = self.compute_pca(event, n_bins_i)
                
            print('---------------------------------------------------------')
            print('                     Bin size = ' + str(int(self.n_bins/n_bins_i)) + 'ms')
            print('---------------------------------------------------------')

            # plot means
            for session_i, session in enumerate(self.sessions):
                plt.figure()
                plt.title('Dimensionality at ' + session + ' for bin size ' + str(int(100/n_bins_i)) + ' ms' )
                plt.imshow(dim_mean[session_i,:,:].T, origin='lower', cmap='viridis', 
                           vmin=np.min(dim_mean), vmax=np.max(dim_mean))
                plt.xlabel('Neuromodulatior Concentration', fontsize=18)
                plt.ylabel('Stimulation Frequency', fontsize=18)
                plt.xticks(np.arange(10), self.nmods, fontsize=11)
                plt.yticks(np.arange(10), self.freqs, fontsize=11)
                plt.colorbar()
                fname = 'SuppFig4_nbins_' + str(n_bins_i) + '_' + session
                plt.savefig(figures+fname, format='svg')
                plt.show()

            plt.figure()
            plt.title('Difference ' + self.sessions[1] + ' - ' + self.sessions[0] + ' for bin size ' + str(int(100/n_bins_i)) + ' ms' )
            plt.imshow(dim_mean[1,:,:].T - dim_mean[0,:,:].T, origin='lower', 
                       cmap='seismic', vmin=-10, vmax=10)
            plt.xlabel('Neuromodulatior Concentration', fontsize=18)
            plt.ylabel('Stimulation Frequency', fontsize=18)
            plt.xticks(np.arange(10), self.nmods, fontsize=11)
            plt.yticks(np.arange(10), self.freqs, fontsize=11)
            plt.colorbar()
            fname = 'SuppFig4_nbins_' + str(n_bins_i) + '_diff'
            plt.savefig(figures+fname, format='svg')
            plt.show()
    
         
    def make_supplementary_figure_5(self):
        '''
        makes supplementary figure 5: computes regression or loads data from file
        '''
        # directories
        regression_data_dir = self.regression_data_dir
        
        # name of figure
        fig = 'SuppFig5_ca_'
        
        # suffix for storing data and figures, all neurons
        name = 'ca'
        
        # all neurons
        neuron_start = 0
        neuron_stop = 150
        
        # compute regression, or load data if already exists
        # returns filenames where data is stored
        file_r_sq, file_y, file_y_hat = self.compute_regression(neuron_start,
                                                                neuron_stop,
                                                                name)

        # make subfigures
        self._make_regression_rasters_plot(file_r_sq, fig, name)
        self._make_r_sq_vs_perc_conn_plot(file_r_sq, fig)
        self._make_example_random_walk_plot(file_y, file_y_hat, fig)
        
        
    def make_supplementary_figure_6(self):
        '''
        makes supplementary figure 6: computes regression or loads data from file
        '''
        # directories
        regression_data_dir = self.regression_data_dir
        
        # name of figure
        fig = 'SuppFig6_nc_'
        
        # suffix for storing data and figures, all neurons
        name = 'nc'
        
        # all neurons
        neuron_start = 150
        neuron_stop = 1600
        
        # compute regression, or load data if already exists
        # returns filenames where data is stored
        file_r_sq, file_y, file_y_hat = self.compute_regression(neuron_start,
                                                                neuron_stop,
                                                                name)

        # make subfigures
        self._make_regression_rasters_plot(file_r_sq, fig, name)
        self._make_r_sq_vs_perc_conn_plot(file_r_sq, fig)
        self._make_example_random_walk_plot(file_y, file_y_hat, fig)
    
    
    def _compute_pca_dimensionality_for_statistics(self, n_bins):
        # parameters
        neuron_start = 0
        neuron_stop = 1600
        var_thr = 0.7
        n_sessions = len(self.sessions)

        f = 60
        event = 'recall'
        event_i = self.events.index(event)

        size = (self.n_E, n_bins)

        nms = [0.0, 0.06, 0.12, 0.18]
        labels = ['none', 'low', 'mod', 'high']

        # compute PCA for each binned raster and check how many dimensions needed for var_thr percent variance
        #dim = np.zeros((n_sessions, len(nms), self.n_total_trials))
        dim = {}

        # go through each session, nmod, freq
        for session_i, session in enumerate(self.sessions):
            for nm_i, nm in enumerate(nms):

                # load event data
                event_data = self.load_event_data_obj(session, nm, f)

                # dummy array to hold dimensionality
                _dim = np.zeros(self.n_total_trials)

                for trial in range(self.n_total_trials):

                    # time stamps and ids to be binned
                    ts = (event_data.ts[trial, event_i] - 
                          self.t_events[session][event])
                    ids = event_data.ids[trial, event_i]

                    # bin data
                    ts_binned = self.bin_data(ts, n_bins)

                    # make raster
                    raster = self.make_raster(ts_binned, ids, size)

                    # do PCA
                    pca = PCA(n_components=n_bins)
                    pca.fit(raster.T)

                    k = 0
                    cumul_var = np.sum(pca.explained_variance_ratio_[:k])
                    while (k <= n_bins) & (cumul_var <= var_thr):
                        k += 1
                        cumul_var = np.sum(pca.explained_variance_ratio_[:k])

                    _dim[trial] = k  

                    #dim[session_i, nm_i, trial] = k
                dim[session, labels[nm_i]] = _dim
                
        np.save(self.pca_data_dir + 'pca_dim_for_stats_' + str(n_bins) + '.npy' , dim, allow_pickle=True)
    
    def compute_statistics_pca(self, n_bins, name):
        '''
        computes (or loads) and prints statistics on pca data
        '''
        labels = ['none', 'low', 'mod', 'high']
        all_files = [f for f in os.listdir(self.pca_data_dir) if not f.startswith('.')]
        
        file = 'pca_dim_for_stats_' + str(n_bins) + '.npy'
        if file in all_files:
            dim = np.load(self.pca_data_dir + file, allow_pickle=True)
            dim = dim.tolist()
        else:
            self._compute_pca_dimensionality_for_statistics(n_bins)
            dim = np.load(self.pca_data_dir + file, allow_pickle=True)
            dim = dim.tolist()
            
        # alpha
        alpha = 0.001
        
        # compute means
        mean_none_10s = np.mean(dim['10s', 'none'])
        mean_none_8h = np.mean(dim['8h', 'none'])
        mean_low_10s = np.mean(dim['10s', 'low'])
        mean_low_8h = np.mean(dim['8h', 'low'])
        mean_mod_10s = np.mean(dim['10s', 'mod'])
        mean_mod_8h = np.mean(dim['8h', 'mod'])
        mean_high_10s = np.mean(dim['10s', 'high'])
        mean_high_8h = np.mean(dim['8h', 'high'])

        # compute mann whitney u
        u_none = stats.mannwhitneyu(dim['8h', 'none'], dim['10s', 'none'])
        u_low = stats.mannwhitneyu(dim['8h', 'low'], dim['10s', 'low'])
        u_mod = stats.mannwhitneyu(dim['8h', 'mod'], dim['10s', 'mod'])
        u_high = stats.mannwhitneyu(dim['8h', 'high'], dim['10s', 'high'])
        
        # compare low vs moderate and low vs high at 8 hour recall
        #u_low_vs_mod = stats.mannwhitneyu(dim['8h', 'low'], dim['8h', 'mod'])
        #u_low_vs_high = stats.mannwhitneyu(dim['8h', 'low'], dim['8h', 'high'])
        
        print('_______________________________________________________________________________________________________________________')
        print(name + ', PCA dimensionality, 8h vs 10s, bin size = ' + str(int(100/n_bins)) + 'ms')
        print('_______________________________________________________________________________________________________________________')
        
        print(tabulate([[str(int(100/n_bins)) + 'ms', labels[0] + ' vs. ' + labels[0], mean_none_10s, mean_none_8h, mean_none_8h - mean_none_10s, 
                         'U=' + str(u_none[0]), 
                         'p=' + ('{:0.2e}'.format(u_none[1]) if u_none[1] < alpha else str(u_none[1].round(3))), 
                         '*' if u_none[1].round(3) <= alpha else 'n.s.'],
                        [str(int(100/n_bins)) + 'ms', labels[1] + ' vs. ' + labels[1], mean_low_10s, mean_low_8h, mean_low_8h - mean_low_10s, 
                         'U=' + str(u_low[0]), 
                         'p=' + ('{:0.2e}'.format(u_low[1]) if u_low[1] < alpha else str(u_low[1].round(3))), 
                         '*' if u_low[1].round(3) <= alpha else 'n.s.'],
                        [str(int(100/n_bins)) + 'ms', labels[2] + ' vs. ' + labels[2], mean_mod_10s, mean_mod_8h, mean_mod_8h - mean_mod_10s, 
                         'U=' + str(u_mod[0]),
                         'p=' + ('{:0.2e}'.format(u_mod[1]) if u_mod[1] < alpha else str(u_mod[1].round(3))), '*' if u_mod[1].round(3) <= alpha else 'n.s.'],
                        [str(int(100/n_bins)) + 'ms', labels[3] + ' vs. ' + labels[3], mean_high_10s, mean_high_8h, mean_high_8h - mean_high_10s, 
                         'U=' + str(u_high[0]),
                         'p=' + ('{:0.2e}'.format(u_high[1]) if u_high[1] < alpha else str(u_high[1].round(3))), '*' if u_high[1].round(3) <= alpha else 'n.s.'],
                       ], 
                       headers=['bin size', 'neuromodulator', 'mean 10s', 'mean 8h', 'mean 8h - mean 10s', 'mann whitney u', 'p-value', 'sig. at $\\alpha$ = 0.001'],
                       numalign='center', stralign='center', tablefmt='latex_raw'))
        
        print('\n')
    
    
    def compute_statistics_pca_compare_low_with_mod_and_high_at_8h(self):
        '''
        computes and prints statistics on pca data
        '''
        # parameters
        neuron_start = 0
        neuron_stop = 1600
        var_thr = 0.7
        n_sessions = len(self.sessions)

        f = 60
        event = 'recall'
        event_i = self.events.index(event)
        
        nms = [0.0, 0.06, 0.12, 0.18]
        labels = ['none', 'low', 'mod', 'high']
        
        all_files = [f for f in os.listdir(self.pca_data_dir) if not f.startswith('.')]
        
        table = []
        
        for n_bins in [100, 50, 25]:
            file = 'pca_dim_for_stats_' + str(n_bins) + '.npy'
            if file in all_files:
                dim = np.load(self.pca_data_dir + file, allow_pickle=True)
                dim = dim.tolist()
            else:
                self._compute_pca_dimensionality_for_statistics(n_bins)
                dim = np.load(self.pca_data_dir + file, allow_pickle=True)
                dim = dim.tolist()

            size = (self.n_E, n_bins)

            # alpha
            alpha = 0.001

            # compute means
            mean_low_8h = np.mean(dim['8h', 'low'])
            mean_mod_8h = np.mean(dim['8h', 'mod'])
            mean_high_8h = np.mean(dim['8h', 'high'])
            
            # compare low vs moderate and low vs high at 8 hour recall
            u_low_vs_mod = stats.mannwhitneyu(dim['8h', 'low'], dim['8h', 'mod'])
            u_low_vs_high = stats.mannwhitneyu(dim['8h', 'low'], dim['8h', 'high'])
            
            
            table.append([str(int(100/n_bins)) + 'ms', labels[1] + ' vs. ' + labels[2], mean_low_8h, mean_mod_8h, mean_mod_8h - mean_low_8h, 
                          'U=' + str(u_low_vs_mod[0]), 
                          'p=' + ('{:0.2e}'.format(u_low_vs_mod[1]) if u_low_vs_mod[1] < alpha else str(u_low_vs_mod[1].round(3))), 
                         '*' if u_low_vs_mod[1].round(3) <= alpha else 'n.s.'])
            table.append([str(int(100/n_bins)) + 'ms', labels[1] + ' vs. ' + labels[3], mean_low_8h, mean_high_8h, mean_high_8h - mean_low_8h, 
                          'U=' + str(u_low_vs_high[0]), 
                          'p=' + ('{:0.2e}'.format(u_low_vs_high[1]) if u_low_vs_high[1] < alpha else str(u_low_vs_high[1].round(3))), 
                          '*' if u_low_vs_high[1].round(3) <= alpha else 'n.s.'])
            
        
        print('_______________________________________________________________________________________________________________________')
        print('PCA dimensionality, low vs mod/high at 8h')
        print('_______________________________________________________________________________________________________________________')
        
        print(tabulate(table, headers=['bin size', 'neuromodulator', 'mean low', 'mean mod/high', 'mod/high - low', 'mann whitney u', 'p-value', 'sig. at $\\alpha$ = 0.001'], numalign='center', stralign='center', tablefmt='latex_raw'))
        
        #print(u_low_vs_mod)
        #print(u_low_vs_high)
        
        print('\n')
        
        
        
    def compute_statistics_stability(self):
        '''
        computes and prints statistics for spike time stability
        '''
        
        temporal_stability_data_dir = self.temporal_stability_data_dir
        file_stab = 'stability.npy' 
        file_ctrl = 'stability_ctrl.npy' 
        
        stability = np.load(temporal_stability_data_dir + file_stab)
        stability_ctrl = np.load(temporal_stability_data_dir + file_ctrl)
        
        labels = ['None', 'Low', 'Moderate', 'High']
        nms = [0.0, 0.06, 0.12, 0.18]
        n1=0
        n2=1600
        
        stab = np.zeros((2,4,50))
        for i in range(len(nms)):
            stab[0, i, :] = np.nanmean(stability[i, :, n1:n2], axis=1)
            stab[1, i, :] = np.nanmean(stability_ctrl[i, :, n1:n2], axis=1)

        # Check normality
        print('__________________________________________')
        print('Normality (Shapiro-Wilk)')
        print('__________________________________________')
        normal = np.zeros((2,4))
        for i in range(len(nms)):
            normal[0, i] = stats.shapiro(stab[0, i, :])[1]
            normal[1, i] = stats.shapiro(stab[1, i, :])[1]
            print(labels[i])
            print('  -real:      p=' + str(normal[0, i]))
            print('  -shuffled:  p=' + str(normal[1, i]))

        # Anova
        cond_labels = np.concatenate((np.repeat(['real'], 200), np.repeat(['shuffled'], 200)))
        nm_labels = np.tile(np.repeat(['none', 'low', 'mod', 'high'], 50), 2)

        df = pd.DataFrame({'cond': cond_labels,
                           'nmod': nm_labels,
                           'stability': stab.reshape((400,), order='C')
                         })

        model = ols('stability ~ C(cond) + C(nmod) + C(cond):C(nmod)', data=df).fit()
        table = sm.stats.anova_lm(model, typ=3)

        
        
        print('__________________________________________')
        print('ANOVA')
        print('__________________________________________')
        print(table)

        
        print('__________________________________________')
        print('Post-hoc tests')
        print('__________________________________________')
        # post hoc
        
        print('__________________________________________')
        print('Real vs. shuffled data')
        print('__________________________________________')
        # differences between real and shuffled
        print('none vs none: ', stats.ttest_ind(stab[0, 0, :], stab[1, 0, :]))
        print('low vs low: ', stats.ttest_ind(stab[0, 1, :], stab[1, 1, :]))
        print('mod vs mod: ', stats.ttest_ind(stab[0, 2, :], stab[1, 2, :]))
        print('high vs high: ', stats.ttest_ind(stab[0, 3, :], stab[1, 3, :]))
        
        print('__________________________________________')
        print('Real data - comparing neuromodulator levels')
        print('__________________________________________')
        # differences within real data
        # low vs mod, mod vs high, low vs high
        print('none vs low: ', stats.ttest_ind(stab[0, 0, :], stab[0, 1, :]))
        print('none vs mod: ', stats.ttest_ind(stab[0, 0, :], stab[0, 2, :]))
        print('none vs high: ', stats.ttest_ind(stab[0, 0, :], stab[0, 3, :]))
        print('low vs mod: ', stats.ttest_ind(stab[0, 1, :], stab[0, 2, :]))
        print('mod vs high: ', stats.ttest_ind(stab[0, 2, :], stab[0, 3, :]))
        print('low vs high: ', stats.ttest_ind(stab[0, 1, :], stab[0, 3, :]))
        
    
    def compute_statistics_stability_vs_weight(self):
        '''
        computes and prints correlations between spike time stability and 
        strength of weights arriving from the core
        '''
        
        print('________________________________________________________________________________')
        print('Correlations between spike time stability and strength of weights from assembly')
        print('________________________________________________________________________________')
        
        labels = ['None', 'Low', 'Moderate', 'High']
        nms = [0.0, 0.06, 0.12, 0.18]
        freq = 60

        # get weights
        weight = self.extract_outgoing_weights(nms, freq)

        stability = np.load(self.temporal_stability_data_dir + 'stability.npy')
        num_spks = np.load(self.temporal_stability_data_dir + 'n_spks_learn3.npy')

        n1 = 150
        n2 = 1600
        tr0 = 0
        tr1 = 50

        n_thr_min = 0

        for nm_i, nm in enumerate(nms):

            w = weight[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))
            n = num_spks[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))
            s = stability[nm_i, tr0:tr1, n1:n2].reshape((tr1-tr0)*(n2-n1))

            not_nan_inds = ~np.isnan(w) & ~np.isnan(s)
            n_cutoff_inds = (n > n_thr_min) #& (n < n_thr_max)

            w = w[not_nan_inds & n_cutoff_inds] # & pos_s_inds]
            n = n[not_nan_inds & n_cutoff_inds]
            s = s[not_nan_inds & n_cutoff_inds] # & pos_s_inds]

            if nm != 0:
                print(labels[nm_i])
                print('   corr: ', stats.spearmanr(w,s))
                print('   n: ', len(w))

    
    def compute_temporal_sequence_stats(self, nm, label):
        '''
        computes and prints mann whitney u results for temporal sequence storage
        '''
        
        freq_i = 5
        nm_i = list(self.nmods).index(nm)

        file_r_sq = 'r_sq_r_learn_all.npy'

        # directories
        regression_data_dir = self.regression_data_dir

        # load data
        r_sq_r_learn = np.load(regression_data_dir + file_r_sq) 

        # session_i, nm_i, f_j, trial, output, p_i
        r_10s = r_sq_r_learn[0,nm_i,freq_i,:,:,:].mean(axis=(1,2)).reshape(50,)
        r_8h = r_sq_r_learn[1,nm_i,freq_i,:,:,:].mean(axis=(1,2)).reshape(50,)

        # means
        r_mean_10s = r_10s.mean()
        r_mean_8h = r_8h.mean()

        mean_diff = r_mean_8h - r_mean_10s

        # stds
        r_std_10s = r_10s.std()
        r_std_8h = r_8h.std()

        # mann whitney u
        u = stats.mannwhitneyu(r_8h, r_10s)

        print('_________________________________')
        print(label + ' neuromodulator level')
        print('_________________________________')
        print('Mean difference (8h - 10s): ' + str(mean_diff))
        print('Mann-Whitney U: ' + str(u))
    
    def compute_statistics(self):
        '''
        computes and prints all statistics from the paper
        '''
        print('=========================================================================')
        print('Computing PCA statistics (Figure 5 and Supplementary Figure S4)....')
        print('=========================================================================')
        print('\n')
        self.compute_statistics_pca(n_bins=self.n_bins, name='Figure 5')
        self.compute_statistics_pca(n_bins=50, name='Supplementary Figure 4')
        self.compute_statistics_pca(n_bins=25, name='Supplementary Figure 4')
        self.compute_statistics_pca_compare_low_with_mod_and_high_at_8h()
        
        print('\n')
        print('=========================================================================')
        print('Computing spike time stability statistics (Figure 6)....')
        print('=========================================================================')
        print('\n')
        self.compute_statistics_stability()
        self.compute_statistics_stability_vs_weight()
        
        print('\n')
        print('=========================================================================')
        print('Computing temporal sequences storage statistics (Figure 7)....')
        print('=========================================================================')
        print('\n')
        self.compute_temporal_sequence_stats(0.06, 'Low')
        self.compute_temporal_sequence_stats(0.18, 'High')