# Memory consolidation in recurrent spiking neural networks


## Outline

This package serves to simulate recurrent spiking neural networks with calcium-based synaptic plasticity and synaptic tagging and capture.
The generic C++ program code and build scripts for specific simulations are located in the directory __simulation-code/__. Building will create binaries in the directory __simulation-bin/__, which contains scripts to run the specific simulations. Building has most recently been tested with g++ 10.5.0 and boost 1.77.0. Using a Linux system is recommended.

The directory __analysis/__ contains Python scripts serving to analyze the data produced by the simulations.

The directory __notebooks/__ contains Jupyter notebooks serving to reproduce data with a graphical user interface.

The package that is provided here has been developed and used for a number of publications (see the list [here](PUBLICATIONS.md)).
Please cite accordingly if you use parts of the code or the model for your research (BibTeX code can be found [here](BIBTEX.md)). 
Also note that the simulation code provided here contains some features that have not been used in any publications yet. Please feel free to [contact us](mailto:mail@jlubo.net) for any questions or comments!


## Simulation code

The code for running the main simulation is located in the directory __simulation-code/__. It comes with build scripts for specific cases.

### Files

* __NetworkMain.cpp__ - contains the main function initializing network simulations
* __NetworkSimulation.cpp__ - class performing network simulations
* __Network.cpp__ - class describing the network
* __Neuron.cpp__ - class describing one neuron
* __Stimulus.cpp__ - class describing a stimulus
* __StimulusProtocols.cpp__ - class to define specific stimulus protocols
* __Definitions.hpp__ - general definitions
* __SpecialCases.hpp__ - definitions for special simulations (see this to reproduce results of the studies mentioned above)
* __Tools.cpp__ - collection of utility functions
* __Plots.cpp__ - collection of plotting functions employing *gnuplot*
* __plotFunctions.py__ - collection of plotting functions employing *Matplotlib*

### Compiling and linking

The simulation code comes with shell scripts (in addition to the included Makefile) to build it for different purposes, related to the studies mentioned above. The subdirectories of __simulation-code/__ contain the following scripts:

* __build_scripts_paper1/__:

	* __compile_2N1S__ - compiles the code for simulating the effect of basic plasticity induction protocols at a single synapse
	* __compile_IRS__ - compiles the code for network simulations of learning, consolidation, recall, and the effect of intermediate stimulation on a memory representation 
	* __compile_recall_varied_size__ - compiles the code for network simulations of learning, consolidation, and recall of a memory representation of certain size
	
	
* __build_scripts_paper2/__:

	* __compile_activation__ - compiles the code for network simulations to investigate the spontaneous activation of assemblies in the absence of plasticity and in the presence of background noise
	* __compile_organization__ - compiles the code for network simulations of learning and consolidating three memory representations in different organizational paradigms
	* __compile_organization_IC__ - compiles the code for network simulations of learning and consolidating three memory representations in different organizational paradigms, with intermediate consolidation
	* __compile_organization_noLTD__ - compiles the code for network simulations of learning and consolidating three memory representations in different organizational paradigms, without LTD
	* __compile_organization_randweight__ - compiles the code for network simulations of learning and consolidating three memory representations in different organizational paradigms, with randomly initialized weights
	* __compile_recall__ - compiles the code for network simulations to investigate in the absence of plasticity the recall of different assemblies


* __build_scripts_paper3/__:

	* __compile_nm_psth_c__ - compiles the code for network simulations of learning, neuromodulator-dependent consolidation, and recall of a memory representation of 150 neurons
	

* __build_scripts_misc/__:
	* __compile__ - compiles the code as it is (without setting any specific preprocessor definition for a particular simulation)
	* __compile_2N1S__ - compiles the code for simulating the effect of basic plasticity induction protocols at a single synapse
	* __compile_2N1S_conv__ - compiles the code for testing the convergence of the membrane potential
	* __compile_2N1S_Li2016__ - compiles the code for simulating the effect of basic plasticity induction protocols at a single synapse, with the same model as in Li, Kulvicius, Tetzlaff, PLOS ONE, 2016
	* __compile_2N1S_basic_early__ - compiles the code for a simple example of the induction of early-phase plasticity by a few pre-defined spikes
	* __compile_2N1S_basic_late__ - compiles the code for a simple example of the induction of late-phase plasticity by prolonged substantial stimulation of one neuron
	* __compile_activation_attractors__ - compiles the code for network simulations to investigate the spontaneous activation of assemblies in the absence of plasticity and in the presence of background noise or in the presence of 1 Hz/5 Hz oscillatory input to the inhibitory population
	* __compile_CA200__ - compiles the code for network simulations of learning, consolidation, and recall of a memory representation of 200 neurons
	* __compile_max_activity__ - compiles the code for a network that has one neuron spiking at maximal activity
	* __compile_onespike__ - compiles the code for a network that is stimulated with a single pulse to evoke one spike in one neuron
	* __compile_organization_attractors__ - compiles the code for a network that learns and consolidates three attractor memory representations in different organizational paradigms
	* __compile_PFreq__ - compiles the code for simulations serving to characterize plasticity regimes depending on pre- and post-synaptic firing rate
	* __compile_recall_varied_size__ - compiles the code for network simulations of learning, consolidation, and recall of a memory representation of certain size
	* __compile_smallnet_ou__ - compiles the code for a small network that is stimulated with an Ornstein-Uhlenbeck current

### Running a simulation
To run a simulation, execute the binary file with or without command line options (as defined in __NetworkMain.cpp__, e.g., via one of the following shell scripts).
Please note that in addition, there are preprocessor options that can be set before compiling (see, for example, __NetworkSimulation.cpp__ or one of the __compile\*__ scripts) but that cannot be changed during
runtime.

The binaries and run scripts for the studies mentioned above are located in subdirectories of __simulation-bin/__. Please note that some of these scripts trigger a cascade of many simulations by using the `screen` command (which has to be installed). This may cause less powerful machines to take very long or to run into memory issues. In those cases, you might consider to run simulations separately.

* __run_binary_paper1/__:

	* __run_2N1S__ - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity (see, for example, Sajikumar et al., J Neurosci, 2005)
	* __run_IRS__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, apply intermediate stimulation, let the memory representation consolidate, and recall after 8 hours
	* __run_recall_150_full__ - learn a memory representation, let it consolidate, and recall after 8 hours (no fast-forwarding, takes very long)
	* __run_recall_varied_inhibition__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours; do this for varied inhibition parameters
	* __run_recall_varied_size__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours; for varied pattern size
	* __connections.txt__ - the default connectivity matrix used in this paper; if this file is absent, the simulation program will automatically generate a new connectivity structure


* __run_binary_paper2/__:

	* __run_activation*__ - simulate the activity in a previously consolidated network for 3 minutes without plasticity (it is required to run the according __run_learn_cons*__ script beforehand)
	* __run_learn_cons__ - subsequently learn 3 memory representations and let them consolidate for 8 hours
	* __run_learn_cons_interleaved__ - learn 3 memory representations in an interleaved manner and let them consolidate for 8 hours (it is required to run __run_learn_cons__ beforehand)
	* __run_learn_cons_IC__ - subsequently learn 3 memory representations; intermediate consolidation for 8 hours after learning each individual assembly
	* __run_learn_cons_noLTD__ - subsequently learn 3 memory representations and let them consolidate for 8 hours; without LTD
	* __run_learn_cons_randweight__ - subsequently learn 3 memory representations and let them consolidate for 8 hours; with randomly initialized weights
	* __run_priming_and_activation__ - prime one of the assemblies in a previously consolidated network at a certain time and then simulate the activity for 3 minutes without plasticity (it is required to run __run_learn_cons_interleaved__ beforehand)
	* __run_recall__ - apply recall stimuli to the assemblies in a previously consolidated network and in a control network (it is required to run __run_learn_cons__ beforehand)


* __run_binary_paper3/__:

	* __run_raster_10s__ - learn a memory representation (with varied stimulation strength) and recall after 10 seconds
	* __run_raster_8h__ - learn a memory representation (with varied stimulation strength), let it consolidate (with varied level of neuromodulation), and recall after 8 hours
	* __run_nm_timing__ - learn a memory representation, let it consolidate (with low or high level of neuromodulation and varied timing), and recall after 8 hours


* __run_binary_misc/__:

	* __run_2N1S_basic_early__ - simple example of the induction of early-phase plasticity by a few pre-defined spikes, with stochastic plasticity dynamics (the interactive notebook in __notebooks/simulator_comparison_basic/__ can be used for the same purpose)
	* __run_2N1S_basic_early_det__ - simple example of the induction of early-phase plasticity by a few pre-defined spikes, with deterministic plasticity dynamics
	* __run_2N1S_basic_late__ - simple example of the induction of late-phase plasticity by prolonged substantial stimulation of one neuron (the interactive notebook in __notebooks/simulator_comparison_basic/__ can be used for the same purpose)
	* __run_2N1S_conv__ - tests the convergence of the neuronal membrane potential following current stimulation
	* __run_2N1S_facilitated_spiking__ - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity, with facilitated spiking (hence, with more postsynaptic spikes)
	* __run_2N1S_Li2016__ - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity, with the same model as in Li, Kulvicius, Tetzlaff, PLOS ONE, 2016
	* __run_activation_attractors__ - simulate the activity in a previously consolidated network for 3 minutes without plasticity (it is required to run __run_learn_cons_attractors__ beforehand)
	* __run_benchmark__ - pipeline for benchmarks of runtime and memory usage; can be used with different paradigms ('CA200', '2N1S\_basic\_late', ...)
	* __run_learn_cons_attractors__ - subsequently learn 3 attractor memory representations; consolidate for 8 hours after learning each assembly
	* __run_defaultnet_onespike_exc__ - test case to study the transmission of a single spike of an excitatory neuron in a network
	* __run_defaultnet_onespike_inh__ - test case to study the transmission of a single spike of an inhibitory neuron in a network
	* __run_PFreq__ - network simulations to characterize plasticity regimes depending on pre- and post-synaptic firing rate
	* __run_smallnet2_det_max_activity__ -  simulating a small network of 4 excitatory neurons with deterministic dynamics, where one neuron fires at maximum rate
	* __run_smallnet3_8h-recall__ - simulating a small network of 4 excitatory and 1 inhibitory neurons, where one excitatory neuron receives typical "learning" stimulation, and "recall" stimulation after 8 hours (the interactive notebook in __notebooks/simulator_comparison_basic/__ can be used for the same purpose)
	* __run_smallnet3_10s-recall__ - simulating a small network of 4 excitatory and 1 inhibitory neurons, where one excitatory neuron receives typical "learning" stimulation, and "recall" stimulation after 10 seconds (the interactive notebook in __notebooks/simulator_comparison_basic/__ can be used for the same purpose)
	* __run_recall_varied_size__ - learn a memory representation and recall after 10 seconds; learn a memory representation, let it consolidate, and recall after 8 hours; do this for different pattern sizes (the interactive notebook in __notebooks/simulator_comparison_memory_recall/__ can be used for the same purpose)
	* __runner_recall_varied_size_desktop__ - run a number of trials of __run_recall_varied_size__ via screen
	* __runner_recall_varied_size_gwdg-medium_ - run a number of trials of __run_recall_varied_size__ on the SLURM-managed GWDG SCC cluster
	* __track_allocated_memory__ - script to track the memory usage of a given process


## Analysis scripts

The following scripts, located in __analysis/__, serve to process and analyze the data produced by the simulation code. They were tested to run with Python 3.7.3, NumPy 1.20.1, SciPy 1.6.0, and pandas 1.0.3.
Please note that some of the script files are interdependent.
Also note that not all script files and functions have to be used to reproduce the results of any single study mentioned above (also see section "Interactive scripts" below).

### Files

* __adjacencyFunctions.py__ - functions to analyze the connectivity and weights in a network (used to compute mean and standard deviation of early- and late-phase weights)
* __analyzeWeights.py__ - routine that runs functions to investigate the synaptic weight structure of networks (reads from `[timestamp]_net_[time].txt` files produced by the simulation program)
* __assemblyAttractorStatistics.py__ - determines the statistics of the activation of attractor cell assemblies (considering exclusive activation and transitions between attractors)
* __assemblyAvalancheStatistics.py__ - determines the statistics of avalanche occurrence within cell assemblies
* __averageFileColumnsAdvanced.py__ - averages data columns across files (for example, average over multiple weight traces or probability distributions)
* __averageWeights-py__ - averages across multiple weight matrices
* __calculateMIa.py__ - calculates the mutual information from two firing rate distributions (given either by `[timestamp]_net_[time].txt` files or by arrays)
* __calculateQ.py__ - calculates the pattern completion coefficient Q for an input-defined cell assembly (either from the firing rate distribution given by a `[timestamp]_net_[time].txt` file or from mean firing rates)
* __computeRateFromSpikes.py__ - computes the firing rate over time via fixed time windows from spike raster data
* __extractAverageQMI.py__ - averages over trials of Q and MI data
* __extractParamsFiringRates.py__ - recursively extracts the mean firing rates of neuronal subpopulations along with the simulation parameters from directories containing spike raster data (can be used to process many datasets)
* __extractParamsMeanWeights.py__ - recursively extracts the mean weights across neuronal subpopulations along with the simulation parameters from directories containing `[timestamp]_net_[time].txt` files (can be used to process many datasets)
* __extractParamsProteins.py__ - recursively extracts the mean protein amount across core and non-core neurons along with the simulation parameters from directories containing `[timestamp]_mean_weight.txt` files (can be used to process many datasets)
* __extractParamsQMI.py__ - recursively extracts the Q and MI measures along with the simulation parameters from directories containing `[timestamp]_net_[time].txt` files (can be used to process many datasets)
* __extractParamsQMIfromSpikes.py__ - extracts simulation parameters and firing rates and/or Q and MI measures from spike raster data; optionally recursively searches directories for suitable data files (can be used to process many datasets)
* __frequencyAnalysisSpikeRaster.py__ - computes the frequency spectrum of spike raster data
* __meanCorrelations.py__ - computes firing rate correlations for neuron pairs from spike raster data and averages over subpopulations of the network
* __nmAnalysisClass.py__ - generates temporal traces and analyzes their consolidation; produces most of the final data for paper 3 (["Neuromodulator-dependent synaptic tagging..."](https://doi.org/10.1038/s41598-022-22430-7); an interactive frontend for this can be found in __notebooks/__)
* __numberOfSpikesInBins.py__ - computes the distribution of spikes per time bin from cell assembly time series data
* __overlapParadigms.py__ - defines paradigms of overlapping cell assemblies
* __plotQMICoreSize.py__ - functions to plot Q and MI values (from 10s- and 8h-recall) over assembly core size
* __plotSimResultsComparisonMeanSEM.py__ - plots a comparison of traces to compare across different simulators (adopted from [here](https://github.com/jlubo/simulator_comparison))
* __utilityFunctions.py__ - diverse utility functions, e.g., to read firing rate, early- and late-phase weight data from `[timestamp]_net_[time].txt` files produced by the simulation program
* __valueDistributions.py__ - functions to analyze and plot weight and firing rate distributions


## Interactive scripts

The following subfolders of __notebooks/__ contain interactive Jupyter notebooks serving to reproduce figures for a related paper. They were tested to run with Jupyter Lab 3.2.8 and 3.4.3.

### Files

* __lehr_luboeinski_tetzlaff_2022/__ - reproduces figures of paper 3 (["Neuromodulator-dependent synaptic tagging..."](https://doi.org/10.1038/s41598-022-22430-7)); raw data has to be used from provided sources or to be generated using the simulation code
* __simulator_comparison_basic/__ - runs basic simulations (basic early- and late phase plasticity, `smallnet3` dynamics) to compare with other simulators, in particular, with [Arbor](https://github.com/jlubo/arbor_network_consolidation) or [Brian 2](https://github.com/jlubo/brian_synaptic_plasticity_stc)
* __simulator_comparison_memory_recall/__ - runs simulations of memory recall (after 10 seconds and after 8 hours) to compare with other simulators, in particular, with [Arbor](https://github.com/jlubo/arbor_network_consolidation) or [Brian 2](https://github.com/jlubo/brian_network_plasticity)
