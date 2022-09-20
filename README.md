# Memory consolidation in recurrent spiking neural networks


## Outline

This package serves to simulate recurrent spiking neural networks with calcium-based synaptic plasticity and synaptic tagging and capture.
The generic C++ code and build scripts for specific simulations are located in the directory __simulation-code/__. Building will create binaries in the directory __simulation-bin/__, which contains scripts to run the specific simulations. Building has been tested with g++ 10.3.0 and boost 1.77.0 on Linux kernel 5.13.0. Besides that, pre-compiled binaries are available as part of the latest release.

The directory __analysis/__ contains Python scripts serving to analyze the data produced by the simulations.

The package that is provided here has been developed and used for the following publications:

1. Luboeinski, J., Tetzlaff, C. Memory consolidation and improvement by synaptic tagging and capture in recurrent neural networks. Commun. Biol. 4, 275 (2021). https://doi.org/10.1038/s42003-021-01778-y ---
   _original model underlying the simulation code, investigation of synaptic consolidation and improvement of a single memory representation_

2. Luboeinski, J., Tetzlaff, C. Organization and priming of long-term memory representations with two-phase plasticity. Cogn. Comput. (2022). https://doi.org/10.1007/s12559-022-10021-7 ---
   _extension of the model, investigation of functional consequences of interactions between multiple memory representations in different paradigms_

3. Lehr, A.B., Luboeinski, J., Tetzlaff, C. Neuromodulator-dependent synaptic tagging and capture retroactively controls neural coding in spiking neural networks. Sci. Rep., under review (2022). ---
   _extension of the model, investigation of the effects of neuromodulator-dependent protein synthesis on memory consolidation_

4. Luboeinski, J. The Role of Synaptic Tagging and Capture for Memory Dynamics in Spiking Neural Networks \[Dissertation\]. University of Göttingen (2021). https://doi.org/10.53846/goediss-463 ---
   _extension of the model, further investigation of multiple memory representations with attractor dynamics as well as characterization of plasticity regimes depending on pre- and post-synaptic firing rate_

Please cite accordingly if you use parts of the software package or the model for your research (BibTeX code is found [here](BIBTEX.md)). Note that the simulation code provided here contains some features that have not been used in any publication yet. Please feel free to contact us for any questions.


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
	* __compile_sizes__ - compiles the code for network simulations of learning, consolidation, and recall of a memory representation of certain size
	
	
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
	* __compile_2N1S_conv__ - compiles the code for testing the convergence of the membrane potential
	* __compile_2N1S_Li2016__ - compiles the code for simulating the effect of basic plasticity induction protocols at a single synapse, with the same model as in Li, Kulvicius, Tetzlaff, PLOS ONE, 2016
	* __compile_2N1S_minimal__ - compiles the code for a minimal example of the induction of early-phase plasticity by a few pre-defined spikes
	* __compile_activation_attractors__ - compiles the code for network simulations to investigate the spontaneous activation of assemblies in the absence of plasticity and in the presence of background noise or in the presence of 1 Hz/5 Hz oscillatory input to the inhibitory population
	* __compile_CA150__ - compiles the code for network simulations of learning, consolidation, and recall of a memory representation of 150 neurons
	* __compile_max_activity__ - compiles the code for a network that has one neuron spiking at maximal activity
	* __compile_onespike__ - compiles the code for a network that is stimulated with a single pulse to evoke one spike in one neuron
	* __compile_organization_attractors__ - compiles the code for network simulations of learning and consolidating three attractor memory representations in different organizational paradigms
	* __compile_PFreq__ - compiles the code for network simulations to characterize plasticity regimes depending on pre- and post-synaptic firing rate

The simulation is run by executing the binary file with or without command line options (as defined in __NetworkMain.cpp__, e.g., via one of the following shell scripts).
Please note that in addition there are preprocessor options that can be set before compiling (e.g., in __NetworkSimulation.cpp__) but cannot be changed during
runtime.

The binaries and run scripts for the studies mentioned above are located in subdirectories of __simulation-bin/__. Please note that some of these scripts trigger a cascade of many simulations by using the `screen` command. This may cause less powerful machines to take very long or to run into memory issues. In those cases, you might consider to run simulations separately.

* __run_binary_paper1/__:

	* __run_2N1S__ - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity (see, for example, Sajikumar et al., J Neurosci, 2005)
	* __run_full__ - learn a memory representation, let it consolidate, and recall after 8 hours (no fast-forwarding, takes very long)
	* __run_IRS__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, apply intermediate stimulation, let the memory representation consolidate, and recall after 8 hours
	* __run_varied_inhibition__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours; do this for varied inhibition parameters
	* __run_varied_size__ - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours; do this for varied cell assembly size
	* __connections.txt__ - the default connectivity matrix used for this paper; if this file is absent, the simulation program will automatically generate a new network structure

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

	* __run_2N1S_conv__ - tests the convergence of the neuronal membrane potential following current stimulation
	* __run_2N1S_Li2016__ - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity, with the same model as in Li, Kulvicius, Tetzlaff, PLOS ONE, 2016
	* __run_2N1S_minimal__ - minimal example of the induction of early-phase plasticity by a few pre-defined spikes, with stochastic plasticity dynamics
	* __run_2N1S_minimal_det__ - minimal example of the induction of early-phase plasticity by a few pre-defined spikes, with deterministic plasticity dynamics
	* __run_activation_attractors__ - simulate the activity in a previously consolidated network for 3 minutes without plasticity (it is required to run __run_learn_cons_attractors__ beforehand)
	* __run_learn_cons_attractors__ - subsequently learn 3 attractor memory representations; consolidate for 8 hours after learning each assembly
	* __run_PFreq__ - network simulations to characterize plasticity regimes depending on pre- and post-synaptic firing rate
	* __connections.txt__ - the default connectivity matrix used for paper 1
	
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
* __extractParamsFiringRates.py__ - recursively extracts the mean firing rates of neuronal subpopulations along with the simulation parameters from directories containing spike raster data (can be used to process many datasets)
* __extractParamsMeanWeights.py__ - recursively extracts the mean weights across neuronal subpopulations along with the simulation parameters from directories containing `[timestamp]_net_[time].txt` files (can be used to process many datasets)
* __extractParamsProteins.py__ - recursively extracts the mean protein amount across core and non-core neurons along with the simulation parameters from directories containing `[timestamp]_mean_weight.txt` files (can be used to process many datasets)
* __extractParamsQMI.py__ - recursively extracts the Q and MI measures along with the simulation parameters from directories containing `[timestamp]_net_[time].txt` files (can be used to process many datasets)
* __extractParamsQMIfromSpikes.py__ - recursively extracts the Q and MI measures along with the simulation parameters from directories containing spike raster data (can be used to process many datasets)
* __frequencyAnalysisSpikeRaster.py__ - computes the frequency spectrum of spike raster data
* __meanCorrelations.py__ - computes firing rate correlations for neuron pairs from spike raster data and averages over subpopulations of the network
* __nmAnalysisClass.py__ - generates temporal traces and analyzes their consolidation; produces most of the final data for paper 3 (an interactive frontend for this is found in __notebooks/__)
* __numberOfSpikesInBins.py__ - computes the distribution of spikes per time bin from cell assembly time series data
* __overlapParadigms.py__ - defines paradigms of overlapping cell assemblies
* __utilityFunctions.py__ - diverse utility functions, e.g., to read firing rate, early- and late-phase weight data from `[timestamp]_net_[time].txt` files produced by the simulation program
* __valueDistributions.py__ - functions to analyze and plot weight and firing rate distributions


## Interactive scripts

The following subfolders of __notebooks/__ contain interactive Jupyter notebooks serving to reproduce the figures of a related paper. They were tested to run with Jupyter Lab 3.2.8 and 3.4.3. Raw data has to be used from provided sources or to be generated using the simulation code.

### Files

* __lehr_luboeinski_tetzlaff_2022/__ - reproduces figures of paper 3
