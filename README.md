# Memory consolidation in recurrent spiking neural networks


## Outline

This package serves to simulate recurrent spiking neural networks with calcium-based synaptic plasticity and synaptic tagging and capture. 
The C++ code for these simulations is located in the directory 'simulation-code/'. 

Pre-compiled binaries of the simulation code for Linux (tested with Kernel version 5.4) can be found in 'simulation-bin/', compiled and linked with g++ in version 7.4.0 and boost in version 1.65.1.

Furthermore, the package contains Python scripts to analyze the data produced by the simulations. These scripts are located in the directory 'analysis/'. 

Please kindly cite our papers if you use the code in your research:

1. Luboeinski, J., Tetzlaff, C. Memory consolidation and improvement by synaptic tagging and capture in recurrent neural networks. Commun. Biol. 4, 275 (2021). https://doi.org/10.1038/s42003-021-01778-y

2. Luboeinski, J., Tetzlaff, C. Organization and priming of long-term memory representations with two-phase plasticity. bioRxiv (2021).
	
The first paper presents the original model that underlies the simulation code and derives findings for the synaptic consolidation of a single memory representation.
The second paper extends the model to investigate the interaction of multiple memory representations in different paradigms.
Additionally, the code contains features that have not been used in publications yet. Please feel free to contact us on questions about the code and on further investigations that can be done with it.


## Simulation code

### Files

* 'NetworkMain.cpp' - main function that initializes network simulations
* 'NetworkSimulation.cpp' - class performing network simulations
* 'Network.cpp' - class describing the network
* 'Neuron.cpp' - class describing one neuron
* 'Stimulus.hpp' - class describing a stimulus
* 'StimulusProtocols.hpp' - class to define specific stimulus protocols
* 'Definitions.hpp' - general definitions
* 'SpecialCases.hpp' - definitions for special simulations (see this to reproduce results of the papers mentioned above)
* 'Tools.hpp' - collection of utility functions
* 'Plots.hpp' - collection of plotting functions employing *gnuplot*
* 'plotFunctions.py' - collection of plotting functions employing *Matplotlib*

### Compiling and linking

To build the simulation code, in addition to the included Makefile, the code comes with shell scripts for different purposes. For each paper, the related scripts are located in a specific subdirectory of 'simulation-code/':
* 'build\_scripts\_paper1':
	* 'compile_sizes' - compiles the code for performing network simulations to learn, consolidate, and recall a memory representation of a certain size
	* 'compile_2N1S' - compiles the code for applying basic plasticity induction protocols to a single synapse
	* 'compile_IRS' - compiles the code for performing network simulations to learn and consolidate a memory representation, apply intermediate stimulation, and recall

* 'build\_scripts\_paper2':
	* 'compile_organization' - compiles the code to learn and consolidate three memory representations in different organizational paradigms
	* 'compile\_organization_noLTD' - compiles the code to learn and consolidate three memory representations in different organizational paradigms, without LTD
	* 'compile_activation' - compiles the code to investigate the spontaneous activation in a network in the absence of plasticity
	* 'compile_recall' - compiles the code to investigate the recall of different assemblies in the absence of plasticity

### Running the simulation

The simulation is run by executing the binary file with or without command line options (as defined in 'NetworkMain.cpp', e.g., via one of the following shell scripts).
Please note that there are additional preprocessor options that have to be set in 'NetworkSimulation.cpp' before compiling and cannot be changed during
runtime. 

The binaries and run scripts for the papers mentioned above are located in specific subdirectories of 'simulation-bin/':
* 'run\_scripts\_paper1':
	* 'run_sizes' - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours
	* 'run_IRS' - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, apply intermediate stimulation, let the memory representation consolidate, and recall after 8 hours
	* 'run_full' - learn a memory representation, let it consolidate, and recall after 8 hours (no fast-forwarding, takes very long)
	* 'run_2N1S' - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity
	* 'connections.txt' - the default connectivity matrix used in this paper; if this file is absent, the simulation program will automatically generate a new network structure

* 'run\_scripts\_paper2':
	* 'run\_learn\_cons' - subsequently learn 3 memory representations and let them consolidate for 8 hours
	* 'run\_learn\_cons\_noLTD' - subsequently learn 3 memory representations and let them consolidate for 8 hours, without LTD
	* 'run\_activation' - simulate the activity in a previously consolidated network for 3 minutes without plasticity (required to run 'run\_learn\_cons' first)
	* 'run\_priming\_and\_activation' - prime one of the assemblies in a previously consolidated network at a certain time and then simulate the activity for 3 minutes without plasticity (required to run 'run\_learn\_cons' first)
	* 'run\_recall' - apply recall stimuli to the assemblies in a previously consolidated network and in a control network (required to run 'run\_learn\_cons' first)

## Analysis scripts

The following scripts serve to process and analyze the data produced by the simulation code. They were tested to run with Python 3.7.3, NumPy 1.20.1, SciPy 1.6.0, and pandas 1.0.3. 
Note that some of the script files depend on others. 
Also note that to reproduce the results of only one of the papers mentioned above, not all script files and functions will be required.

### Files

* 'adjacencyFunctions.py' - functions to analyze the connectivity and weights (used to compute mean and standard deviation of early- and late-phase weights)
* 'analyzeWeights.py' - routine running functions to investigate the synaptic weight structure in a network
* 'assemblyAvalancheStatistics.py' - determines the statistics of avalanche occurrence within the assemblies in a network
* 'averageFileColumnsAdvanced.py' - averages data columns across files (used to average over multiple weight traces or probability distributions, for example)
* 'averageWeights-py' - averages across multiple weight matrices
* 'calculateMIa.py' - calculates the mutual information from two firing rate distributions
* 'calculateQ.py' - calculates the pattern completion coefficient Q for an input-defined cell assembly from a firing rate distribution
* 'extractParamsQMI.py' - recursively extracts the Q and MI measures along with the simulation parameters from directories containing simulation data (intended to process many datasets to produce raster plots)
* 'frequencyAnalysisSpikeRaster.py' - computes the frequency spectrum of spike raster data
* 'meanCorrelations.py' - computes firing-rate correlations for neuron pairs from spike raster data and averages over subpopulations of the network
* 'numberOfSpikesInBins.py' - computes the distribution of spikes per time bin from cell assembly time series data
* 'overlapParadigms.py' - defines paradigms of overlapping cell assemblies
* 'utilityFunctions.py' - diverse utility functions, e.g., to read firing-rate, early- and late-phase weight data from '[timestamp]\_net\_[time].txt' files produced by the simulation program
* 'valueDistributions.py' - functions to analyze and plot weight and firing-rate distributions

