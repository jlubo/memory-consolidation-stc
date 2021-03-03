# Memory consolidation in recurrent spiking neural networks


## Outline

This package serves to simulate recurrent spiking neural networks with calcium-based synaptic plasticity and synaptic tagging and capture. 
The C++ code for these simulations is located in the directory 'simulation-code/'. 

Pre-compiled binaries of the simulation code for Linux (tested with Kernel version 5.4) can be found in 'simulation-bin/', compiled and linked with g++ in version 7.4.0 and boost in version 1.65.1.

Furthermore, the package contains Python scripts to analyze the data produced by the simulations. These scripts are located in the directory 'analysis/'. 

If you make use of the code or the binaries, please cite our paper:

	Luboeinski, J., Tetzlaff, C. Memory consolidation and improvement by synaptic tagging and capture in recurrent neural networks. Commun. Biol. 4, 275 (2021). 
	https://doi.org/10.1038/s42003-021-01778-y
	
The paper presents the model that underlies the simulation code provided here, as well as findings derived from the model.
However, the code contains additional features that have not been used in publications yet. Please feel free
to contact us on questions about the code and on further investigations that can be done with it.


## Simulation code

### Files

* 'NetworkMain.cpp' - main function that initializes network simulations
* 'NetworkSimulation.cpp' - class performing network simulations
* 'Network.cpp' - class describing the network
* 'Neuron.cpp' - class describing one neuron
* 'Stimulus.hpp' - class describing a stimulus
* 'StimulusProtocols.hpp' - class to define specific stimulus protocols
* 'Definitions.hpp' - general definitions
* 'SpecialCases.hpp' - definitions for special simulations
* 'Tools.hpp' - collection of utility functions
* 'Plots.hpp' - collection of plotting functions employing *gnuplot*
* 'plotFunctions.py' - collection of plotting functions employing *Matplotlib*

### Compiling and linking

The simulation code comes with shell scripts to compile and link the simulation code for different purposes:
* 'compile' (has the same effect as the included Makefile) - compiles the code for performing network simulations to learn, consolidate and recall a memory representation
* 'compile_2N1S' - compiles the code for applying basic plasticity protocols to a single synapse
* 'compile_irs' - compiles the code for performing network simulations to learn and consolidate a memory representation, apply intermediate stimulation and recall

### Running the simulation

The simulation is run by executing the binary with or without command line options as defined in 'NetworkMain.cpp' (e.g., via one of the following shell scripts).
Please note that there are other options that have to be set in 'NetworkSimulation.cpp' before compiling and cannot be changed during
runtime.

The directory 'simulation-bin/' contains the following sample shell scripts:
* 'run' - learn a memory representation, let it consolidate, and recall after 8 hours
* 'run2' - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, let the memory representation consolidate, and recall after 8 hours
* 'run3' - learn a memory representation, save the network state, and recall after 10 seconds; load the network state, apply intermediate stimulation, let the memory representation consolidate, and recall after 8 hours
* 'run_2N1S' - reproduce single-synapse data resulting from basic induction protocols for synaptic plasticity

The file 'connections.txt' contains the default connectivity matrix used in Luboeinski and Tetzlaff, Commun. Biol., 2020. If this file is absent, the simulation program will automatically generate a new network structure.


## Analysis scripts

The following scripts serve to process and analyze the data produced by the simulation code. They were tested to run with Python 3.7.3 and NumPy 1.16.4. Note that some of the script files depend on others.

### Files

* 'adjacencyFunctions.py' - functions to analyze the connectivity and weights (used to compute mean and standard deviation of early- and late-phase weights)
* 'averageFileColumnsAdvanced.py' - averages data columns across files (used to average over multiple weight traces)
* 'averageWeights-py' - averages across multiple weight matrices
* 'calculateMIa.py' - calculates the mutual information from two firing rate distributions
* 'calculateQ.py' - calculates the pattern completion coefficient Q for an input-defined cell assembly from a firing rate distribution
* 'extractParamsQMI.py' - recursively extracts the Q and MI measures along with the simulation parameters from directories containing simulation data
* 'meanCorrelations.py' - computes firing rate correlations for neuron pairs from spike raster data and averages over subpopulations of the network
* 'probDistributions.py' - functions to analyze and plot weight and activity distributions
* 'readWeightData.py' - reads the firing rate data and early- and late-phase weight matrices from '[timestamp]\_net\_[time].txt' files produced by the simulation program

