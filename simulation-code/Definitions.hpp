/***************************************************
 *** Definitions for network simulation settings ***
 ***************************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

// Generic definitions
#define NETWORKSIMULATION // important for Tools.cpp
#define EPSILON 10e-12 // very small number that is counted as zero
#define ON 1 // to switch on a function
#define OFF 0 // to switch off a function

// STIM_TYPE
#define POISSON_STIMULATION 1 // enables stimulation with Poisson-like spiking
#define DET_STIMULATION 2 // enables stimulation with deterministic pattern
#define GAUSS_STIMULATION 3 // enables stimulation with Gaussian current
#define OU_STIMULATION 4 // enables stimulation with Ornstein Uhlenbeck current

// NEURON_MODEL
#define LIF 1 // usage of Leaky Integrate-and-Fire model for neurons
#define MAT2 2 // usage of Multi-Adaptive Threshold model with two time constants for neurons

// NEURON TYPE
#define TYPE_EXC 1 // excitatory
#define TYPE_INH 2 // inhibitory

// SYNAPSE_MODEL
#define DELTA 1 // usage of delta synapses (synaptic time constant -> 0)
#define MONOEXP 2 // usage of synapses with monoexponentially decaying PSPs (finite synaptic time constant)

// PROTEIN POOL SETTINGS
#define POOLS_C 1 // one common pool
#define POOLS_PD 2 // one pool for potentiation and one for depression
#define POOLS_PCD 3 // one pool for potentiation, one for depression and one common pool

// SPIKE_PLOTTING
#define NUMBER 1 // plotting the spike number over time bins
#define RASTER 2 // plotting the spikes of all neurons
#define NUMBER_AND_RASTER 3 // plotting spike number and spike raster

// CONNECTION_PLOT
#define CREATE_SCRIPT 1 // creates a gnuplot file for a plot showing the detailed interneuronal connections within the excitatory network
#define CREATE_PLOT 2 // creates the gnuplot file

// CELL ASSEMBLIES
#define FIRST 1 // simply use the first block of neurons in the network as the assembly (would equal a hypothetical "OVERLAP100")
#define SECOND 2 // simply use the second distinct block of neurons in the network as the assembly (would equal a hypothetical "OVERLAP0")
#define OVERLAP10_2ND 3 // use a second block of neurons as the assembly, overlapping by 10% with the first assembly
#define OVERLAP15_2ND 4 // use a second block of neurons as the assembly, overlapping by 15% with the first assembly
#define OVERLAP20_2ND 5 // use a second block of neurons as the assembly, overlapping by 20% with the first assembly
#define THIRD 6 // simply use the third distinct block of neurons in the network as the assembly
#define OVERLAP10_3RD 7 // use a third block of neurons as the assembly, overlapping by 5% with the first assembly exclusively, by 5% with the second assembly exclusively, and by 5% with both
#define OVERLAP10_3RD_NO_ABC 8 // use a third block of neurons as the assembly, overlapping by 10% with the first assembly exclusively, and by 10% with the second assembly exclusively
#define OVERLAP10_3RD_NO_AC_NO_ABC 9 // use a third block of neurons as the assembly, overlapping by 10% with the first assembly exclusively
#define OVERLAP10_3RD_NO_BC_NO_ABC 10 // use a third block of neurons as the assembly, overlapping by 10% with the second assembly exclusively
#define OVERLAP15_3RD 11 // use a third block of neurons as the assembly, overlapping by 7.5% with the first assembly exclusively, by 7.5% with the second assembly exclusively, and by 7.5% with both
#define OVERLAP15_3RD_NO_ABC 12 // use a third block of neurons as the assembly, overlapping by 15% with the first assembly exclusively, and by 15% with the second assembly exclusively
#define OVERLAP15_3RD_NO_AC_NO_ABC 13 // use a third block of neurons as the assembly, overlapping by 15% with the first assembly exclusively
#define OVERLAP15_3RD_NO_BC_NO_ABC 14 // use a third block of neurons as the assembly, overlapping by 15% with the second assembly exclusively
#define OVERLAP20_3RD 15 // use a third block of neurons as the assembly, overlapping by 10% with the first assembly exclusively, by 10% with the second assembly exclusively, and by 10% with both
#define OVERLAP20_3RD_NO_ABC 16 // use a third block of neurons as the assembly, overlapping by 20% with the first assembly exclusively, and by 20% with the second assembly exclusively
#define OVERLAP20_3RD_NO_AC_NO_ABC 17 // use a third block of neurons as the assembly, overlapping by 20% with the first assembly exclusively
#define OVERLAP20_3RD_NO_BC_NO_ABC 18 // use a third block of neurons as the assembly, overlapping by 20% with the second assembly exclusively
#define RAND 19 // use randomly selected neurons as the assembly

// PLASTICITY
#define CALCIUM 1 // use the calcium model as plasticity mechanism (early phase only)
#define CALCIUM_AND_STC 2 // use the calcium model with synaptic and capture as plasticity mechanism
#define STDP 3 // use an STDP rule as plasticity mechanism (early phase only)
#define STDP_AND_STC 4 // use an STDP rule with synaptic and capture as plasticity mechanism
