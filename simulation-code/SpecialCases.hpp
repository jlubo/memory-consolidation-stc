/***************************************************
 *** Definitions for special network simulations ***
 ***************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

// Special case to basic induction protocols for synaptic plasticity, with the same model as Li et al., 2016
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE
#ifdef TWO_NEURONS_ONE_SYNAPSE
	#undef STIM_TYPE
	#define STIM_TYPE POISSON_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL MAT2 // uses Multi-Adaptive Threshold Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL DELTA // uses delta synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C // uses protein pool setting C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF // switches off stipulation
	#undef SAVE_NET_STATE
	#define SAVE_NET_STATE OFF // switches off saving the network state
#endif


// Special case to basic induction protocols for synaptic plasticity, with monoxeponential synapses and LIF model
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE_ALT
#ifdef TWO_NEURONS_ONE_SYNAPSE_ALT
	#define TWO_NEURONS_ONE_SYNAPSE
	#undef STIM_TYPE
	#define STIM_TYPE POISSON_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF // uses Leaky Integrate-and-Fire Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP // uses monoexponential synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C // uses protein pool setting C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF // switches off stipulation
	#undef SAVE_NET_STATE
	#define SAVE_NET_STATE OFF // switches off saving the network state
#endif


// Special case to measure the weight between two neurons as a function of their firing rates (uses prot_learn and prot_recall to stimulate)
//#define PLASTICITY_OVER_FREQ
#ifdef PLASTICITY_OVER_FREQ
	#define TWO_NEURONS_ONE_SYNAPSE
	#undef STIM_TYPE
	#define STIM_TYPE POISSON_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL MAT2 // uses Multi-Adaptive Threshold Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP // uses delta synapses
	#undef PLASTICITY
	#define PLASTICITY OFF // switches on plasticity
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF // switches off stipulation
#endif

// Special case for seeking mean input current
//#define SEEK_I_0 0.5 // if defined, I_const will be varied and the I_const value that leads to the defined mean frequency value (e.g., 0.5 Hz) in the absence of plasticity and stimulation will be determined
#ifdef SEEK_I_0
	#undef PLASTICITY
	#define PLASTICITY OFF // switches off plasticity
#endif

// Special case for using the learning stimulation to apply an intermediate recall stimulus
//#ifdef INTERMEDIATE_RECALL
