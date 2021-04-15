/***************************************************
 *** Definitions for special network simulations ***
 ***************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

// Special case to measure the plasticity between two neurons as a function of their firing rates (uses prot_learn and prot_recall to stimulate)
//#define PLASTICITY_OVER_FREQ
#ifdef PLASTICITY_OVER_FREQ
	#define TWO_NEURONS_ONE_SYNAPSE
	// here, NEURON_MODEL and SYNAPSE_MODEL options are rather irrelevant because pre- and postsynaptic neuron are modeled as Poisson neurons
#endif

// Special case of basic induction protocols for synaptic plasticity, with the same model as Li et al., 2016
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
	#define STIPULATE_CA OFF // switches off stipulation of a cell assembly
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef LTP_FR_THRESHOLD
	#define LTP_FR_THRESHOLD OFF
	#undef LTD_FR_THRESHOLD
	#define LTD_FR_THRESHOLD OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN ON
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM OFF
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
	#undef SAVE_NET_STATE
	#define SAVE_NET_STATE OFF // switches off saving the network state
#endif


// Special case of basic induction protocols for synaptic plasticity, with monoxeponential synapses and LIF model
// as in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1038/s42003-021-01778-y
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
	#define STIPULATE_CA OFF // switches off stipulation of a cell assembly
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef LTP_FR_THRESHOLD
	#define LTP_FR_THRESHOLD OFF
	#undef LTD_FR_THRESHOLD
	#define LTD_FR_THRESHOLD OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN ON
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM OFF
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
	#undef SAVE_NET_STATE
	#define SAVE_NET_STATE OFF // switches off saving the network state
#endif

// Special case for seeking mean input current
//#define SEEK_I_0 0.5 // if defined, I_const will be varied and the I_const value that leads to the defined mean frequency value (e.g., 0.5 Hz) in the absence of plasticity and stimulation will be determined
#ifdef SEEK_I_0
	#undef PLASTICITY
	#define PLASTICITY OFF // switches off plasticity
#endif

// Simulations to learn and consolidate organizational paradigms in Luboeinski and Tetzlaff, 2021, "Organization and priming of long-term memory representations with two-phase plasticity"
//#define ORGANIZATION_P2
#ifdef ORGANIZATION_P2
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE CORE_SHAPE_CMD // is to be set via compiler option, see shell script "compile_organization"
	#undef CORE_SIZE
	#define CORE_SIZE 600
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef LTP_FR_THRESHOLD
	#define LTP_FR_THRESHOLD 40
	#undef LTD_FR_THRESHOLD
	#define LTD_FR_THRESHOLD OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN ON
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM ON // important for priming
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
#endif


// Simulations to learn and consolidate organizational paradigms in Luboeinski and Tetzlaff, 2021, "Organization and priming of long-term memory representations with two-phase plasticity"
// without LTD
//#define ORGANIZATION_NOLTD_P2
#ifdef ORGANIZATION_NOLTD_P2
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE CORE_SHAPE_CMD // is to be set via compiler option, see shell script "compile_organization"
	#undef CORE_SIZE
	#define CORE_SIZE 600
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef LTP_FR_THRESHOLD
	#define LTP_FR_THRESHOLD 40
	#undef LTD_FR_THRESHOLD
	#define LTD_FR_THRESHOLD 1000 // effectively switches off LTD
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN ON
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM ON
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
#endif

// Simulations to test recall in Luboeinski and Tetzlaff, 2021, "Organization and priming of long-term memory representations with two-phase plasticity"
//#define RECALL_P2
#ifdef RECALL_P2
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY OFF // no plasticity
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE CORE_SHAPE_CMD // is to be set via compiler option, see shell script "compile_organization"
	#undef CORE_SIZE
	#define CORE_SIZE 600
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN OFF
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM OFF
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
#endif

// Simulations to investigate the spontaneous activation of assemblies in Luboeinski and Tetzlaff, 2021, "Organization and priming of long-term memory representations with two-phase plasticity"
//#define ACTIVATION_P2
#ifdef ACTIVATION_P2
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY OFF // no plasticity
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN OFF
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM OFF
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD OFF
	#undef OSCILL_INP
	#define OSCILL_INP OFF
#endif

// Special case for simulations of memory consolidation and recall with intermediate recall in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1038/s42003-021-01778-y
// (using the learning stimulation to apply an intermediate recall stimulus)
//#define INTERMEDIATE_RECALL_LT2021
#ifdef INTERMEDIATE_RECALL_P1
	#define MEMORY_CONSOLIDATION_P1
	#define CORE_SIZE_CMD 150
#endif

// Simulations of memory consolidation and recall in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1038/s42003-021-01778-y
//#define MEMORY_CONSOLIDATION_P1
#ifdef MEMORY_CONSOLIDATION_P1
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE FIRST
	#undef CORE_SIZE
	#define CORE_SIZE CORE_SIZE_CMD // can be set via compiler option, see shell script "compile_sizes"
	#undef COND_BASED_SYN
	#define COND_BASED_SYN OFF
	#undef SYN_SCALING
	#define SYN_SCALING OFF
	#undef DENDR_SPIKES
	#define DENDR_SPIKES OFF
	#undef LTP_FR_THRESHOLD
	#define LTP_FR_THRESHOLD OFF
	#undef LTD_FR_THRESHOLD
	#define LTD_FR_THRESHOLD OFF
	#undef FF_AFTER_LEARN
	#define FF_AFTER_LEARN ON
	#undef FF_AFTER_STIM
	#define FF_AFTER_STIM OFF
	#undef FF_AFTER_NETLOAD
	#define FF_AFTER_NETLOAD ON
	#undef OSCILL_INP
	#define OSCILL_INP OFF
#endif

