/***************************************************
 *** Definitions for special network simulations ***
 ***************************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

// Special case of a network stimulated with a single pulse to evoke one spike in neuron 6
//#define ONESPIKE_EXC
#ifdef ONESPIKE_EXC
	#define ONESPIKE
	#warning "Special case: ONESPIKE_EXC"
#endif

// Special case of a network stimulated with a single pulse to evoke one spike in neuron 1615
//#define ONESPIKE_INH
#ifdef ONESPIKE_INH
	#define ONESPIKE
	#warning "Special case: ONESPIKE_INH"
#endif

// Special case of a network stimulated with a single pulse to evoke one spike in one neuron (similar configuration as MEMORY_CONSOLIDATION_P1)
#ifdef ONESPIKE
	#undef STIM_TYPE
	#define STIM_TYPE DET_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE FIRST
	#undef CORE_SIZE
	#define CORE_SIZE 1
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
	#warning "Special case: ONESPIKE"
#endif

// Special case of a network that has one neuron spiking at maximal activity (similar configuration as MEMORY_CONSOLIDATION_P1; also shares some features with ONESPIKE_EXC)
//#define MAX_ACTIVITY_NEURON
#ifdef MAX_ACTIVITY_NEURON
	#undef STIM_TYPE
	#define STIM_TYPE DET_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
	#undef PROTEIN_POOLS
	#define PROTEIN_POOLS POOLS_C
	#undef STIPULATE_CA
	#define STIPULATE_CA OFF
	#undef CORE_SHAPE
	#define CORE_SHAPE FIRST
	#undef CORE_SIZE
	#define CORE_SIZE 1
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
	#warning "Special case: MAX_ACTIVITY_NEURON"
#endif

// Special case to measure the plasticity between two neurons as a function of their firing rates (uses prot_learn and prot_recall to stimulate)
// as in Luboeinski, 2021, doctoral thesis
//#define PLASTICITY_OVER_FREQ
#ifdef PLASTICITY_OVER_FREQ
	#define TWO_NEURONS_ONE_SYNAPSE // sets general flag of single synapse simulations
	// here, NEURON_MODEL and SYNAPSE_MODEL options are rather irrelevant because pre- and postsynaptic neuron are modeled as Poisson neurons
	#warning "Special case: PLASTICITY_OVER_FREQ"
#endif

// Special case of basic induction protocols for synaptic plasticity, with monoxeponential synapses and LIF model
// as in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1038/s42003-021-01778-y
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE_P1
#ifdef TWO_NEURONS_ONE_SYNAPSE_P1
	#define TWO_NEURONS_ONE_SYNAPSE // sets general flag of single synapse simulations
	#undef STIM_TYPE
	#define STIM_TYPE POISSON_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF // uses Leaky Integrate-and-Fire Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP // uses monoexponential synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: TWO_NEURONS_ONE_SYNAPSE_P1"
#endif



// Special case of basic induction protocols for synaptic plasticity, with the same model as Li et al., 2016, https://doi.org/10.1371/journal.pone.0161679
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE_LI2016
#ifdef TWO_NEURONS_ONE_SYNAPSE_LI2016
	#define TWO_NEURONS_ONE_SYNAPSE // sets general flag of single synapse simulations
	#undef STIM_TYPE
	#define STIM_TYPE POISSON_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL MAT2 // uses Multi-Adaptive Threshold Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL DELTA // uses delta synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: TWO_NEURONS_ONE_SYNAPSE_LI2016"
#endif

// Special case of minimal example for early-phase plasticity, induced by a few pre-defined spikes
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE_MIN
#ifdef TWO_NEURONS_ONE_SYNAPSE_MIN
	#define TWO_NEURONS_ONE_SYNAPSE // sets general flag of single synapse simulations
	#undef STIM_TYPE
	#define STIM_TYPE DET_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF // uses Leaky Integrate-and-Fire Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP // uses monoexponential synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: TWO_NEURONS_ONE_SYNAPSE_MIN"
#endif

// Special case of neuron that is first stimulated with a constant current such that its membrane converges to a depolarized state, and then 
// relaxes back to the reversal potential 
// --> neuron 0 is stimulated via neuron 1, which is depolarized following the stimulus protocol
// --> changes some global variables in int main()
//#define TWO_NEURONS_ONE_SYNAPSE_CONV
#ifdef TWO_NEURONS_ONE_SYNAPSE_CONV
	#define TWO_NEURONS_ONE_SYNAPSE // sets general flag of single synapse simulations
	#undef STIM_TYPE
	#define STIM_TYPE DET_STIMULATION // uses Poisson-like spikes
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF // uses Leaky Integrate-and-Fire Model
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP // uses monoexponential synapses
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC // switches on plasticity
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: TWO_NEURONS_ONE_SYNAPSE_CONV"
#endif

// Special case for seeking mean input current
//#define SEEK_I_0 0.5 // if defined, I_const will be varied and the I_const value that leads to the defined mean frequency value (e.g., 0.5 Hz) in the absence of plasticity and stimulation will be determined
#ifdef SEEK_I_0
	#undef PLASTICITY
	#define PLASTICITY OFF // switches off plasticity
	#warning "Special case: SEEK_I_0"
#endif

// Simulations of consolidation of (topologically) spatial and temporal patterns with neuromodulator-dependent STC in Lehr, Luboeinski, Tetzlaff, 2022
//#define MEMORY_CONSOLIDATION_NM_STC_P3
#ifdef MEMORY_CONSOLIDATION_NM_STC_P3
	#define MEMORY_CONSOLIDATION_P1
	#define CORE_SIZE_CMD 150
	#warning "Special case: MEMORY_CONSOLIDATION_NM_STC_P3"
#endif

// Simulations to learn and consolidate organizational paradigms of attractor memories as in Luboeinski, 2021, doctoral thesis
//#define ORGANIZATION_ATTR
#ifdef ORGANIZATION_ATTR
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#define LTP_FR_THRESHOLD 80
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
	#warning "Special case: ORGANIZATION_ATTR"
#endif

// Simulations to investigate the activation of assemblies by oscillating inhibition (at 1 Hz, if dt=0.2s) as in Luboeinski, 2021, doctoral thesis
//#define ACTIVATION_OSC_1HZ
#ifdef ACTIVATION_OSC_1HZ
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
	#define OSCILL_INP 5000 // timesteps per oscillation period
	#warning "Special case: ACTIVATION_OSC_1HZ"
#endif


// Simulations to investigate the activation of assemblies by oscillating inhibition (at 5 Hz, if dt=0.2s) as in Luboeinski, 2021, doctoral thesis
//#define ACTIVATION_OSC_5HZ
#ifdef ACTIVATION_OSC_5HZ
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
	#define OSCILL_INP 1000 // timesteps per oscillation period
	#warning "Special case: ACTIVATION_OSC_5HZ"
#endif

// Simulations to learn and consolidate organizational paradigms in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1101/2021.04.15.439982
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
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: ORGANIZATION_P2"
#endif

// Simulations to learn and consolidate organizational paradigms in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1101/2021.04.15.439982
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
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: ORGANIZATION_NOLTD_P2"
#endif

// Simulations to learn and consolidate organizational paradigms in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1101/2021.04.15.439982
// with randomly drawn initial weights
//#define ORGANIZATION_RANDW_P2
#ifdef ORGANIZATION_RANDW_P2
	#undef STIM_TYPE
	#define STIM_TYPE OU_STIMULATION
	#undef NEURON_MODEL
	#define NEURON_MODEL LIF
	#undef SYNAPSE_MODEL
	#define SYNAPSE_MODEL MONOEXP
	#undef PLASTICITY
	#define PLASTICITY CALCIUM_AND_STC
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS ON
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
	#warning "Special case: ORGANIZATION_RANDW_P2"
#endif

// Simulations to test recall in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1101/2021.04.15.439982
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
	#warning "Special case: RECALL_P2"
#endif

// Simulations to investigate the spontaneous activation of assemblies in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1101/2021.04.15.439982 (also for attractor assemblies in Luboeinski, 2021, doctoral thesis)
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
	#warning "Special case: ACTIVATION_P2"
#endif

// Special case for simulations of memory consolidation and recall with intermediate recall in Luboeinski and Tetzlaff, 2021, https://doi.org/10.1038/s42003-021-01778-y
// (using the learning stimulation to apply an intermediate recall stimulus)
//#define INTERMEDIATE_RECALL_P1
#ifdef INTERMEDIATE_RECALL_P1
	#define MEMORY_CONSOLIDATION_P1
	#define CORE_SIZE_CMD 150
	#warning "Special case: INTERMEDIATE_RECALL_P1"
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
	#undef RAND_INIT_WEIGHTS
	#define RAND_INIT_WEIGHTS OFF
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
	#warning "Special case: MEMORY_CONSOLIDATION_P1"
#endif

