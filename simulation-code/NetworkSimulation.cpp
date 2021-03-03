/*************************************************************************
 *** Simulation of a spiking neural network with early- and late-phase ***
 ***                         long-term plasticity                      ***
 *************************************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

#include <iostream>
#include <fstream>
#include <unistd.h> // for chdir()
#include "Definitions.hpp"

using namespace std;

// Simulation options (cf. Definitions.hpp)
#define STIM_TYPE        OU_STIMULATION // the type of stimulation that enters the neurons (not to be confused with the stimulus protocol)
#define STIM_PREPROC     OFF // defines if stimuli are pre-processed (at cost of memory, can speed up, but can also slow down runtime!)
#define NEURON_MODEL     LIF // definition of the neuron model
#define SYNAPSE_MODEL    MONOEXP // the synapse model that is used
#define PLASTICITY       CALCIUM_AND_STC // defines what type of plasticity shall be used (or OFF)
#define PROTEIN_POOLS    POOLS_C // the protein pool setting for synaptic consolidation/plasticity-related proteins
#define STIPULATE_CA     OFF // if ON: stipulate a cell assembly with strong interconnections at the beginning of the learning stimulus, no actual learning is required
#define CORE_SIZE        150 // size of the cell assembly
#define COND_BASED_SYN   OFF // if ON: use conductance-based synapses that may be different for inhibitory and excitatory neurons
#define SYN_SCALING      OFF // if ON: use synaptic scaling following Tetzlaff et al., 2013
#define DENDR_SPIKES     OFF // if ON: use dendritic spikes as described by Jahnke et al., 2015
#define FF_AFTER_LEARN   ON // if ON: employ fast-forward computing after learning stimulus
#define FF_AFTER_STIM    OFF // if ON: employ fast-forward computing as soon as all stimulation has ended
#define FF_AFTER_NETLOAD ON // if ON: employ fast-forward computing after loading a previous network state

// Output options (cf. Definitions.hpp)
#define SPIKE_PLOTTING           NUMBER_AND_RASTER // defines what information about spiking dynamics is saved and plotted
#define PRINT_CONNECTIONS        ON // if ON: output of connectivity matrix of excitatory subnetwork using Network::printConnections()
#define NET_OUTPUT_PRECISION     7 // precision of numbers (also of h_0) used for network plots
#define OUTPUT_PRECISION         9 // precision of numbers used for plots over time
#define SAVE_NET_STATE           ON // if ON: output of whole simulation data before recall

#include "SpecialCases.hpp"
#include "Tools.hpp"
#include "Network.cpp"
#include "Plots.hpp"
#include "StimulusProtocols.hpp"

/*** NetworkSimulation class ***
 * simulates a network of neurons, has instance of Network class */
class NetworkSimulation {

private:

/*** Simulation parameters ***/
const int Nl; // number of neurons in one line (row or column) of the excitatory population
const int Nl_inh; // number of neurons in one line (row or column) of the inhibitory population
const double dt; // s, one time step for numerical integration
const double t_max;  // s, total duration of simulation
const double pc; // connection probability for unidirectional neuron connections
double tau_syn; // s, synaptic time constant
double w_ei; // E->I coupling strength in units of h_0
double w_ie; // I->E coupling strength in units of h_0
double w_ii; // I->I coupling strength in units of h_0
const double t_wfr; // s, size of the time window for computing instantaneous firing rates
const int wfr; // size of the time window for computing instantaneous firing rates in time steps
string prot_learn; // the stimulation protocol for training
string prot_recall; // the stimulation protocol for recall
double recall_fraction; // recall stimulus is applied to this fraction of the original assembly
double w_stim; // nA s, synaptic weight for incoming stimulation to neurons
int N_stim; // the number of hypothetical synapses per neuron that are used for stimulation (each of strength w_stim)
bool ff_enabled; // specifies if fast-forward mode can be used
Network net;  // the network
double z_max;  // maximum late-phase coupling strength

/*** Output parameters ***/
vector<int> exc_neuron_output {};
vector<int> inh_neuron_output {};
vector<synapse> synapse_output {};
int output_period; // number of timesteps to pass for the next data output (if set to 1, most detailed output is obtained)
int net_output_period; // number of timesteps to pass for the next network plot
vector<int> net_output {}; // vector of times for selected network plots

#ifdef SEEK_I_0
double *seekic; // pointer to a variable to communicate with the main(...) function while seeking I_0
#endif

/*** saveParams ***
 * Saves the crucial parameters in a file *
 * - str: string containing additional information */
void saveParams(string str)
{
	ofstream f (dateStr("_PARAMS.txt"));
	f << dateStr() << endl; // time stamp
	f << endl;

	// Parameters
	f << "Simulation parameters:" << endl;
	f << "dt = " << dt << " s" << endl;
	f << "t_max = " << t_max << " s (" << int(ceil(t_max / dt)) << " steps)" << endl;
	f << "t_wfr = " << t_wfr << " s" << endl;
	f << "stimulation type = "
#if STIM_TYPE == POISSON_STIMULATION
	  << "Poisson";
#elif STIM_TYPE == DET_STIMULATION
	  << "deterministic";
#elif STIM_TYPE == OU_STIMULATION
	  << "OU";
#elif STIM_TYPE == GAUSS_STIMULATION
	  << "Gauss";
#endif
	f << endl;
#if STIPULATE_CA == OFF
	f << "learning stimulus = " << prot_learn << endl;
#else
	f << "learning stimulus = STIP" << endl;
#endif
	f << "recall stimulus = " << prot_recall << endl;
	f << "core = first " << CORE_SIZE << " neurons" << endl;
	f << "recall fraction = " << recall_fraction << endl;
	f << "N_stim = " << N_stim << endl;

	net.saveNetworkParams(&f);
	f << endl;

	// Additional information
	f << "Purpose: " << str << endl;
	f.close();
}

/*** addToParamsFile ***
 * Adds a text to the parameter text file *
 * - str: string containing information to be added, like the elapsed time */
void addToParamsFile(string str)
{
	ofstream f (dateStr("_PARAMS.txt"), ofstream::out | ofstream::app);
	f << str << endl;
	f.close();
}

/*** instFiringRates ***
 * Computes the instantaneous firing rates of all the network neurons and prints them to a given data file *
 * - txt_net_tprime: pointer to the data file *
 * - jprime: timestep for which the firing rates shall be calculated */
void instFiringRates(ofstream* txt_net_tprime, int jprime)
{
	double tprime = jprime*dt;
	double t_wfr_eff;

	if (jprime < wfr/2.)
		t_wfr_eff = t_wfr/2. + tprime;
	else
		t_wfr_eff = t_wfr;

	for (int m=0; m < pow2(Nl); m++)
	{
		int num_spikes = 0; // number of spikes in time window t_wfr_eff
		bool removed = false; // specified if old, now irrelevant spikes have been removed
		int sp = 1;
		int hist_size = net.getSpikeHistorySize(m);

		while ( sp <= hist_size )
		{
			int spt = net.getSpikeTime(sp, m);

			if (spt >= jprime-wfr/2.)
			{
				if (!removed) // is entered for first spike after jprime-wfr/2.
				              // for removal, (net_output_period + wfr/2) has to be larger than t_syn_delay_steps and t_Ca_delay_steps (which is usually the case)
				{
					net.removeSpikes(1, sp-1, m);
					sp = 0;
					hist_size = net.getSpikeHistorySize(m);
					net.resetLastSpikeIndex(m);
					removed = true;
				}
				if (spt < jprime+wfr/2.)
					num_spikes++;
				else // beyond time window
					break;
			}

			sp++;
		}

		double rate = num_spikes / t_wfr_eff;
		*txt_net_tprime << fixed << rate;

		if ((m+1) % Nl != 0) // still in the same neuron row
			*txt_net_tprime << "\t\t";
		else
			*txt_net_tprime << endl; // next neuron row begins

	}

	// Plot everything
	*txt_net_tprime << endl;
	txt_net_tprime->close();
	delete txt_net_tprime;
}

/*** netOutput ***
 * Returns true if a given timestep has been selected for network output *
 * - j: time given in timesteps *
 * - return: true if j is in net_output vector */
bool netOutput(int j)
{
	for (int i=0; i<net_output.size(); i++)
	{
		if (net_output[i] == j)
			return true;
	}
	return false;
}

public:


#ifdef TWO_NEURONS_ONE_SYNAPSE
/*** getPlasticityType ***
 * Returns the kind of plasticity evoked by the stimulus *
 * - return: 0 for ELTP, 1 for ELTP with tag, 2 for LLTP, 3 for ELTD, 4 for ELTD with tag, 5 for LLTD, -1 else */
int getPlasticityType()
{
	return net.getPlasticityType();
}

/*** getMaxDev ***
 * Returns the maximum deviation from h_0 *
 * - return: the maximum deviation from h_0*/
double getMaxDev()
{
	return net.getMaxDev();
}
#endif

/*** simulate ***
 * Runs the network simulation *
 * - working_dir: working directory *
 * - first_sim: tells if this is the first simulation run by the main(...) function *
 * - _purpose: a short text describing the purpose of this simulation */
int simulate(string working_dir, bool first_sim, string _purpose)
{
// ==============================================================================================================================
	// Initialize

	// Start time measurement
	timeMeasure(true);

	// Constants
	const int n = int(ceil(t_max / dt)); // number of time steps
	const int tenth_sec = int(round(0.1/dt)); // a tenth of asecond in timesteps
	const double h_0 = net.getInitialWeight(); w_stim = h_0; // initial synaptic weight; set coupling strength for stimulation

	const string separator = getSeparator(); cout << separator << endl; // string of characters for a separator in command line

	// Neuronal and synaptic output
#ifdef TWO_NEURONS_ONE_SYNAPSE
	exc_neuron_output = vector<int> {0,1};
	inh_neuron_output = vector<int> {};
	synapse_output = vector<synapse> {synapse(1,0)};
#elif defined SEEK_I_0
	exc_neuron_output = vector<int> {};
	inh_neuron_output = vector<int> {};
#else
	if (Nl == 20)
	{
		exc_neuron_output = vector<int> {0,1,2};
		inh_neuron_output = vector<int> {400};
		synapse_output = vector<synapse> {synapse(0,1),synapse(0,50),synapse(0,100),synapse(400,0)};
	}
	else if (Nl == 35)
	{
		//exc_neuron_output = vector<int> {608,609,610};
		exc_neuron_output = vector<int> {0,1,2};
		inh_neuron_output = vector<int> {1225};
		//synapse_output = vector<synapse> {synapse(608,609),synapse(609,608),synapse(609,610)};
		synapse_output = vector<synapse> {synapse(0,1),synapse(0,50),synapse(0,100)};
	}
	else if (Nl == 40)
	{
		//exc_neuron_output = vector<int> {816,817,818};
		//exc_neuron_output = vector<int> {cNN(20,21),cNN(23,21),cNN(30,21)}; // three neurons: one "as", one "ans" and one "e"
		//exc_neuron_output = vector<int> {1, 53, 260, 660, 777, 940, 941};
		exc_neuron_output = vector<int> {6, 68};
		//inh_neuron_output = vector<int> {1615, 1710, 1750};
		inh_neuron_output = vector<int> {1615};
		//synapse_output = vector<synapse> {synapse(777,260), synapse(777,940), synapse(777,941),synapse(660,260),
		//                                  synapse(660,940), synapse(660,941), synapse(782,1),  synapse(782,53)};
		synapse_output = vector<synapse> {synapse(6,68)};
	}
	else if (Nl == 50)
	{
		exc_neuron_output = vector<int> {1,640,1300};
		inh_neuron_output = vector<int> {2500};
		synapse_output = vector<synapse> {synapse(1,9),synapse(1,640),synapse(1,1300)};
	}
#endif

#ifndef TWO_NEURONS_ONE_SYNAPSE
	// Read connectivity matrix from file
	if (copyFile("../connections.txt", "connections.txt"))
	{
		int r = net.readConnections(string("connections.txt"));
		if (r == 1)
			cout << "Importing connections failed - file could not be opened." << endl;
		else if (r == 0)
			cout << "Importing connections failed - dimension mismatch." << endl;
		else
			cout << "Connections successfully imported." << endl;
		remove("connections.txt");
	}

	// Read excitatory->excitatory coupling strengths from file
	if (copyFile("../coupling_strengths.txt", "coupling_strengths.txt"))
	{
		int r = net.readCouplingStrengths(string("coupling_strengths.txt"));
		if (r == 1)
			cout << "Importing coupling strengths failed - file could not be opened." << endl;
		else if (r == 0)
			cout << "Importing coupling strengths failed - dimension mismatch." << endl;
		else
			cout << "Coupling strengths successfully imported." << endl;
		remove("coupling_strengths.txt");
	}

	// Load state of the whole network
	int tb_start = net.loadNetworkState("../saved_state.txt") + 1;

#ifdef SEEK_I_0
	// Create working directory for this I_0 value
	const string path = working_dir + "/I_0 " + dtos(net.getConstCurrent(1,1), 6); // path to working directory
	system(concat("mkdir -p \"", path + "\"").c_str());

	// Try to change directory
	if (chdir(path.c_str()) == -1) {
		showChDirErrMessage();
		cout << separator << endl;
		return -1;
	}
	// Copy Python plot script
	if (!copyFile("../plotFunctions.py", "plotFunctions.py"))
	{
		cout << "Python script 'plotFunctions.py' not found!" << endl;
		return -2;
	}
#endif

	// Create directory for network plots and copy Python plot script
	system("mkdir -p \"network_plots\"");

#else // if TWO_NEURONS_ONE_SYNAPSE is defined

	int tb_start = 0;

#endif
	// Output with general information
	// if this is the first simulation run by the main(...) function (first_sim = true), use time stamp from NetworkMain.cpp, else, set a new time stamp
	cout << "\x1b[33mNetwork simulation with N_exc = " << pow2(Nl) << ", N_inh = " << pow2(Nl_inh)
		  << ", t_max = " << t_max << " s (" << dateStr("", !first_sim) << ")\x1b[0m" << endl;
	cout << "Learning protocol: " << prot_learn << endl;
	cout << "Recall protocol: " << prot_recall << endl;
	cout << "Description: " << _purpose << endl; //rectangle, \x1b[34mt_pulse = " << t_pulse << " s, \x1b[35mfrequency = " << frequency << " Hz\x1b[0m" << endl;
	cout << "Connectivity: \x1b[32mp_c = " << pc << ", \x1b[31mtau_syn = "
#if SYNAPSE_MODEL == DELTA
		  << 0
#elif SYNAPSE_MODEL == MONOEXP
		  << tau_syn*1000.
#endif
		  << " ms, \x1b[36mw_ei = " << w_ei << ", w_ie = " << w_ie << ", w_ii = " << w_ii << "\x1b[0m" << endl;
	cout << "Other parameters: \x1b[35mI_0 = " << net.getConstCurrent(1,1) << " nA, sigma_WN = " << net.getSigma(1,1) << " nA s^1/2, \x1b[37mN_stim = " << N_stim << "\x1b[0m" << endl;
	saveParams(_purpose); // create parameter file

	// Declarations and initializations
	double max_firing_rate = 0.0; // contains the maximum firing rate over all exc. neurons and stimulus amplitudes
	double min_firing_rate = numeric_limits<double>::max(); // contains the minimum firing rate over all exc. neurons and stimulus amplitudes
	double max_firing_rate_inh = 0.0; // contains the maximum firing rate over all inh. neurons and stimulus amplitudes
	double min_firing_rate_inh = numeric_limits<double>::max(); // contains the minimum firing rate over all inh. neurons and stimulus amplitudes
	ofstream txt_fr(dateStr("_fr_exc.txt")); // data for creating plot of mean neuronal firing rates in excitatory population
	ofstream txt_fr_inh(dateStr("_fr_inh.txt")); // data for creating plot of mean neuronal firing rates in inhibitory population
	//TODO: create gpl files in Plots.hpp
	ofstream gpl_fr("fr_exc_map.gpl"); // gnuplot file for map plot of mean neuronal firing rates in excitatory population
	ofstream gpl_fr_inh("fr_inh_map.gpl"); // gnuplot file for map plot of mean neuronal firing rates in inhibitory population
	ofstream gpl_cn("inc_connection_map.gpl"); // gnuplot file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream txt_cn(dateStr("_inc_connection.txt")); // data file for creating a color plot of the numbers of incoming interneuronal connections
	ofstream gplscript("gpl"); // shell script for calling all the gnuplot scripts (to re-generate the plots, e.g., for another gnuplot version)
#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
	ofstream txt_spike_number(dateStr("_spike_number.txt")); // data file containing the total number of spikes in the last output_period bins over time
#endif
#if SPIKE_PLOTTING == RASTER || SPIKE_PLOTTING == NUMBER_AND_RASTER
	ofstream txt_spike_raster(dateStr("_spike_raster.txt")); // data file containing spike times of all neurons
#endif
	ofstream txt_data(dateStr("_data.txt")); // data file containing V, h, z, Ca etc. of selected neurons and synapses
	ofstream txt_mean(dateStr("_mean_weight.txt")); // data file containing mean weights and weight standard deviations of network subpopulations
	ofstream* txt_net_t; // pointer to a data file containing all neuronal firing rates and synaptic weights at time t
	ofstream* txt_net_tprime; // pointer to a data file containing all neuronal firing rates and synaptic weights at time tprime

	ofstream logf("log.txt"); // log file containing information about data processing

	bool data_to_plot = false; // indicates if there is still network data that has to be plotted
	bool STC = true; // indicates if late-phase dynamics occur

	writePalViridis(); // create palette file for gnuplot color plots

	// learning stimulation
	Stimulus st_learn = createStimulusFromProtocols(prot_learn, "", dt, w_stim, N_stim, tau_syn); // create Stimulus object containing learning stimulation only
	// recall stimulation
	Stimulus st_recall = createStimulusFromProtocols("", prot_recall, dt, w_stim, N_stim, tau_syn); // create Stimulus object containing recall stimulation only

#if STIPULATE_CA == ON
	Stimulus st_full = st_recall; // only stipulated cell assembly: recall stimulation is all (effective) stimulation there is
#elif defined TWO_NEURONS_ONE_SYNAPSE
	Stimulus st_full = st_learn; // no network, no recall: just a dummy for st_full has to be defined
#else
	// learning + recall stimulation
	Stimulus st_full = createStimulusFromProtocols(prot_learn, prot_recall, dt, w_stim, N_stim, tau_syn); // create Stimulus object containing learning and recall stimulation
#endif
	int tb_stim_start = st_learn.getStimulationStart(); // time bin in which all stimulation begins
	int tb_stim_end = st_full.getStimulationEnd(); // time bin in which all stimulation ends
	int tb_stim_learn_end = st_learn.getStimulationEnd(); // time bin in which learning stimulation ends
	int tb_recall_start = st_recall.getStimulationStart(); // time bin in which recall stimulation begins
	net.setStimulationEnd(tb_stim_end);

	if (tb_start > 0)
		logf << "Network state loaded, starting simulation from t = " << dtos(tb_start*dt, 3) << " s" << endl;

	if (st_full.getStimulationDuration() > 0) // stimulation information in logfile
	{
		logf << "Stimulation occurs between ";

		if (st_learn.getStimulationDuration() > 0) // there is a learning stimulus
			logf << dtos(tb_stim_start*dt, 3);
		else // there is no learning stimulus
			logf << dtos(tb_recall_start*dt, 3);

		logf << " s and " << dtos(tb_stim_end*dt, 3) << " s" << endl;
	}
	else
		logf << "No stimulation" << endl;

#if STIPULATE_CA == ON
	st_learn = Stimulus(dt); // unset learning stimulus
#endif

#ifdef PLASTICITY_OVER_FREQ
	Stimulus st_learn2 = createStimulusFromProtocols(prot_recall, "", dt, w_stim, N_stim, tau_syn); // create Stimulus object for stimulating second neuron
#endif

	double p = 0.0; // percentage of process completeness

	int jprime = 0; // time step in which the last network plot data were collected

#if PRINT_CONNECTIONS == ON
	net.printConnections(dateStr("_connections.txt"));
#endif

	// Network output
	if (net_output_period < 0)	
		net_output_period = n+1; // disable periodic network output
	int stim_net_output_period = n+1; //tenth_sec; // number of timesteps per which network output is generated during stimulation
	int ten_secs_before_recall = int(round(28800./dt));
	int ten_secs_after_recall = int(round(28820./dt));
	int ten_mins_after_recall = int(round(29410./dt));
	int one_hour_after_recall = int(round(32410./dt));
	net_output = vector<int> {tb_stim_start, tb_stim_start+tenth_sec, tb_stim_learn_end-tenth_sec, tb_recall_start, tb_recall_start+tenth_sec, ten_secs_before_recall,
	                          ten_secs_after_recall, ten_mins_after_recall, one_hour_after_recall}; // selected times at which network output is generated during stimulation
	//net_output = vector<int> {tb_stim_start, tb_stim_start+tenth_sec, tb_stim_learn_end-tenth_sec, tb_recall_start, tb_recall_start+tenth_sec}; // selected times at which network output is generated during stimulation
	for (int ci = 0; ci < net_output.size(); ci++)
	{
		if (net_output[ci] % output_period != 0)
		{
			int net_out_new = net_output[ci] + (net_output[ci] % output_period);

			cout << "WARNING: 'net_output[" << ci << "] = " << net_output[ci] << "' cannot be produced with the current setting of 'output_period'. It will be adjusted to 'net_output[" << ci << "] = " << net_out_new << "'." << endl;
			net_output[ci] = net_out_new;
		}
	}
	int rich_comp_buffer; // buffer for rich computation
	if (stim_net_output_period > n)
		rich_comp_buffer = int(10./dt) + wfr/2.;
	else
		rich_comp_buffer = stim_net_output_period + wfr/2.;

	// ==============================================================================================================================
		// Check if files have been opened properly and if yes, set precision

		if ( !txt_fr.is_open() || !gpl_fr.is_open() || !txt_fr_inh.is_open() || !gpl_fr_inh.is_open() || !gplscript.is_open()
				 || !gpl_cn.is_open() || !txt_cn.is_open() || !txt_data.is_open() || !txt_mean.is_open()
	#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
				 || !txt_spike_number.is_open()
	#endif
	#if SPIKE_PLOTTING == RASTER || SPIKE_PLOTTING == NUMBER_AND_RASTER
				 || !txt_spike_raster.is_open()
	#endif
			 )
		{
			cout << "Unable to open file!" << endl << separator << endl;
			return -1;
		}

		txt_data.precision(OUTPUT_PRECISION);
		txt_mean.precision(OUTPUT_PRECISION);
		txt_spike_number.precision(OUTPUT_PRECISION);
		txt_spike_raster.precision(OUTPUT_PRECISION);


// ==============================================================================================================================
	// Output script for connection plot

	int total_c_count_exc_exc = 0; // total number of exc.->exc. connections
	int total_c_count_inh_exc = 0; // total number of inh.->exc. connections
	int total_c_count_exc_inh = 0; // total number of exc.->inh. connections
	int total_c_count_inh_inh = 0; // total number of inh.->inh. connections

	for (int m=0; m<pow2(Nl)+pow2(Nl_inh); m++) // consecutive neuron numbering
	{
		const int c_count_exc = net.getNumberIncoming(TYPE_EXC, m); // number of incoming excitatory connections to neuron m
		const int c_count_inh = net.getNumberIncoming(TYPE_INH, m); // number of incoming inhibitory connections to neuron m

		if (m < pow2(Nl)) // m is an excitatory neuron
		{
			const int k = row(m); // get row number of neuron m
			const int l = col(m); // get column number of neuron m

			total_c_count_exc_exc += c_count_exc; // E->E
			total_c_count_inh_exc += c_count_inh; // I->E

			txt_cn << fixed << k << "\t\t" << l << "\t\t" << c_count_exc+c_count_inh << endl;

			if (l == Nl)
				txt_cn << endl; // empty line for color plot
		}
		else // m is an inhibitory neuron
		{
			total_c_count_exc_inh += c_count_exc; // E->I
			total_c_count_inh_inh += c_count_inh; // I->I
		}
	}

	txt_cn.close();
	createNetworkColorPlot(gpl_cn, Nl, -1.0, 3, "inc_connection", "", true, "# of incoming connections");
	gplscript << "gnuplot inc_connection_map.gpl" << endl;

	cout.precision(4);
	cout << fixed << "E->E connectivity: " << total_c_count_exc_exc / double(pow2(Nl*Nl)-pow2(Nl)) * 100 << " % (expected: " << pc*100 << " %)" << endl;
	cout << fixed << "I->E connectivity: " << total_c_count_inh_exc / double(pow2(Nl*Nl_inh)) * 100 << " % (expected: " << pc*100 << " %)" << endl;
	cout << fixed << "E->I connectivity: " << total_c_count_exc_inh / double(pow2(Nl*Nl_inh)) * 100 << " % (expected: " << pc*100 << " %)" << endl;
	cout << fixed << "I->I connectivity: " << total_c_count_inh_inh / double(pow2(Nl_inh*Nl_inh)-pow2(Nl_inh)) * 100 << " % (expected: " << pc*100 << " %)" << endl;

// ==============================================================================================================================
	// Actual simulation
	int compmode = 1; // indicates the computation mode
#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
	int total_spike_num = 0; // total number of spikes in currently considered output_period interval
#endif
	// Write headline to data file
	txt_data << "#Time\t\t"; // time
	for (int nn=0; nn < exc_neuron_output.size(); nn++)
	{
		txt_data << fixed
		         << "V(" << exc_neuron_output[nn] << ")    \t\t" // exc. neuron voltage
#if NEURON_MODEL == LIF
		         << "I_tot(" << exc_neuron_output[nn] << ")\t\t" // exc. neuron total input current
#elif NEURON_MODEL == MAT2
		         << "theta(" << exc_neuron_output[nn] << ")\t\t" // exc. neuron threshold
#endif
#if PROTEIN_POOLS == POOLS_C
		         << "p^C(" << exc_neuron_output[nn] << ")\t\t"; // exc. neuron unspecific protein amount
#elif PROTEIN_POOLS == POOLS_PD
		         << "p^LTP(" << exc_neuron_output[nn] << ")\t\t" // exc. neuron LTP protein amount
		         << "p^LTD(" << exc_neuron_output[nn] << ")\t\t"; // exc. neuron LTD protein amount
#elif PROTEIN_POOLS == POOLS_PCD
		         << "p^LTP(" << exc_neuron_output[nn] << ")\t\t" // exc. neuron LTP protein amount
			 << "p^C(" << exc_neuron_output[nn] << ")\t\t" // exc. neuron unspecific protein amount
		         << "p^LTD(" << exc_neuron_output[nn] << ")\t\t"; // exc. neuron LTD protein amount
#endif
	}

	for (int nn=0; nn < inh_neuron_output.size(); nn++)
		txt_data << fixed
		         << "V(" << inh_neuron_output[nn] << ")    \t\t" // inh. neuron voltage
#if NEURON_MODEL == LIF
		         << "I_tot(" << inh_neuron_output[nn] << ")\t\t"; // inh. neuron total input current
#elif NEURON_MODEL == MAT2
		         << "theta(" << inh_neuron_output[nn] << ")\t\t"; // inh. neuron threshold
#endif

	for (int sn=0; sn < synapse_output.size(); sn++)
	{
		txt_data << fixed
		         << "h(" << synapse_output[sn].presyn_neuron << "," << synapse_output[sn].postsyn_neuron << ")\t\t" // early-phase synaptic strength
		         << "z(" << synapse_output[sn].presyn_neuron << "," << synapse_output[sn].postsyn_neuron << ")\t\t" // late-phase synaptic strength
		         << "Ca(" << synapse_output[sn].presyn_neuron << "," << synapse_output[sn].postsyn_neuron << ")"; // synaptic calcium amount
		if (sn < synapse_output.size()-1)
			txt_data << "\t\t";
	}
	txt_data << "\r\n";

  	// Write headline to mean weight file
	txt_mean << "#Time\t\tEarly mean CA\t\tEarly SD CA\t\tLate mean CA\t\tLate SD CA\t\tProtein mean CA\t\tProtein SD CA\t\t"
	         << "Early mean control\t\tEarly SD control\t\tLate mean control\t\tLate SD control\t\tProtein mean control\t\tProtein SD control\r\n";

	// Time loop
	for (int j = tb_start; j <= n; j++)
	{
		int spike_num; // number of spikes that have occurred in this timestep in the whole network

		// update percentage
		double p_new = double(round(double(j-tb_start) / double(n-tb_start) * 1000.0)) / 10.0; // round to first decimal place
		if (p_new > p)
		{
			p = p_new;
			printf("\rProgress: %.1f %% completed.", p);
			fflush(stdout);
		}

		// in rich state (full, detailed computation)
		if (compmode == 1)
		{
			// Calculate next step for Network
#if SPIKE_PLOTTING == RASTER || SPIKE_PLOTTING == NUMBER_AND_RASTER
			spike_num = net.processTimeStep(j, &txt_spike_raster);
#else
			spike_num = net.processTimeStep(j);
#endif

			if (j == tb_start && st_full.isSet()) // in the very first timestep
			{
				// Set stimulation
#ifdef TWO_NEURONS_ONE_SYNAPSE
				// only stimulate one neuron
				net.setSingleNeuronStimulus(1, st_learn);

	#ifdef PLASTICITY_OVER_FREQ
				// stimulate another neuron
				net.setSingleNeuronStimulus(0, st_learn2);
	#endif

#elif defined INTERMEDIATE_RECALL
				if (tb_start == 0) // real learning stimulus
					net.setBlockStimulus(st_learn, CORE_SIZE); // learning stimulus for the first CORE_SIZE neurons
				else // "fake learning stimulus" used for intermediate recall
					net.setRandomStimulus(st_learn, int(round(recall_fraction*CORE_SIZE)), &logf, 0, CORE_SIZE); // intermediate recall stimulus for randomly selected recall_fraction*CORE_SIZE neurons
				// "real recall stimulus" is assigned right before it starts, to prevent interference
#else
				net.setBlockStimulus(st_learn, CORE_SIZE); // learning stimulus for the first CORE_SIZE neurons
#endif
				// Plot stimulation
				//st_learn.plotAll("learning_stimulation");
				//st_full.plotAll("full_stimulation"); // may use a lot of disk space if recall and learning are far apart

			}

			if (j == tb_stim_start-1) // one step before learning stimulus
			{
#ifndef INTERMEDIATE_RECALL
				net.resetPlasticity(true, true, true, true); // set all plastic changes to zero (reset changes from control stimulus for instance)
				//net.reset(); // reset network (only for the case that a control recall was used before)
#endif
#if STIPULATE_CA == ON
				net.stipulateFirstNeuronsAssembly(CORE_SIZE); // stipulate assembly weights
#endif
			}

			if (j == tb_recall_start-1) // one step before recall stimulus
			{
#if SAVE_NET_STATE == ON
				net.saveNetworkState("saved_state.txt", j);
				if (!copyFile("saved_state.txt", "../saved_state0.txt"))
					throw runtime_error(string("Network state file could not be copied to upper directory."));
				else
					remove("saved_state.txt");
#endif
				net.setBlockStimulus(st_recall, int(round(recall_fraction*CORE_SIZE))); // recall stimulation for the first recall_fraction*CORE_SIZE neurons
			}

#if FF_AFTER_LEARN == ON || FF_AFTER_STIM == ON || FF_AFTER_NETLOAD == ON
			if ( 
#if FF_AFTER_LEARN == ON
			     (st_learn.isSet() && j == tb_stim_learn_end + rich_comp_buffer)  // end of learning stimulation plus buffer reached
#endif
#if FF_AFTER_STIM == ON
			     || (st_full.isSet() && j == tb_stim_end + rich_comp_buffer)  // end of all stimulation plus buffer reached
#endif
#if FF_AFTER_NETLOAD == ON
			     || (j == tb_start && tb_start > 0) // network state loaded (or TODO: always in the beginning, i.e., also for tb_start = 0?)
#endif
			   )
			{
				vector<double> ps_end = net.getProteinSynthesisEnd();

				if (ff_enabled && j != tb_recall_start)
				{
					logf << "FF mode entered at t = " << dtos(j*dt, 3) << " s" << endl;
					compmode = 2; // enter fast-forward mode
				}

				logf << "Tags will have vanished at t = " << dtos(net.getTagVanishTime(), 3) << " s" << endl;

#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
		         	logf << "Protein synthesis (C) will have ended at t = " << dtos(ps_end[1], 3) << " s" << endl;
#endif
#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
				logf << "Protein synthesis (P) will have ended at t = " << dtos(ps_end[0], 3) << " s" << endl
				     << "Protein synthesis (D) will have ended at t = " << dtos(ps_end[2], 3) << " s" << endl;
#endif

				net.resetPlasticity(false, false, true, false); // sets all calcium values to zero

			}
#endif

		}

		// in fast-forward mode (only computation of late-phase dynamics)
		else if (compmode == 2)
		{
			double delta_t = 50.; // if completely numerical [OBSOLETE]: use timesteps of not more than 1 min - that is small enough to not visibly cut off the peak of the protein curve
			int new_j = int(round(floor(j*dt + delta_t)/dt)); // use floor function to ensure that FF steps end with full seconds

			if (new_j > tb_stim_start-1 && j < tb_stim_start-1) // step in which learning stimulus begins
			{
				new_j = tb_stim_start-2;
				delta_t = (new_j - j) * dt;

				compmode = 1; // go back to rich computation mode
			}
			else if (new_j > tb_recall_start-1 && j < tb_recall_start-1) // step in which recall stimulus begins
			{
				new_j = tb_recall_start-2;
				delta_t = (new_j - j) * dt;

				compmode = 1; // go back to rich computation mode
			}
			else if (new_j > n) // last step
			{
				new_j = n;
				delta_t = (new_j - j) * dt;
			}

			logf << "FF comp. until t = " << new_j*dt << " s" << endl;

			if (!net.processTimeStep_FF(j, delta_t, &logf))
			{
				if (STC)
				{
					logf << "No more late-phase dynamics" << endl;
					STC = false;
				}
			}
			else
			{
				STC = true;
			}	

			if (compmode == 1)
				logf << "Rich mode entered at t = " << new_j*dt << " s" << endl;

			j = new_j;
		}

#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
		total_spike_num += spike_num;
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// OUTPUT FOR PLOTS ///////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////

		if (j % output_period == 0 || compmode == 2) // use steps in intervals of output_period
   		{

#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
			// Write spike data to file
			txt_spike_number << j*dt << "\t\t" << total_spike_num << endl; // time and number of spikes within interval output_period*dt
			total_spike_num = 0; // reset total spike count
#endif

			// Output of selected neuronal and synaptic data to file
			txt_data << j*dt << "\t\t"; // time

			for (int nn=0; nn < exc_neuron_output.size(); nn++)
			{
				txt_data << fixed
				         << net.getVoltage(exc_neuron_output[nn]) << "\t\t" // exc. neuron voltage
#if NEURON_MODEL == LIF
				         << net.getCurrent(exc_neuron_output[nn]) << "\t\t" // exc. neuron total input current
#elif NEURON_MODEL == MAT2
				         << net.getThreshold(exc_neuron_output[nn]) << "\t\t" // exc. neuron threshold
#endif
#if PROTEIN_POOLS == POOLS_C
				         << net.getCProteinAmount(exc_neuron_output[nn]) << "\t\t"; // exc. neuron unspecific protein amount
#elif PROTEIN_POOLS == POOLS_PD

				         << net.getPProteinAmount(exc_neuron_output[nn]) << "\t\t" // exc. neuron LTP protein amount
				         << net.getDProteinAmount(exc_neuron_output[nn]) << "\t\t"; // exc. neuron LTD protein amount
#elif PROTEIN_POOLS == POOLS_PCD
				         << net.getPProteinAmount(exc_neuron_output[nn]) << "\t\t" // exc. neuron LTP protein amount
				         << net.getCProteinAmount(exc_neuron_output[nn]) << "\t\t" // exc. neuron unspecific protein amount
				         << net.getDProteinAmount(exc_neuron_output[nn]) << "\t\t"; // exc. neuron LTD protein amount
#endif
			}

			for (int nn=0; nn < inh_neuron_output.size(); nn++)
				txt_data << fixed
				         << net.getVoltage(inh_neuron_output[nn]) << "\t\t" // inh. neuron voltage
#if NEURON_MODEL == LIF
				         << net.getCurrent(inh_neuron_output[nn]) << "\t\t"; // inh. neuron total input current
#elif NEURON_MODEL == MAT2
				         << net.getThreshold(inh_neuron_output[nn]) << "\t\t"; // inh. neuron threshold
#endif

			for (int sn=0; sn < synapse_output.size(); sn++)
			{
				txt_data << fixed
				         << net.getEarlySynapticStrength(synapse_output[sn]) << "\t\t" // early-phase synaptic strength
				         << net.getLateSynapticStrength(synapse_output[sn]) << "\t\t" // late-phase synaptic strength
				         << net.getSynapticCalcium(synapse_output[sn]); // synaptic calcium amount
				if (sn < synapse_output.size()-1)
					txt_data << "\t\t";
			}
			txt_data << "\r\n";

#ifndef TWO_NEURONS_ONE_SYNAPSE
			// Output of mean and std. dev. of the weights within subpopulations
			txt_mean << j*dt << "\t\t"; // time

			// get mean and std. dev. of the early- and late weight and of the protein amount within the cell assembly
			double early_mean = net.getMeanEarlySynapticStrength(CORE_SIZE);
			double late_mean = net.getMeanLateSynapticStrength(CORE_SIZE);
			double prot_mean = net.getMeanCProteinAmount(CORE_SIZE);
			int size_control = pow2(Nl) - CORE_SIZE;

			txt_mean << early_mean << "\t\t" << net.getSDEarlySynapticStrength(early_mean, CORE_SIZE) << "\t\t"
			         << late_mean << "\t\t" << net.getSDLateSynapticStrength(late_mean, CORE_SIZE) << "\t\t"
 			         << prot_mean << "\t\t" << net.getSDCProteinAmount(late_mean, CORE_SIZE) << "\t\t";

			// get mean and std. dev. of the early- and late weight and of the protein amount within the control subpopulation
 			early_mean = net.getMeanEarlySynapticStrength(size_control, CORE_SIZE);
 			late_mean = net.getMeanLateSynapticStrength(size_control, CORE_SIZE);
			prot_mean = net.getMeanCProteinAmount(size_control, CORE_SIZE);

 			txt_mean << early_mean << "\t\t" << net.getSDEarlySynapticStrength(early_mean, size_control, CORE_SIZE) << "\t\t"
 			         << late_mean << "\t\t" << net.getSDLateSynapticStrength(late_mean, size_control, CORE_SIZE) << "\t\t"
 			         << prot_mean << "\t\t" << net.getSDCProteinAmount(late_mean, size_control, CORE_SIZE);

			txt_mean << "\r\n";


			// Output of network plots (usually every net_output_period; under stimulation every stim_net_output_period) and clearing of spike time list
			if ( (j % net_output_period == 0)
			   || ((j % stim_net_output_period == 0) && st_full.stimExists(j-stim_net_output_period)) // creates network plots every stim_net_output_period-th timestep
			                                                                                          // until stim_net_output_period after the end of stimulation
			   || (netOutput(j)) // creates network plots at specified times
			   /*|| (compmode == 2) */) // if in fast-forward computation mode [OBSOLETE: net_output_period should now be chosen to be a multiple of the FF mode timestep]

			{
				if (j < (n - wfr/2.)) // if there are still enough timesteps left to compute the firing rate
				{
					double t = j*dt;
					txt_net_t = new ofstream(string("network_plots/") + dateStr("_net_") + dtos(t,1) + string(".txt"));
					txt_net_t->precision(NET_OUTPUT_PRECISION);

					if (!txt_net_t->is_open())
						cout << "Network data file could not be opened!" << endl;


					// Output of early-phase matrix at time 't'
					for (int m=0; m < pow2(Nl); m++)
					{
						for (int n=0; n < pow2(Nl); n++)
						{
							*txt_net_t << fixed << net.getEarlySynapticStrength(synapse(m,n));

							if (n < pow2(Nl) - 1)
								*txt_net_t << "\t\t";
						}
						*txt_net_t << endl; // next neuron row begins
					}
					*txt_net_t << endl;

					// Output of late-phase matrix at time 't'
					for (int m=0; m < pow2(Nl); m++)
					{
						for (int n=0; n < pow2(Nl); n++)
						{
							*txt_net_t << fixed << net.getLateSynapticStrength(synapse(m,n));

							if (n < pow2(Nl) - 1)
								*txt_net_t << "\t\t";
						}
						*txt_net_t << endl; // next neuron row begins
					}
					*txt_net_t << endl;

					// Output of firing rate matrix at (past) time 'tprime' and creating plot
					if (data_to_plot)
					{
						instFiringRates(txt_net_tprime, jprime);
						createNetworkPlotAveragedWeights(jprime*dt, h_0, Nl, z_max);
						if (jprime == 0)
							createNetworkPlotWeights(jprime*dt, h_0, Nl, z_max);
					}

					txt_net_tprime = txt_net_t; // save file handle for next plot step
					jprime = j; // save time step for next plot step
					data_to_plot = true;
				}
			}
#endif
		}

	} // end of for(j)

	cout << endl;

#ifndef TWO_NEURONS_ONE_SYNAPSE
	// Finish remaining network output by printing firing rate matrix for (past) time 'tprime' and creating plot
	if (data_to_plot)
	{
		instFiringRates(txt_net_tprime, jprime);
		createNetworkPlotAveragedWeights(jprime*dt, h_0, Nl, z_max);
	}
#endif

// ==============================================================================================================================

	// Compute mean firing rate and standard deviation of the firing rate for both exc. and inh. population
	double mfr = 0.; // mean firing rate of the whole excitatory population
	double mfr_inh = 0.; // mean firing rate of the whole inhibitory population
	double sdfr = 0.; // standard deviation of the firing rate of the whole excitatory population
	double sdfr_inh = 0.; // standard deviation of the firing rate of the whole inhibitory population
	int sdfr_part1 = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
	int sdfr_part2 = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
	int sdfr_part1_inh = 0; // part 1 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts
	int sdfr_part2_inh = 0; // part 2 of the standard deviation (due to Steiner's translation theorem) - integer because it contains only spike counts

	for (int m=0; m<pow2(Nl)+pow2(Nl_inh); m++)
	{
		int nu = net.getSpikeCount(m); // get number of spikes of neuron m

		if (m < pow2(Nl)) // neuron m is in excitatory population
		{
			mfr += double(nu);
			sdfr_part1 += pow2(nu);
			sdfr_part2 += nu;
		}
		else // neuron m is in inhibitory population
		{
			mfr_inh += double(nu);
			sdfr_part1_inh += pow2(nu);
			sdfr_part2_inh += nu;
		}
	}

	mfr /= (pow2(Nl) * t_max); // compute the mean firing rate of the exc. population from the number of spikes
	if (Nl_inh > 0)
		mfr_inh /= (pow2(Nl_inh) * t_max); // compute the mean firing rate of the inh. population from the number of spikes
	sdfr = sqrt( ( double(sdfr_part1) - double(pow2(sdfr_part2)) / pow2(Nl) ) /
	             double(pow2(Nl)-1) ) / t_max; // compute standard deviation of the firing rate according to Steiner's translation 																													// deviation of the firing rates according to Steiner's translation theorem
	if (Nl_inh > 0)
		sdfr_inh = sqrt( ( double(sdfr_part1_inh) - double(pow2(sdfr_part2_inh)) / pow2(Nl_inh) ) /
		                 double(pow2(Nl_inh)-1) ) / t_max; // compute standard deviation of the firing rate according to Steiner's translation 																										 // theorem
// ==============================================================================================================================
	// Write mean firing rate and standard deviation of the firing rate in log file

	logf << "Mean firing rate (exc. population): " << dtos(mfr, 6, true) << " +- " << dtos(sdfr, 6, true) << endl
	     << "Mean firing rate (inh. population): " << dtos(mfr_inh, 6, true) << " +- " << dtos(sdfr_inh, 6, true) << endl;


#ifdef SEEK_I_0
	*seekic = mfr; // return mean firing rate in "SEEK_I_0" mode
#endif


// ==============================================================================================================================
	// Create mean firing rate map plots

	// Write data files and determine minimum and maximum firing rates
	for (int m=0; m < pow2(Nl)+pow2(Nl_inh); m++)
	{
		double fr = double(net.getSpikeCount(m)) / t_max;

		if (m < pow2(Nl)) // neuron m is in excitatory population
		{
			if (fr > max_firing_rate)
				max_firing_rate = fr;
			if (fr < min_firing_rate)
				min_firing_rate = fr;

			txt_fr << fixed << col(m) << "\t\t" << row(m) << "\t\t" << fr << endl;

			if ((m+1) % Nl == 0)
				txt_fr << endl; // another free line
		}
		else // neuron m is in inhibitory population
		{
			if (fr > max_firing_rate_inh)
				max_firing_rate_inh = fr;
			if (fr < min_firing_rate_inh)
				min_firing_rate_inh = fr;

			int m_eff = m - pow2(Nl);
			txt_fr_inh << fixed << colG(m_eff, Nl_inh) << "\t\t" << rowG(m_eff, Nl_inh) << "\t\t" << fr << endl;

			if ((m_eff+1) % Nl_inh == 0)
				txt_fr_inh << endl; // another free line
		}
	}

	// Close data files
	txt_fr.close();
	txt_fr_inh.close();

	// Create firing rate plot of excitatory population
	createNetworkColorPlot(gpl_fr, Nl, -1, 3, "fr_exc", "", true, "{/Symbol n} / Hz", min_firing_rate, max_firing_rate);
	gplscript << "gnuplot fr_exc_map.gpl" << endl;

	// Create firing rate plot of inhibitory population
	createNetworkColorPlot(gpl_fr_inh, Nl_inh, -1, 3, "fr_inh", "", true, "{/Symbol n} / Hz", min_firing_rate_inh, max_firing_rate_inh);
	gplscript << "gnuplot fr_inh_map.gpl" << endl;

// ==============================================================================================================================
	// Create plot of spike number per output_period
#if SPIKE_PLOTTING == NUMBER || SPIKE_PLOTTING == NUMBER_AND_RASTER
	txt_spike_number.close();

	createSpikeNumberPlot(gplscript, t_max, output_period*dt);
#endif
	// Create spike raster plot
#if SPIKE_PLOTTING == RASTER || SPIKE_PLOTTING == NUMBER_AND_RASTER
	txt_spike_raster.close();

	//if (tb_stim_end > tb_stim_start)
	//	createSpikeRasterPlot(gplscript, 0.9*tb_stim_start*dt, 1.1*tb_stim_end*dt, pow2(Nl), pow2(Nl_inh));
	//else
	//	createSpikeRasterPlot(gplscript, 0., t_max, pow2(Nl), pow2(Nl_inh), recall_fraction*CORE_SIZE, CORE_SIZE);

	createSpikeRasterPlot(gplscript, (t_max > 28809.99 ? 28809.99 : 0.), (t_max > 28810.20 ? 28810.20 : t_max), pow2(Nl), pow2(Nl_inh), recall_fraction*CORE_SIZE, CORE_SIZE);
#endif

// ==============================================================================================================================
	// Create neuron and synapse plots
	txt_data.close();

	createExcNeuronPlot(exc_neuron_output, gplscript, t_max);
#if PROTEIN_POOLS == POOLS_C  // spare one columns for exc. neurons because of protein synthesis
	createInhNeuronPlot(inh_neuron_output, gplscript, t_max, 3*exc_neuron_output.size());
	createSynapsePlot(synapse_output, gplscript, t_max, prot_learn, h_0, net.getThreshold(1,2), net.getThreshold(1,3), net.getThreshold(1,1), net.getThreshold(2,1),
                          3*exc_neuron_output.size()+2*inh_neuron_output.size());
#elif PROTEIN_POOLS == POOLS_PD  // spare two columns for exc. neurons because of protein synthesis
	createInhNeuronPlot(inh_neuron_output, gplscript, t_max, 4*exc_neuron_output.size());
	createSynapsePlot(synapse_output, gplscript, t_max, prot_learn, h_0, net.getThreshold(1,2), net.getThreshold(1,3), net.getThreshold(1,1), net.getThreshold(2,1),
                          4*exc_neuron_output.size()+2*inh_neuron_output.size());
#elif PROTEIN_POOLS == POOLS_PCD  // spare three columns for exc. neurons because of protein synthesis
	createInhNeuronPlot(inh_neuron_output, gplscript, t_max, 5*exc_neuron_output.size());
	createSynapsePlot(synapse_output, gplscript, t_max, prot_learn, h_0, net.getThreshold(1,2), net.getThreshold(1,3), net.getThreshold(1,1), net.getThreshold(2,1),
                          5*exc_neuron_output.size()+2*inh_neuron_output.size());
#endif

// ==============================================================================================================================
	txt_mean.close();
#ifndef TWO_NEURONS_ONE_SYNAPSE
	// Create plots of the mean and std. dev. of the weight and of the protein amount

	createMeanWeightPlotCA(gplscript, t_max, h_0); // in the cell assembly
	createMeanWeightPlotControl(gplscript, t_max, h_0); // in the control subpopulation
#endif
// ==============================================================================================================================
	// Display final information and save parameters
	int tsec = timeMeasure(false);
	string el_time = "Elapsed time: ";

	if (tsec < 60)
        	el_time += dtos(tsec) + " s";
    	else if (tsec < 3600)
        	el_time += dtos(floor(tsec / 60.)) + " m " + dtos(tsec % 60) + " s";
    	else
	{
		if (tsec > 48*3600) // show number of full days if more than 48 hrs. have passed
		{
			int days = int(floor(tsec / double(24*3600)));
			el_time += dtos(days) + " d ";
			tsec -= 24*3600*days;
		}

        	el_time += dtos(floor(tsec / 3600.)) + " h " + dtos(floor((tsec % 3600) / 60.)) + " m " + dtos((tsec % 3600) % 60) + " s";
	}

	cout << el_time << endl;
	addToParamsFile(el_time); // add additional information to parameter file: elapsed time

#ifdef SEEK_I_0
	if (chdir("..") == -1)
		showChDirErrMessage();
#endif

	cout << separator << endl;

// ==============================================================================================================================
	// Free allocated memory and close files that have remained open
	logf.close();
	gplscript.close(); // script is not executed, it is only created in case it is needed later on
	return 0;
}

#ifdef SEEK_I_0
/*** setSeekICVar ***
 * Sets the variable to communicate with the main(...) function for seeking I_0 */
void setSeekICVar(double *_seekic)
{
	seekic = _seekic;
}
#endif


/*** setParams ***
 * Sets the simulation parameters on given values and resets network(s) */
void setParams(double _I_0, double _sigma_WN, double _tau_syn, double _w_ee, double _w_ei, double _w_ie, double _w_ii, 
               string _prot_learn, string _prot_recall, int _output_period, int _net_output_period, int _N_stim, double _theta_p, double _theta_d,
               double _Ca_pre, double _Ca_post, double _theta_pro_P, double _theta_pro_C, double _theta_pro_D, double _recall_fraction)
{
	tau_syn = _tau_syn;
	prot_learn = _prot_learn;
	prot_recall = _prot_recall;
	recall_fraction = _recall_fraction;
	output_period = _output_period;
	net_output_period = _net_output_period;
	N_stim = _N_stim;
	net.setConstCurrent(_I_0);
	net.setSigma(_sigma_WN);
	net.setSynTimeConstant(_tau_syn);
	w_ei = _w_ei;
	w_ie = _w_ie;
	w_ii = _w_ii;
	net.setCouplingStrengths(_w_ee, _w_ei, _w_ie, _w_ii);
	net.setCaConstants(_theta_p, _theta_d, _Ca_pre, _Ca_post);
	net.setPSThresholds(_theta_pro_P, _theta_pro_C, _theta_pro_D);
	net.setSpikeStorageTime(net_output_period + int(round(wfr/2.)));
	net.reset();
}

/*** Constructor ***
 * Sets all parameters on given values and calls constructors for Neuron instances */
NetworkSimulation(int _Nl, int _Nl_inh, double _dt, double _t_max,
                  double _pc, double _sigma_plasticity, double _z_max, double _t_wfr, bool _ff_enabled)
	: Nl(_Nl), Nl_inh(_Nl_inh), dt(_dt), t_max(_t_max), pc(_pc),
	  t_wfr(_t_wfr), wfr(_t_wfr/dt), ff_enabled(_ff_enabled), z_max(_z_max),
	  net(_dt, _Nl, _Nl_inh, _pc, _sigma_plasticity, _z_max)
{

}

};
