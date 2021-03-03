/**************************************************************************************************
 ***     Model of a network of neurons with long-term plasticity between excitatory neurons     ***
 **************************************************************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

#include <random>
#include <sstream>

using namespace std;

#include "Neuron.cpp"

struct synapse // structure for synapse definition
{
	int presyn_neuron; // the number of the presynaptic neuron
	int postsyn_neuron; // the number of the postsynaptic neuron
	synapse(int _presyn_neuron, int _postsyn_neuron) // constructor
	{
		presyn_neuron = _presyn_neuron;
		postsyn_neuron = _postsyn_neuron;
	}
};

/*** Network class ***
 * Represents a network of neurons */
class Network {

#if (PROTEIN_POOLS != POOLS_C && PROTEIN_POOLS != POOLS_PD && PROTEIN_POOLS != POOLS_PCD)
	#error "Unsupported option for PROTEIN_POOLS."
#endif

friend class boost::serialization::access;

private:

/*** Computational parameters ***/
double dt; // s, one time step for numerical simulation
int N; // total number of excitatory plus inhibitory neurons
int t_syn_delay_steps; // constant t_syn_delay converted to time steps
int t_Ca_delay_steps; // constant t_Ca_delay converted to time steps

/*** State variables ***/
vector<Neuron> neurons; // vector of all N neuron instances (first excitatory, then inhibitory)
bool** conn; // the binary connectivity matrix, the main diagonal is zero (because there is no self-coupling)
double** Ca; // the matrix of postsynaptic calcium concentrations
double** h; // the matrix of early-phase coupling strengths
double** z; // the matrix of late-phase coupling strengths
int* last_Ca_spike_index; // contains the indices of the last spikes that were important for calcium dynamics
minstd_rand0 rg; // default uniform generator for random numbers to establish connections (seed is chosen in constructor)
uniform_real_distribution<double> u_dist; // uniform distribution, constructed in Network class constructor
normal_distribution<double> norm_dist; // normal distribution to obtain Gaussian white noise, constructed in Network class constructor
int stimulation_end; // timestep by which all stimuli have ended
double* sum_h_diff; // sum of all early-phase changes for each postsynaptic neuron
double* sum_h_diff_p; // sum of E-LTP changes for each postsynaptic neuron
double* sum_h_diff_d; // sum of E-LTD changes for each postsynaptic neuron

protected:

/*** Physical parameters ***/
int Nl; // number of neurons in one line (row or column) of the exc. population (better choose an odd number, for there exists a "central" neuron)
int Nl_inh; // number of neurons in one line (row or column) of the inh. population (better choose an odd number, for there exists a "central" neuron)
double tau_syn; // s, the synaptic time constant
double t_syn_delay; // s, the synaptic transmission delay for PSPs - has to be at least one time step!
double p_c; // connection probability (prob. that a directed connection exists)
double w_ee; // nC, magnitude of excitatory PSP effecting an excitatory postsynaptic neuron
double w_ei; // nC, magnitude of excitatory PSP effecting an inhibitory postsynaptic neuron
double w_ie; // nC, magnitude of inhibitory PSP effecting an excitatory postsynaptic neuron
double w_ii; // nC, magnitude of inhibitory PSP effecting an inhibitory postsynaptic neuron

/*** Plasticity parameters ***/
double t_Ca_delay; // s, delay for spikes to affect calcium dynamics - has to be at least one time step!
double Ca_pre; // s^-1, increase in calcium current evoked by presynaptic spike
double Ca_post; // s^-1, increase in calcium current evoked by postsynaptic spike
double tau_Ca; // s, time constant for calcium dynamics
double tau_Ca_steps; // time constant for calcium dynamics in timesteps
double tau_h; // s, time constant for early-phase plasticity
double tau_pp; // h, time constant of LTP-related protein synthesis
double tau_pc; // h, time constant of common protein synthesis
double tau_pd; // h, time constant of LTD-related protein synthesis
double tau_z; // min, time constant of consolidation
double gamma_p; // constant for potentiation process
double gamma_d; // constant for depression process
double theta_p; // threshold for calcium concentration to induce potentiation
double theta_d; // threshold for calcium concentration to induce depotentiation
double sigma_plasticity; // nA s, standard deviation of plasticity noise
double alpha_p; // LTP-related protein synthesis rate
double alpha_c; // common protein synthesis rate
double alpha_d; // LTD-related protein synthesis rate
double h_0; // nA, initial value for early-phase plasticity
double theta_pro_p; // nA s, threshold for LTP-related protein synthesis
double theta_pro_c; // nA s, threshold for common protein synthesis
double theta_pro_d; // nA s, threshold for LTD-related protein synthesis
double theta_tag_p; // nA s, threshold for LTP-related tag
double theta_tag_d; // nA s, threshold for LTD-related tag
double z_max; // upper z bound

public:
#ifdef TWO_NEURONS_ONE_SYNAPSE
bool tag_glob; // specifies if a synapse was tagged ever
bool ps_glob; // specifies if protein synthesis ever occurred in any neuron
#endif

double max_dev; // maximum deviation from h_0 (deviation of the synapse with the largest change)
int tb_max_dev; // time bin at which max_dev was encountered
#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
double max_sum_diff; // maximum sum of early-phase changes (sum of the neuron with the most changes)
int tb_max_sum_diff; // time bin at which max_sum_diff was encountered
#endif
#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
double max_sum_diff_p; // maximum sum of LTP early-phase changes (sum of the neuron with the most changes)
int tb_max_sum_diff_p; // time bin at which max_sum_diff_p was encountered
double max_sum_diff_d; // maximum sum of LTD early-phase changes (sum of the neuron with the most changes)
int tb_max_sum_diff_d; // time bin at which max_sum_diff_d was encountered
#endif

/*** rowG (macro) ***
 * Returns the row number for element n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number (general case for a row/column size of d) *
 * - int n: the consecutive element number *
 * - int d: the row/column size */
#define rowG(n, d) ((((n) - ((n) % d)) / d) + 1)

/*** colG (macro) ***
 * Returns the column number for element n (in consecutive numbering), be  *
 * aware that it starts with one, unlike the consecutive number (general case for a row/column size of d) *
 * - int n: the consecutive element number *
 * - int d: the row/column size */
#define colG(n, d) (((n) % d) + 1)

/*** cNN (macro) ***
 * Returns a consecutive number for excitatory neuron (i|j) rather than a pair of numbers like (i|j), be  *
 * aware that it starts with zero, unlike i and j *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
#define cNN(i, j) (((i)-1)*Nl + ((j)-1))

/*** row (macro) ***
 * Returns the row number for excitatory neuron n, be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
#define row(n) (rowG(n, Nl))

/*** col (macro) ***
 * Returns the column number for excitatory neuron n, be  *
 * aware that it starts with one, unlike the consecutive number *
 * - int n: the consecutive neuron number */
#define col(n) (colG(n, Nl))

/*** symm (macro) ***
 * Returns the number of the symmetric element for an element given  *
 * by its consecutive number *
 * - int n: the consecutive element number */
#define symm(n) (cNN(col(n),row(n)))


/*** shallBeConnected ***
 * Draws a uniformly distributed random number from the interval 0.0 to 1.0 and returns, *
 * depending on the connection probability, whether or not a connection shall be established *
 * - int m: consecutive number of presynaptic neuron *
 * - int n: consecutive number of postsynaptic neuron *
 * - return: true if connection shall be established, false if not */
bool shallBeConnected(int m, int n)
{
#ifdef TWO_NEURONS_ONE_SYNAPSE
	// in this paradigm, there is only one synapse from neuron 1 to neuron 0
	if (m == 1 && n == 0)
	{
		neurons[m].addOutgoingConnection(n, TYPE_EXC);
		return true;
	}
#else
	// exc.->exc. synapse
	if (m < pow2(Nl) && n < pow2(Nl))
	{
		if (u_dist(rg) <= p_c) // draw random number
		{
			neurons[n].incNumberIncoming(TYPE_EXC);
			neurons[m].addOutgoingConnection(n, TYPE_EXC);
			return true;
		}
	}

	// exc.->inh. synapse
	else if (m < pow2(Nl) && n >= pow2(Nl))
	{
		if (u_dist(rg) <= p_c) // draw random number
		{
			neurons[n].incNumberIncoming(TYPE_EXC);
			neurons[m].addOutgoingConnection(n, TYPE_INH);
			return true;
		}

	}

	// inh.->exc. synapse
	else if (m >= pow2(Nl) && n < pow2(Nl))
	{
		if (u_dist(rg) <= p_c) // draw random number
		{
			neurons[n].incNumberIncoming(TYPE_INH);
			neurons[m].addOutgoingConnection(n, TYPE_EXC);
			return true;
		}
	}

	// inh.->inh. synapse
	else if (m >= pow2(Nl) && n >= pow2(Nl))
	{
		if (u_dist(rg) <= p_c) // draw random number
		{
			neurons[n].incNumberIncoming(TYPE_INH);
			neurons[m].addOutgoingConnection(n, TYPE_INH);
			return true;
		}
	}
#endif
	return false;
}

/*** areConnected ***
 * Returns whether or not there is a synapse from neuron m to neuron n *
 * - int m: the number of the first neuron in consecutive order *
 * - int n: the number of the second neuron in consecutive order *
 * - return: true if connection from m to n exists, false if not */
bool areConnected(int m, int n) const
{
	if (conn[m][n])
		return true;
	else
		return false;
}

/*** saveNetworkParams ***
 * Saves all the network parameters (including the neuron and channel parameters) to a given file */
void saveNetworkParams(ofstream *f) const
{
	*f << endl;
	*f << "Network parameters:" << endl;
	*f << "N_exc = " << pow2(Nl) << " (" << Nl << " x " << Nl << ")" << endl;
	*f << "N_inh = " << pow2(Nl_inh) << " (" << Nl_inh << " x " << Nl_inh << ")" << endl;
	*f << "tau_syn = "
#if SYNAPSE_MODEL == DELTA
	   << 0
#elif SYNAPSE_MODEL == MONOEXP
	   << tau_syn
#endif
	   << " s" << endl;
	*f << "t_syn_delay = " << t_syn_delay << " s" << endl;
	*f << "h_0 = " << h_0 << " nA s" << endl;
	*f << "w_ee = " << dtos(w_ee/h_0,1) << " h_0" << endl;
	*f << "w_ei = " << dtos(w_ei/h_0,1) << " h_0" << endl;
	*f << "w_ie = " << dtos(w_ie/h_0,1) << " h_0" << endl;
	*f << "w_ii = " << dtos(w_ii/h_0,1) << " h_0" << endl;
	*f << "p_c = " << p_c << endl;
	*f << endl;

	*f << "Plasticity parameters"
#if PLASTICITY == OFF
	   << " <switched off>"
#endif
	   << ": " << endl;
	*f << "t_Ca_delay = " << t_Ca_delay << " s" << endl;
	*f << "Ca_pre = " << Ca_pre << endl;
	*f << "Ca_post = " << Ca_post << endl;
	*f << "tau_Ca = " << tau_Ca << " s" << endl;
	*f << "tau_h = " << tau_h << " s" << endl;
	*f << "tau_pp = " << tau_pp << " h" << endl;
	*f << "tau_pc = " << tau_pc << " h" << endl;
	*f << "tau_pd = " << tau_pd << " h" << endl;
	*f << "tau_z = " << tau_z << " min" << endl;
	*f << "z_max = " << z_max << endl;
	*f << "gamma_p = " << gamma_p << endl;
	*f << "gamma_d = " << gamma_d << endl;
	*f << "theta_p = " << theta_p << endl;
	*f << "theta_d = " << theta_d << endl;
	*f << "sigma_plasticity = " << sigma_plasticity << " nA s" << endl;
	*f << "alpha_p = " << alpha_p << endl;
	*f << "alpha_c = " << alpha_c << endl;
	*f << "alpha_d = " << alpha_d << endl;
	*f << "theta_pro_p = " << theta_pro_p << " nA s" << endl;
	*f << "theta_pro_c = " << theta_pro_c << " nA s" << endl;
	*f << "theta_pro_d = " << theta_pro_d << " nA s" << endl;
	*f << "theta_tag_p = " << theta_tag_p << " nA s" << endl;
	*f << "theta_tag_d = " << theta_tag_d << " nA s" << endl;

	neurons[0].saveNeuronParams(f); // all neurons have the same parameters, take the first one
}

/*** saveNetworkState ***
 * Saves the current state of the whole network to a given file using boost function serialize(...) *
 * - file: the file to read the data from *
 * - tb: current time step */
void saveNetworkState(string file, int tb)
{
	ofstream savefile(file);

	if (!savefile.is_open())
		throw runtime_error(string("Network state could not be saved."));

	boost::archive::text_oarchive oa(savefile);
	oa << tb; // write the current time (in steps) to archive oa
	oa << *this; // write this instance to archive oa

	savefile.close();
}


/*** loadNetworkState ***
 * Load the state of the whole network from a given file using boost function serialize(...); *
 * connectivity matrix 'conn' of the old and the new simulation has to be the same! *
 * - file: the file to read the data from *
 * - return: the simulation time at which the network state was saved (or -1 if nothing was loaded) */
int loadNetworkState(string file)
{
	ifstream loadfile(file);
	int tb;

	if (!loadfile.is_open())
		return -1;

	boost::archive::text_iarchive ia(loadfile);
	ia >> tb; // read the current time (in steps) from archive ia
	ia >> *this; // read this instance from archive ia

	loadfile.close();

	cout << "Network state successfully loaded." << endl;
	return tb;
}

/*** serialize ***
 * Saves all state variables to a file using serialization from boost *
 * - ar: the archive stream *
 * - version: the archive version */
template<class Archive> void serialize(Archive &ar, const unsigned int version)
{
	for (int m=0; m<N; m++)
	{
		ar & neurons[m]; // read/write Neuron instances

		for (int n=0; n<N; n++)
		{
			ar & Ca[m][n]; // read/write matrix of postsynaptic Calcium concentrations
			ar & h[m][n]; // read/write early-phase weight matrix
			ar & z[m][n]; // read/write late-phase weight matrix
		}

		ar & last_Ca_spike_index[m]; // read/write array of the indices of the last spikes that were important for calcium dynamics
		ar & sum_h_diff[m]; // read/write sum of all early-phase changes for each postsynaptic neuron
		ar & sum_h_diff_p[m]; // read/write sum of E-LTP changes for each postsynaptic neuron
		ar & sum_h_diff_d[m]; // read/write sum of E-LTD changes for each postsynaptic neuron
	}
}


/*** processTimeStep ***
 * Processes one timestep (of duration dt) for the network [rich mode / compmode == 1] *
 * - int tb: current time step (for evaluating stimulus and for computing spike contributions) *
 * - ofstream* txt_spike_raster [optional]: file containing spike times for spike raster plot *
 * - return: number of spikes that occurred within the considered timestep in the whole network */
int processTimeStep(int tb, ofstream* txt_spike_raster = NULL)
{
	int spike_count = 0; // number of neurons that have spiked in this timestep
	int st_PSP = tb - t_syn_delay_steps; // presynaptic spike time for evoking PSP in this timestep tb
	int st_CA = tb - t_Ca_delay_steps; // presynaptic spike time for evoking calcium contribution in this timestep tb
	bool STC = false; // specifies if at least one synapse is tagged and receives proteins
	bool ps_neuron = false; // specifies if at least one neuron is exhibiting protein synthesis

	/*******************************************************/
	// compute neuronal dynamics
	for (int m=0; m<N; m++) // loop over neurons (in consecutive order)
	{
		neurons[m].processTimeStep(tb, -1); // computation of individual neuron dynamics

		// add spikes to raster plot and count spikes in this time step
		if (neurons[m].getActivity())
		{
#if SPIKE_PLOTTING == RASTER || SPIKE_PLOTTING == NUMBER_AND_RASTER
			*txt_spike_raster << tb*dt << "\t\t" << m << endl; // add this spike to the raster plot
#endif
			spike_count += 1;
		}

#if COND_BASED_SYN == OFF

	#if SYNAPSE_MODEL == DELTA
		neurons[m].setSynapticCurrent(0.); // reset synaptic current contributions
	#elif SYNAPSE_MODEL == MONOEXP
		neurons[m].setSynapticCurrent(neurons[m].getSynapticCurrent() * exp(-dt/tau_syn)); // exponential decay of previous synaptic current contributions
	#endif

#else

	#if SYNAPSE_MODEL == DELTA
		neurons[m].setExcSynapticCurrent(0.); // reset synaptic current contributions
		neurons[m].setInhSynapticCurrent(0.); // reset synaptic current contributions
	#elif SYNAPSE_MODEL == MONOEXP
		neurons[m].setExcSynapticCurrent(neurons[m].getExcSynapticCurrent() * exp(-dt/tau_syn)); // exponential decay of previous synaptic current contributions
		neurons[m].setInhSynapticCurrent(neurons[m].getInhSynapticCurrent() * exp(-dt/tau_syn)); // exponential decay of previous synaptic current contributions
		//[simple Euler: neurons[m].setExcSynapticCurrent(neurons[m].getExcSynapticCurrent() * (1.-dt/tau_syn)); ]
		//[simple Euler: neurons[m].setInhSynapticCurrent(neurons[m].getInhSynapticCurrent() * (1.-dt/tau_syn)); ]
	#endif

#endif // COND_BASED_SYN == ON

		// Protein dynamics (for neuron m)
#if PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

	#ifdef TWO_NEURONS_ONE_SYNAPSE
		ps_glob = ps_glob || (sum_h_diff[m] >= theta_pro_c);
	#endif

	#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
		double pa_p = neurons[m].getPProteinAmount();
		double pa_d = neurons[m].getDProteinAmount();

		if (sum_h_diff_p[m] > max_sum_diff_p)
		{
			max_sum_diff_p = sum_h_diff_p[m];
			tb_max_sum_diff_p = tb;
		}

		if (sum_h_diff_d[m] > max_sum_diff_d)
		{
			max_sum_diff_d = sum_h_diff_d[m];
			tb_max_sum_diff_d = tb;
		}

		ps_neuron = ps_neuron || (sum_h_diff_p[m] >= theta_pro_p) || (sum_h_diff_d[m] >= theta_pro_d);

		pa_p = pa_p * exp(-dt/(tau_pp * 3600.)) + alpha_p * step(sum_h_diff_p[m] - theta_pro_p) * (1. - exp(-dt/(tau_pp * 3600.)));
		pa_d = pa_d * exp(-dt/(tau_pd * 3600.)) + alpha_d * step(sum_h_diff_d[m] - theta_pro_d) * (1. - exp(-dt/(tau_pd * 3600.)));
		// [simple Euler: pa += (- pa + alpha * step(sum_h_diff - theta_pro)) * (dt / (tau_p * 3600.));]

		sum_h_diff_p[m] = 0.;
		sum_h_diff_d[m] = 0.;
	#else
		double pa_p = 0.;
		double pa_d = 0.;
	#endif

	#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
		double pa_c = neurons[m].getCProteinAmount();

		if (sum_h_diff[m] > max_sum_diff)
		{
			max_sum_diff = sum_h_diff[m];
			tb_max_sum_diff = tb;
		}

		ps_neuron = ps_neuron || (sum_h_diff[m] >= theta_pro_c);

		pa_c = pa_c * exp(-dt/(tau_pc * 3600.)) + alpha_c * step(sum_h_diff[m] - theta_pro_c) * (1. - exp(-dt/(tau_pc * 3600.))); // === ESSENTIAL ===

		sum_h_diff[m] = 0.;
	#else
		double pa_c = 0.;
	#endif

		neurons[m].setProteinAmounts(pa_p, pa_c, pa_d);

#endif // PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC
	}


	/*******************************************************/
	// compute synaptic dynamics
	for (int m=0; m<N; m++) // loop over presynaptic neurons (in consecutive order)
	{
		bool delayed_PSP; // specifies if a presynaptic spike occurred t_syn_delay ago

		// go through presynaptic spikes for PSPs; start from most recent one
		delayed_PSP = neurons[m].spikeAt(st_PSP);
		//delayed_PSP = neurons[m].getActivity(); // in case no synaptic delay is used

#if PLASTICITY == CALCIUM || PLASTICITY == CALCIUM_AND_STC
		bool delayed_Ca = false; // specifies if a presynaptic spike occurred t_Ca_delay ago

		if (m < pow2(Nl)) // plasticity only for exc. -> exc. connections
		{
			// go through presynaptic spikes for calcium contribution; start from last one that was used plus one
			for (int k=last_Ca_spike_index[m]; k<=neurons[m].getSpikeHistorySize(); k++)
			{
				int st = neurons[m].getSpikeTime(k);

				if (st >= st_CA)
				{
					if (st == st_CA) // if presynaptic spike occurred t_Ca_delay ago
					{
						delayed_Ca = true; // presynaptic neuron fired t_Ca_delay ago
						last_Ca_spike_index[m] = k + 1; // next time, start with the next possible spike
					}
					break;
				}
			}
		}
#endif

		/*******************************************************/

		for (int in=0; in<neurons[m].getNumberOutgoing(); in++) // loop over postsynaptic neurons
		{
			int n = neurons[m].getOutgoingConnection(in); // get index (in consecutive order) of postsynaptic neuron
			double h_dev; // the deviation of the early-phase weight from its resting state

			// Synaptic current
			if (delayed_PSP) // if presynaptic spike occurred t_syn_delay ago
			{
				if (neurons[m].getType() == TYPE_EXC)
				{
					double psc; // the postsynaptic current

					if (neurons[n].getType() == TYPE_EXC) // E -> E
					{
						psc = h[m][n] + h_0 * z[m][n];
						neurons[n].increaseExcSynapticCurrent(psc);
					}
					else // E -> I
					{
						psc = w_ei;
						neurons[n].increaseExcSynapticCurrent(psc);
					}
#if DENDR_SPIKES == ON
					neurons[n].updateDendriteInput(psc); // contribution to dendritic spikes
#endif
				}
				else
				{
					if (neurons[n].getType() == TYPE_EXC) // I -> E
					{

						neurons[n].increaseInhSynapticCurrent(w_ie);
					}
					else // I -> I
					{
						neurons[n].increaseInhSynapticCurrent(w_ii);
					}
				}
			}

			// Long-term plasticity
			if (m < pow2(Nl) && n < pow2(Nl)) // plasticity only for exc. -> exc. connections
			{
#if PLASTICITY == CALCIUM || PLASTICITY == CALCIUM_AND_STC
				// Calcium dynamics
				Ca[m][n] *= exp(-dt/tau_Ca); // === ESSENTIAL ===

				if (delayed_Ca) // if presynaptic spike occurred t_Ca_delay ago
					Ca[m][n] += Ca_pre;

				if (neurons[n].getActivity()) // if postsynaptic spike occurred in previous time step
					Ca[m][n] += Ca_post;

				// E-LTP/-LTD
				if (Ca[m][n] >= theta_p)  // if there is E-LTP
				{
					double noise = sigma_plasticity * sqrt(tau_h) * sqrt(2) * norm_dist(rg) / sqrt(dt); // division by sqrt(dt) was not in Li et al., 2016
					double C = 0.1 + gamma_p + gamma_d;
					double hexp = exp(-dt*C/tau_h);
					h[m][n] = h[m][n] * hexp + (0.1*h_0 + gamma_p + noise) / C * (1.- hexp);
					// [simple Euler: h[m][n] += ((0.1 * (h_0 - h[m][n]) + gamma_p * (1-h[m][n]) - gamma_d * h[m][n] + noise)*(dt/tau_h));]

					if (abs(h[m][n] - h_0) > abs(max_dev))
					{
						max_dev = h[m][n] - h_0;
						tb_max_dev = tb;
					}
				}
				else if (Ca[m][n] >= theta_d) // if there is E-LTD
				{
					double noise = sigma_plasticity * sqrt(tau_h) * norm_dist(rg) / sqrt(dt); // division by sqrt(dt) was not in Li et al., 2016
					double C = 0.1 + gamma_d;
					double hexp = exp(-dt*C/tau_h);
					h[m][n] = h[m][n] * hexp + (0.1*h_0 + noise) / C * (1.- hexp);
					// [simple Euler: h[m][n] += ((0.1 * (h_0 - h[m][n]) + gamma_d * h[m][n] + noise)*(dt/tau_h));]

					if (abs(h[m][n] - h_0) > abs(max_dev))
					{
						max_dev = h[m][n] - h_0;
						tb_max_dev = tb;
					}
				}
				else // if early-phase weight just decays
				{
					double hexp = exp(-dt*0.1/tau_h);
					h[m][n] = h[m][n] * hexp + h_0 * (1.- hexp);
					// [simple Euler: h[m][n] += ((0.1 * (h_0 - h[m][n]))*(dt/tau_h));]
				}

				h_dev = h[m][n] - h_0;
	#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
				if (h_dev > 0.)
					sum_h_diff_p[n] += h_dev; // sum of early-phases changes (for LTP-related protein synthesis)
				else if (h_dev < 0.)
					sum_h_diff_d[n] -= h_dev; // sum of early-phases changes (for LTD-related protein synthesis)
	#endif
	#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
				sum_h_diff[n] += abs(h_dev); // sum of early-phases changes (for protein synthesis)
	#endif
#endif // PLASTICITY == CALCIUM || PLASTICITY == CALCIUM_AND_STC

#if PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

				// L-LTP/-LTD
				if (h_dev >= theta_tag_p) // LTP
				{
	#if PROTEIN_POOLS == POOLS_PCD
					double pa = neurons[n].getPProteinAmount()*neurons[n].getCProteinAmount(); // LTP protein amount times common protein amount from previous time step
	#elif PROTEIN_POOLS == POOLS_PD
					double pa = neurons[n].getPProteinAmount(); // LTP protein amountfrom previous time step
	#elif PROTEIN_POOLS == POOLS_C
					double pa = neurons[n].getCProteinAmount(); // common protein amount from previous time step
	#endif
	#ifdef TWO_NEURONS_ONE_SYNAPSE
					tag_glob = true;
	#endif
					if (pa > EPSILON)
					{
						//double zexp = exp(-dt / (tau_z * 60.) * pa*(step1 + step2));
						double zexp = exp(-dt / (tau_z * 60. * z_max) * pa);
						z[m][n] = z[m][n] * zexp + z_max * (1. - zexp);
						STC = true;
					}
				}
				else if (-h_dev >= theta_tag_d) // LTD
				{
	#if PROTEIN_POOLS == POOLS_PCD
					double pa = neurons[n].getDProteinAmount()*neurons[n].getCProteinAmount(); // LTD protein amount times common protein amount from previous time step
	#elif PROTEIN_POOLS == POOLS_PD
					double pa = neurons[n].getDProteinAmount(); // LTD protein amountfrom previous time step
	#elif PROTEIN_POOLS == POOLS_C
					double pa = neurons[n].getCProteinAmount(); // common protein amount from previous time step
	#endif
	#ifdef TWO_NEURONS_ONE_SYNAPSE
					tag_glob = true;
	#endif
					if (pa > EPSILON)
					{
						//double zexp = exp(-dt / (tau_z * 60.) * pa*(step1 + step2));
						double zexp = exp(-dt / (tau_z * 60.) * pa);
						z[m][n] = z[m][n] * zexp - 0.5 * (1. - zexp);
						STC = true;
					}
				}
				// [simple Euler: z[m][n] += (pa * ((1 - z[m][n]) * step1 - (0.5 + z[m][n]) * step2) * (dt / (tau_z * 6.0 * 10)));]


#endif // PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

#if PLASTICITY == STDP
				double Ap = 1.2, Am = -1.0;
				double tau_syn_stdp = 3e-3;
				double tau_p_stdp = 1e-3;
				double tau_m_stdp = 20e-3;
				double tau_pt_stdp = tau_syn_stdp * tau_p_stdp / (tau_syn_stdp + tau_p_stdp);
				double tau_mt_stdp = tau_syn_stdp * tau_m_stdp / (tau_syn_stdp + tau_m_stdp);
				double eta = 12e-2;

				// if presynaptic neuron m spiked in previous time step
				if (delayed_PSP)
				{
					int last_post_spike = neurons[n].getSpikeHistorySize();

					if (last_post_spike > 0)
					{
						int tb_post = neurons[n].getSpikeTime(last_post_spike);
						double stdp_delta_t = (tb + 1 - tb_post) * dt;

						h[m][n] += eta * exp(-abs(stdp_delta_t) / tau_syn_stdp) * (Ap * (1. + abs(stdp_delta_t)/tau_pt_stdp) + Am * (1. + abs(stdp_delta_t)/tau_mt_stdp));

						if (h[m][n] < 0.)
							h[m][n] = 0.1;//h_0 / 100.; // set to a very small value
					}

				}

				// if postsynaptic neuron n spiked in previous time step
				bool delayed_PSP2 = false;
				for (int k=neurons[n].getSpikeHistorySize(); k>0; k--)
				{
					int st = neurons[n].getSpikeTime(k);

					if (st <= st_PSP) // spikes that have already arrived
					{
						if (st == st_PSP) // if presynaptic spike occurred t_syn_delay ago
						{
							delayed_PSP2 = true; // presynaptic neuron fired t_syn_delay ago
						}
						break;
					}
				}
				if (delayed_PSP2)
				{
					int last_pre_spike = neurons[m].getSpikeHistorySize();

					if (last_pre_spike > 0)
					{
						int tb_pre = neurons[n].getSpikeTime(last_pre_spike);
						double stdp_delta_t = (tb_pre - tb - 1) * dt;

						h[m][n] += eta * (Ap * exp(-abs(stdp_delta_t) / tau_p_stdp) + Am * exp(-abs(stdp_delta_t) / tau_m_stdp));

						if (h[m][n] < 0.)
							h[m][n] = 0.1;//h_0 / 100.; // set to a very small value
					}
				}
#endif // PLASTICITY == STDP

			} // plasticity within excitatory population

#if SYN_SCALING == ON

			h[m][n] += ( eta_ss * pow2(h[m][n] / g_star) * (- r[n]) ) * dt;
#endif

		} // loop over postsynaptic neurons

	} // loop over presynaptic neurons

	return spike_count;
}

/*** processTimeStep_FF ***
 * Processes one timestep for the network only computing late-phase observables [fast-forward mode / compmode == 2] *
 * - int tb: current time step (for printing purposes only) *
 * - double delta_t: duration of the fast-forward timestep *
 * - ofstream* logf: pointer to log file handle (for printing interesting information) *
 * - return: true if late-phase dynamics are persisting, false if not */
int processTimeStep_FF(int tb, double delta_t, ofstream* logf)
{
	bool STC = false; // specifies if at least one synapse is tagged and receives proteins
	bool ps_neuron = false; // specifies if at least one neuron is exhibiting protein synthesis

	/*******************************************************/
	// compute neuronal dynamics
	for (int m=0; m<N; m++) // loop over neurons (in consecutive order)
	{
		// Protein dynamics (for neuron m) - computation from analytic functions
#if PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

	#ifdef TWO_NEURONS_ONE_SYNAPSE
		ps_glob = ps_glob || (sum_h_diff[m] >= theta_pro_c);
	#endif

	#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
		double pa_p = neurons[m].getPProteinAmount();
		double pa_d = neurons[m].getDProteinAmount();
		double p_synth_end_p = 0.;
		double p_synth_end_d = 0.;

		// Potentiation pool
		if (sum_h_diff_p[m] > theta_pro_p) // still rising
		{
			ps_neuron = ps_neuron || true;

			p_synth_end_p = tau_h / 0.1 * log(sum_h_diff_p[m] / theta_pro_p);

			if (delta_t < p_synth_end_p) // rising phase only
			{
				pa_p = pa_p * exp(-delta_t/(tau_pp * 3600.)) + alpha_p * (1. - exp(-delta_t/(tau_pp * 3600.)));
			}
			else // rising phase transitioning to declining phase
			{
				pa_p = pa_p * exp(-p_synth_end_p/(tau_pp * 3600.)) + alpha_p * (1. - exp(-p_synth_end_p/(tau_pp * 3600.)));
				pa_p = pa_p * exp(-(delta_t-p_synth_end_p)/(tau_pp * 3600.));

				*logf << "Protein synthesis (P) ending in neuron " << m << " (t = " << p_synth_end_p + tb*dt << " s)" << endl;
			}
		}
		else // declining phase only
		{
			pa_p = pa_p * exp(-delta_t/(tau_pp * 3600.));
		}

		// Depression pool
		if (sum_h_diff_d[m] > theta_pro_d) // still rising
		{
			ps_neuron = ps_neuron || true;

			p_synth_end_d = tau_h / 0.1 * log(sum_h_diff_d[m] / theta_pro_d);

			if (delta_t < p_synth_end_d) // rising phase only
			{
				pa_d = pa_d * exp(-delta_t/(tau_pd * 3600.)) + alpha_d * (1. - exp(-delta_t/(tau_pd * 3600.)));
			}
			else // rising phase transitioning to declining phase
			{
				pa_d = pa_d * exp(-p_synth_end_d/(tau_pd * 3600.)) + alpha_d * (1. - exp(-p_synth_end_d/(tau_pd * 3600.)));
				pa_d = pa_d * exp(-(delta_t-p_synth_end_d)/(tau_pd * 3600.));

				*logf << "Protein synthesis (D) ending in neuron " << m << " (t = " << p_synth_end_d + tb*dt << " s)" << endl;
			}
		}
		else // declining phase only
		{
			pa_d = pa_d * exp(-delta_t/(tau_pd * 3600.));

		}

		sum_h_diff_p[m] = 0.;
		sum_h_diff_d[m] = 0.;
	#else
		double pa_p = 0.;
		double pa_d = 0.;
	#endif
	#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
		double pa_c = neurons[m].getCProteinAmount();
		double p_synth_end_c = 0.;

		// Common pool
		if (sum_h_diff[m] > theta_pro_c) // still rising
		{
			ps_neuron = ps_neuron || true;

			p_synth_end_c = tau_h / 0.1 * log(sum_h_diff[m] / theta_pro_c);

			if (delta_t < p_synth_end_c) // rising phase only
			{
				pa_c = pa_c * exp(-delta_t/(tau_pc * 3600.)) + alpha_c * (1. - exp(-delta_t/(tau_pc * 3600.)));
			}
			else // rising phase transitioning to declining phase
			{
				pa_c = pa_c * exp(-p_synth_end_c/(tau_pc * 3600.)) + alpha_c * (1. - exp(-p_synth_end_c/(tau_pc * 3600.)));
				pa_c = pa_c * exp(-(delta_t-p_synth_end_c)/(tau_pc * 3600.));

				*logf << "Protein synthesis (C) ending in neuron " << m << " (t = " << p_synth_end_c + tb*dt << " s)" << endl;
			}
		}
		else // declining phase only
		{
			pa_c = pa_c * exp(-delta_t/(tau_pc * 3600.));
		}

		sum_h_diff[m] = 0.;
	#else
		double pa_c = 0.;
	#endif

		neurons[m].setProteinAmounts(pa_p, pa_c, pa_d);

#endif // PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC
	}


	/*******************************************************/
	// compute synaptic dynamics
	for (int m=0; m<N; m++) // loop over presynaptic neurons (in consecutive order)
	{
		for (int in=0; in<neurons[m].getNumberOutgoing(); in++) // loop over postsynaptic neurons
		{
			int n = neurons[m].getOutgoingConnection(in); // get index (in consecutive order) of postsynaptic neuron
			double h_dev; // the deviation of the early-phase weight from its resting state

			// Long-term plasticity
			if (m < pow2(Nl) && n < pow2(Nl)) // plasticity only for exc. -> exc. connections
			{
#if PLASTICITY == CALCIUM || PLASTICITY == CALCIUM_AND_STC

				// early-phase weight just decays
				double hexp = exp(-delta_t*0.1/tau_h);
				h[m][n] = h[m][n] * hexp + h_0 * (1.- hexp);
				// [simple Euler: h[m][n] += ((0.1 * (h_0 - h[m][n]))*(delta_t/tau_h));]

				h_dev = h[m][n] - h_0;
	#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
				if (h_dev > 0.)
					sum_h_diff_p[n] += h_dev; // sum of early-phases changes (for LTP-related protein synthesis)
				else if (h_dev < 0.)
					sum_h_diff_d[n] -= h_dev; // sum of early-phases changes (for LTD-related protein synthesis)
	#endif
	#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
				sum_h_diff[n] += abs(h_dev); // sum of early-phases changes (for protein synthesis)
	#endif
#endif // PLASTICITY == CALCIUM || PLASTICITY == CALCIUM_AND_STC

#if PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

				// L-LTP/-LTD
				if (h_dev >= theta_tag_p) // LTP
				{
	#if PROTEIN_POOLS == POOLS_PCD
					double pa = neurons[n].getPProteinAmount()*neurons[n].getCProteinAmount(); // LTP protein amount times common protein amount from previous time step
	#elif PROTEIN_POOLS == POOLS_PD
					double pa = neurons[n].getPProteinAmount(); // LTP protein amountfrom previous time step
	#elif PROTEIN_POOLS == POOLS_C
					double pa = neurons[n].getCProteinAmount(); // common protein amount from previous time step
	#endif
	#ifdef TWO_NEURONS_ONE_SYNAPSE
					tag_glob = true;
	#endif
					if (pa > EPSILON)
					{
						//double zexp = exp(-delta_t / (tau_z * 60.) * pa*(step1 + step2));
						double zexp = exp(-delta_t / (tau_z * 60. * z_max) * pa);
						z[m][n] = z[m][n] * zexp + z_max * (1. - zexp);
						STC = true;
					}
				}
				else if (-h_dev >= theta_tag_d) // LTD
				{
	#if PROTEIN_POOLS == POOLS_PCD
					double pa = neurons[n].getDProteinAmount()*neurons[n].getCProteinAmount(); // LTD protein amount times common protein amount from previous time step
	#elif PROTEIN_POOLS == POOLS_PD
					double pa = neurons[n].getDProteinAmount(); // LTD protein amountfrom previous time step
	#elif PROTEIN_POOLS == POOLS_C
					double pa = neurons[n].getCProteinAmount(); // common protein amount from previous time step
	#endif
	#ifdef TWO_NEURONS_ONE_SYNAPSE
					tag_glob = true;
	#endif
					if (pa > EPSILON)
					{
						//double zexp = exp(-delta_t / (tau_z * 60.) * pa*(step1 + step2));
						double zexp = exp(-delta_t / (tau_z * 60.) * pa);
						z[m][n] = z[m][n] * zexp - 0.5 * (1. - zexp);
						STC = true;
					}
				}
				// [simple Euler: z[m][n] += (pa * ((1 - z[m][n]) * step1 - (0.5 + z[m][n]) * step2) * (delta_t / (tau_z * 6.0 * 10)));]


#endif // PLASTICITY == CALCIUM_AND_STC || PLASTICITY == STDP_AND_STC

			} // plasticity within excitatory population

		} // loop over postsynaptic neurons

	} // loop over presynaptic neurons

	if (!STC) // no late-phase dynamics can take place anymore
	{
		return false;
	}

	return true;
}

/*** getSumDiff ***
 * Returns the sum of absolute values of weight differences to the initial value for a specific neuron *
 * - n: number of the neuron to be considered *
 * - return: sum_h_diff for a specific neuron */
double getSumDiff(int n)
{
	return sum_h_diff[n];
}

#ifdef TWO_NEURONS_ONE_SYNAPSE
/*** getPlasticityType ***
 * Returns the kind of plasticity evoked by the stimulus *
 * - return: 0 for ELTP, 1 for ELTP with tag, 2 for LLTP, 3 for ELTD, 4 for ELTD with tag, 5 for LLTD, -1 else */
int getPlasticityType()
{
	if (max_dev > EPSILON) // LTP
	{
		if (!tag_glob)
			return 0;
		else
		{
			if (!ps_glob)
				return 1;
			else
				return 2;
		}
	}
	else if (max_dev < -EPSILON) // LTD
	{
		if (!tag_glob)
			return 3;
		else
		{
			if (!ps_glob)
				return 4;
			else
				return 5;
		}
	}

	return -1;
}

/*** getMaxDev ***
 * Returns the maximum deviation from h_0 *
 * - return: max_dev*/
double getMaxDev()
{
	return max_dev;
}
#endif

/*** getTagVanishTime ***
 * Returns the time by which all tags will have vanished, based on *
 * the largest early-phase deviation from the mean h_0 (max_dev) *
 * - return: the time difference */
double getTagVanishTime()
{
	double tag_vanish = tau_h / 0.1 * log(abs(max_dev) / min(theta_tag_p, theta_tag_d));

	if (abs(max_dev) > EPSILON && tag_vanish > EPSILON)
		return tag_vanish + tb_max_dev*dt;
	else
		return 0.;
}

/*** getProteinSynthesisEnd ***
 * Returns the time (for every pool) by which all protein synthesis will halt, based on the *
 * largest sum of early-phase deviations from the mean h_0 (max_sum_diff*) *
 * - return: the times for the different pools (P,C,D) in a vector */
vector<double> getProteinSynthesisEnd()
{
	vector<double> ret(3,0.);

#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
	if (max_sum_diff > theta_pro_c)
		ret[1] = tau_h / 0.1 * log(max_sum_diff / theta_pro_c) + tb_max_sum_diff*dt;
#endif

#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
	if (max_sum_diff_p > theta_pro_p)
		ret[0] = tau_h / 0.1 * log(max_sum_diff_p / theta_pro_p) + tb_max_sum_diff_p*dt;

	if (max_sum_diff_d > theta_pro_d)
		ret[2] = tau_h / 0.1 * log(max_sum_diff_d / theta_pro_d) + tb_max_sum_diff_d*dt;
#endif
	return ret;
}

/*** getThreshold ***
 * Returns a specified threshold value (for tag or protein synthesis) *
 * - plast: the type of plasticity (1: LTP, 2: LTD)
 * - which: the type of threshold (1: early-phase calcium treshold, 2: tagging threshold, 3: protein synthesis threshold)
 * - return: the threshold value */
double getThreshold(int plast, int which)
{
	if (plast == 1) // LTP
	{
		if (which == 1) // early-phase calcium treshold
			return theta_p;
		else if (which == 2) // tagging threshold
			return theta_tag_p;
		else // protein synthesis threshold
			return theta_pro_p;
	}
	else // LTD
	{
		if (which == 1) // early-phase calcium treshold
			return theta_d;
		else if (which == 2) // tagging threshold
			return theta_tag_d;
		else // protein synthesis threshold
			return theta_pro_d;
	}
}


/*** setRhombStimulus ***
 * Sets a spatially rhomb-shaped firing rate stimulus in the exc. population *
 * - Stimulus& _st: shape of one stimulus period *
 * - int center: the index of the neuron in the center of the rhomb *
 * - int radius: the "radius" of the rhomb in neurons (radius 3 corresponds to a rhomb containing 25 neurons, radius 4 to 41, radius 5 to 61, radius 9 to 181) */
void setRhombStimulus(Stimulus& _st, int center, int radius)
{
	for (int i=-radius; i<=radius; i++)
	{
		int num_cols = (radius-abs(i));

		for (int j=-num_cols; j<=num_cols; j++)
		{
			neurons[center+i*Nl+j].setCurrentStimulus(_st); // set temporal course of current stimulus for given neuron
		}
	}

	setStimulationEnd(_st.getStimulationEnd());
}

/*** setRhombPartialRandomStimulus ***
 * Sets stimulation for randomly drawn neurons out of a rhomb shape in the exc. population *
 * - Stimulus& _st: shape of one stimulus period *
 * - int center: the index of the neuron in the center of the rhomb *
 * - int radius: the "radius" of the rhomb in neurons (radius 3 corresponds to a rhomb containing 25 neurons) *
 * - double fraction: the fraction of neurons in the rhomb that shall be stimulated */
void setRhombPartialRandomStimulus(Stimulus& _st, int center, int radius, double fraction)
{
	int total_assembly_size = 2*pow2(radius) + 2*radius + 1; // total number of neurons within the rhomb
	int ind = 0, count = 0;
	uniform_int_distribution<int> rhomb_dist(1, total_assembly_size);
	int* indices; // array of indices (consectuive neuron numbers) for rhomb neurons

	indices = new int[total_assembly_size];

	// gather indices of neurons belonging to the rhomb
	for (int i=-radius; i<=radius; i++)
	{
		int num_cols = (radius-abs(i));

		for (int j=-num_cols; j<=num_cols; j++)
		{
			indices[ind++] = center+i*Nl+j;
		}
	}

	// draw random neurons out of the rhomb
	while(count < fraction*total_assembly_size)
	{
		int chosen_n = rhomb_dist(rg);

		if (indices[chosen_n-1] >= 0) // found a neuron that has not be assigned a stimulus
		{
			neurons[indices[chosen_n-1]].setCurrentStimulus(_st); // set temporal course of current stimulus for given neuron
			count++;
			indices[chosen_n-1] = -1;
		}
	}

	setStimulationEnd(_st.getStimulationEnd());

	delete[] indices;
}

/*** setRhombPartialStimulus ***
 * Sets stimulation for first fraction of neurons out of a rhomb shape in the exc. population *
 * - Stimulus& _st: shape of one stimulus period *
 * - int center: the index of the neuron in the center of the rhomb *
 * - int radius: the "radius" of the rhomb in neurons (radius 3 corresponds to a rhomb containing 25 neurons) *
 * - double fraction: the fraction of neurons in the rhomb that shall be stimulated */
void setRhombPartialStimulus(Stimulus& _st, int center, int radius, double fraction)
{
	int count = int(fraction * (2*radius*(radius+1)+1));

	for (int i=-radius; i<=radius; i++)
	{
		int num_cols = (radius-abs(i));

		for (int j=-num_cols; j<=num_cols; j++)
		{
			neurons[center+i*Nl+j].setCurrentStimulus(_st); // set temporal course of current stimulus for given neuron
			count--;
			if (count == 0)
				break;
		}
		if (count == 0)
			break;
	}

	setStimulationEnd(_st.getStimulationEnd());
}

/*** setRandomStimulus ***
 * Sets a randomly distributed firing rate stimulus in the exc. population *
 * - Stimulus& _st: shape of one stimulus period *
 * - int num: the number of neurons to be drawn (i.e., to be stimulated) *
 * - ofstream* f [optional]: handle to a file for output of the randomly drawn neuron numbers *
 * - int range_start [optional]: the lowest neuron number that can be drawn *
 * - int range_end [optional]: one plus the highest neuron number that can be drawn (-1: highest possible) */
void setRandomStimulus(Stimulus& _st, int num, ofstream* f = NULL, int range_start=0, int range_end=-1)
{
	int range_len = (range_end == -1) ? (pow2(Nl) - range_start) : (range_end - range_start); // the number of neurons eligible for being drawn
	bool* stim_neurons = new bool [range_len];
	uniform_int_distribution<int> u_dist_neurons(0, range_len-1); // uniform distribution to draw neuron numbers
	int neurons_left = num;

	if (f != NULL)
		*f << "Randomly drawn neurons for stimulation:" << "   { ";

	for (int i=0; i<range_len; i++)
		stim_neurons[i] = false;

	while (neurons_left > 0)
	{
		int neur = u_dist_neurons(rg); // draw a neuron

		if (!stim_neurons[neur]) // if the stimulus has not yet been assigned to the drawn neuron
		{
			neurons[neur+range_start].setCurrentStimulus(_st); // set temporal course of current stimulus for drawn neuron
			stim_neurons[neur] = true;
			neurons_left--;
			if (f != NULL)
				*f << neur+range_start << ", "; // print stimulated neuron to file
		}
	}

	setStimulationEnd(_st.getStimulationEnd());

	if (f != NULL)
		*f << " }" << endl;
	delete[] stim_neurons;
}

/*** setSingleNeuronStimulus ***
 * Sets a firing rate stimulus for a specified neuron *
 * - int m: number of the neuron to be stimulated *
 * - Stimulus& _st: shape of one stimulus period */
void setSingleNeuronStimulus(int m, Stimulus& _st)
{
	neurons[m].setCurrentStimulus(_st);

	setStimulationEnd(_st.getStimulationEnd());
}

/*** setBlockStimulus ***
 * Sets a stimulus for a given block of n neurons in the network *
 * - Stimulus& _st: shape of one stimulus period *
 * - int n: the number of neurons that shall be stimulated *
 * - int off [optional]: the offset that defines at which neuron number the block begins */
void setBlockStimulus(Stimulus& _st, int n, int off=0)
{
	for (int i=off; i<(n+off); i++)
	{
		neurons[i].setCurrentStimulus(_st); // set temporal course of current stimulus for given neuron
	}

	setStimulationEnd(_st.getStimulationEnd());
}

/*** setConstCurrent ***
 * Sets the constant current for all neurons to a newly defined value
 * - double _I_const: constant current in nA */
void setConstCurrent(double _I_const)
{
	for (int m=0; m<N; m++)
	{
		neurons[m].setConstCurrent(_I_const);
	}
}

#if STIPULATE_CA == ON
/*** stipulateRhombAssembly ***
 * Stipulates a rhomb-shaped cell assembly with strong interconnections *
 * - int center: the index of the neuron in the center of the rhomb *
 * - int radius: the "radius" of the rhomb in neurons (radius 3 corresponds to a rhomb containing 25 neurons, radius 5 to 61 neurons) */
void stipulateRhombAssembly(int center, int radius)
{
	double value = 2*h_0; // stipulated value

	for (int i=-radius; i<=radius; i++)
	{
		int num_cols = (radius-abs(i));

		for (int j=-num_cols; j<=num_cols; j++)
		{
			int m = center+i*Nl+j;

			for (int k=-radius; k<=radius; k++)
			{
				int num_cols = (radius-abs(k));

				for (int l=-num_cols; l<=num_cols; l++)
				{
					int n = center+k*Nl+l;

					if (conn[m][n]) // set all connections within the assembly to this value
						h[m][n] = value;
				}
			}
		}
	}
}


/*** stipulateFirstNeuronsAssembly ***
 * Stipulates an cell assembly consisting of the first neurons with strong interconnections *
 * - int n: the number of neurons that shall be stimulated */
void stipulateFirstNeuronsAssembly(int n)
{
	double value = 2*h_0; // stipulated value

	for (int i=0; i<n; i++)
	{
		for (int j=0; j<n; j++)
		{
			if (conn[i][j]) // set all connections within the assembly to this value
				h[i][j] = value;
		}
	}
}
#endif

/*** setSigma ***
 * Sets the standard deviation for the external input current for all neurons to a newly defined value
 * - double _sigma: standard deviation in nA s^(1/2) */
void setSigma(double _sigma)
{
	for (int m=0; m<N; m++)
	{
		neurons[m].setSigma(_sigma);
	}
}

/*** setSynTimeConstant ***
 * Sets the synaptic time constant *
 * - double _tau_syn: synaptic time constant in s */
void setSynTimeConstant(double _tau_syn)
{
	tau_syn = _tau_syn;
	for (int m=0; m<N; m++)
	{
		neurons[m].setTauOU(tau_syn);
	}
}

/*** setCouplingStrengths ***
 * Sets the synaptic coupling strengths *
 * - double _w_ee: coupling strength for exc. -> exc. connections in units of h_0 *
 * - double _w_ei: coupling strength for exc. -> inh. connections in units of h_0 *
 * - double _w_ie: coupling strength for inh. -> exc. connections in units of h_0 *
 * - double _w_ii: coupling strength for inh. -> inh. connections in units of h_0 */
void setCouplingStrengths(double _w_ee, double _w_ei, double _w_ie, double _w_ii)
{
	w_ee = _w_ee * h_0;
	w_ei = _w_ei * h_0;
	w_ie = _w_ie * h_0;
	w_ii = _w_ii * h_0;
}

/*** getInitialWeight ***
 * Returns the initial exc.->exc. weight (typically h_0) */
double getInitialWeight()
{
	return w_ee;
}

/*** getSynapticCalcium ***
 * Returns the calcium amount at a given synapse *
 * - synapse s: structure specifying pre- and postsynaptic neuron *
 * - return: the synaptic calcium amount */
double getSynapticCalcium(synapse s) const
{
	return Ca[s.presyn_neuron][s.postsyn_neuron];
}

/*** getEarlySynapticStrength ***
 * Returns the early-phase synaptic strength at a given synapse *
 * - synapse s: structure specifying pre- and postsynaptic neuron *
 * - return: the early-phase synaptic strength */
double getEarlySynapticStrength(synapse s) const
{
	return h[s.presyn_neuron][s.postsyn_neuron];
}

/*** getLateSynapticStrength ***
 * Returns the late-phase synaptic strength at a given synapse *
 * - synapse s: structure specifying pre- and postsynaptic neuron *
 * - return: the late-phase synaptic strength */
double getLateSynapticStrength(synapse s) const
{
	return z[s.presyn_neuron][s.postsyn_neuron];
}

/*** getMeanEarlySynapticStrength ***
 * Returns the mean early-phase synaptic strength (averaged over all synapses within the given set of neurons) *
 * - int n: the number of neurons that shall be considered (e.g., n=Nl^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the mean early-phase synaptic strength */
double getMeanEarlySynapticStrength(int n, int off=0) const
{
	double h_mean = 0.;
	int c_number = 0;

	for (int i=off; i < (n+off); i++)
	{
		for (int j=off; j < (n+off); j++)
		{
			if (conn[i][j])
			{
				c_number++;
				h_mean += h[i][j];
			}
			
		}	
	}

	h_mean /= c_number;

	return h_mean;
}

/*** getMeanLateSynapticStrength ***
 * Returns the mean late-phase synaptic strength (averaged over all synapses within the given set of neurons) *
 * - int n: the number of neurons that shall be considered (e.g., n=Nl^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the mean late-phase synaptic strength */
double getMeanLateSynapticStrength(int n, int off=0) const
{
	double z_mean = 0.;
	int c_number = 0;

	for (int i=off; i < (n+off); i++)
	{
		for (int j=off; j < (n+off); j++)
		{
			if (conn[i][j])
			{
				c_number++;
				z_mean += z[i][j];
			}
			
		}
	}

	z_mean /= c_number;

	return z_mean;
}

/*** getSDEarlySynapticStrength ***
 * Returns the standard deviation of the early-phase synaptic strength (over all synapses within the given set of neurons) *
 * - double mean: the mean of the early-phase syn. strength within the given set
 * - int n: the number of neurons that shall be considered (e.g., n=Nl^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the std. dev. of the early-phase synaptic strength */
double getSDEarlySynapticStrength(double mean, int n, int off=0) const
{
	double h_sd = 0.;
	int c_number = 0;

	for (int i=off; i < (n+off); i++)
	{
		for (int j=off; j < (n+off); j++)
		{
			if (conn[i][j])
			{
				c_number++;
				h_sd += pow2(h[i][j] - mean);
			}
			
		}
	}

	h_sd = sqrt(h_sd / c_number);

	return h_sd;
}

/*** getSDLateSynapticStrength ***
 * Returns the standard deviation of the late-phase synaptic strength (over all synapses within the given set of neurons) *
 * - double mean: the mean of the late-phase syn. strength within the given set
 * - int n: the number of neurons that shall be considered (e.g., n=Nl^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the std. dev. of the late-phase synaptic strength */
double getSDLateSynapticStrength(double mean, int n, int off=0) const
{
	double z_sd = 0.;
	int c_number = 0;

	for (int i=off; i < (n+off); i++)
	{
		for (int j=off; j < (n+off); j++)
		{
			if (conn[i][j])
			{
				c_number++;
				z_sd += pow2(z[i][j] - mean);
			}
			
		}	
	}

	z_sd = sqrt(z_sd / c_number);

	return z_sd;
}

/*** getMeanCProteinAmount ***
 * Returns the mean protein amount (averaged over all neurons within the given set) *
 * - int n: the number of neurons that shall be considered (e.g., n=Nl^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the mean protein amount */
double getMeanCProteinAmount(int n, int off=0) const
{
	double pa = 0.;

	for (int i=off; i < (n+off); i++)
	{
		pa += neurons[i].getCProteinAmount();
	}

	pa /= n;

	return pa;
}

/*** getSDCProteinAmount ***
 * Returns the standard deviation of the protein amount (over all neurons within the given set) *
 * - double mean: the mean of the protein amount within the given set
 * - int n: the number of neurons that shall be considered (e.g., n=Nl_exc^2 for all excitatory neurons, or n=N for all neurons) *
 * - int off [optional]: the offset that defines at which neuron number the considered range begins *
 * - return: the std. dev. of the protein amount */
double getSDCProteinAmount(double mean, int n, int off=0) const
{
	double pa_sd = 0.;

	for (int i=off; i < (n+off); i++)
	{
		pa_sd += pow2(neurons[i].getCProteinAmount() - mean);
	}

	pa_sd = sqrt(pa_sd / n);

	return pa_sd;
}


/*** readConnections ***
 * Reads the connectivity matrix from a file, either given by a text-converted numpy array or a plain matrix structure *
 * - file: a text file where the matrix is located *
 * - format: 0 for plain matrix structure, 1 for numpy structure with square brackets
 * - return: 2 - successful, 1 - file could not opened, 0 - dimension mismatch */
int readConnections(string file, int format = 0)
{
	ifstream f(file, ios::in); // the file handle
	string buf; // buffer to read one line

	int m, n = 0; // pre- and postsynaptic neuron
	int initial_brackets = 0; // specifies if the initial brackets have been read yet

	if (!f.is_open()) // check if file was opened successfully
		return 1;

	for (int a=0;a<N;a++) // reset connections of all neurons
		neurons[a].resetConnections();

	if (format == 0)
		m = -1;
	else
		m = 0;

	while (getline(f, buf)) // while end of file has not yet been reached
	{
		if (format == 0 && buf.size() > 0)
		{
			m++;
			n = 0;
		}

		for (int i=0; i<buf.size(); i++) // go through all characters in this line
		{
			if (format == 1 && buf[i] == '[')
			{
				if (initial_brackets < 2) // still reading the initial brackets
					initial_brackets++;
				else
				{
					m++;
					n = 0;
				}
			}
			else if (buf[i] == '1')
			{
				if (m < pow2(Nl) && n < pow2(Nl)) // exc. -> exc.
				{
					neurons[m].addOutgoingConnection(n, TYPE_EXC);
					//cout << m << " -> " << n << " added" << endl;
					neurons[n].incNumberIncoming(TYPE_EXC);
				}
				else if (m < pow2(Nl) && n >= pow2(Nl)) // exc. -> inh.
				{
					neurons[m].addOutgoingConnection(n, TYPE_INH);
					//cout << m << " -> " << n << " added" << endl;
					neurons[n].incNumberIncoming(TYPE_EXC);
				}
				else if (m >= pow2(Nl) && n < pow2(Nl)) // inh. -> exc.
				{
					neurons[m].addOutgoingConnection(n, TYPE_EXC);
					//cout << m << " -> " << n << " added" << endl;
					neurons[n].incNumberIncoming(TYPE_INH);
				}
				else if (m >= pow2(Nl) && n >= pow2(Nl)) // inh. -> inh.
				{
					neurons[m].addOutgoingConnection(n, TYPE_INH);
					//cout << m << " -> " << n << " added" << endl;
					neurons[n].incNumberIncoming(TYPE_INH);
				}

				conn[m][n] = true;
				n++;
			}
			else if (buf[i] == '0')
			{
				conn[m][n] = false;
				n++;
			}
		}
	}
	f.close();

	reset();

	if (m != (N-1) || n != N) // if dimensions do not match
		return 0;

	return 2;
}

/*** printConnections ***
 * Prints the connection matrix to a given file (either in numpy structure or in plain matrix structure) *
 * - file: a text file where the matrix is located *
 * - format: 0 for plain matrix structure, 1 for numpy structure with square brackets
 * - return: 2 - successful, 1 - file could not opened */
int printConnections(string file, int format = 0)
{
	ofstream f(file); // the file handle

	if (!f.is_open()) // check if file was opened successfully
		return 1;

	if (format == 1)
		f << "[";

	for (int m=0;m<N;m++)
	{
		if (format == 1)
			f << "[";
		for (int n=0;n<N;n++)
		{
			if (conn[m][n])
				f << "1 ";
			else
				f << "0 ";
		}
		if (format == 1)
			f << "]";
		f << endl;
	}

	if (format == 1)
		f << "]";
	f.close();

	return 2;
}

/*** printConnections2 ***
 * FOR TESTING THE CONNECTIONS SAVED IN ARRAYS IN NEURONS *
 * - file: a text file where the matrix is located *
 * - format: 0 for plain matrix structure, 1 for numpy structure with square brackets
 * - return: 2 - successful, 1 - file could not opened *
int printConnections2(string file, int format = 0)
{
	ofstream f(file); // the file handle

	if (!f.is_open()) // check if file was opened successfully
		return 1;

	if (format == 1)
		f << "[";

	for (int m=0;m<N;m++)
	{
		if (format == 1)
			f << "[";

		int in=0;
		for (int n=0;n<N;n++)
		{
			if (conn[m][n] && neurons[m].getNumberOutgoing() > in && neurons[m].getOutgoingConnection(in) == n)
			{
				f << "1 ";
				in++;
			}
			else
				f << "0 ";

		}
		if (format == 1)
			f << "]";
		f << endl;
	}

	if (format == 1)
		f << "]";
	f.close();

	return 2;
}*/


/*** printAllInitialWeights ***
 * Prints the connection matrix to a given file (either in numpy structure or in plain matrix structure) *
 * - file: a text file where the matrix is located *
 * - format: 0 for plain matrix structure, 1 for numpy structure with square brackets
 * - return: 2 - successful, 1 - file could not opened */
int printAllInitialWeights(string file, int format = 0)
{
	ofstream f(file); // the file handle

	if (!f.is_open()) // check if file was opened successfully
		return 1;

	if (format == 1)
		f << "[";

	for (int m=0;m<N;m++)
	{
		if (format == 1)
			f << "[";
		for (int n=0;n<N;n++)
		{
			// Output of all initial weights
			if (conn[m][n])
			{
				if (m < pow2(Nl) && n < pow2(Nl)) // exc. -> exc.
					f << h[m][n] << " ";
				else if (m < pow2(Nl) && n >= pow2(Nl)) // exc. -> inh.
					f << w_ei << " ";
				else if (m >= pow2(Nl) && n < pow2(Nl)) // inh. -> exc.
					f << w_ie << " ";
				else if (m >= pow2(Nl) && n >= pow2(Nl)) // inh. -> inh.
					f << w_ii << " ";

			}
			else
				f << 0. << " ";
		}
		if (format == 1)
			f << "]";
		f << endl;
	}

	if (format == 1)
		f << "]";
	f.close();

	return 2;
}

/*** readCouplingStrengths ***
 * Reads the excitatory coupling strengths from a text file that contains a matrix for early-phase weights and a matrix *
 * for late-phase weights, each terminated by a blank line *
 * - file: a text file where the matrix is located *
 * - return: 2 - successful, 1 - file could not opened, 0 - dimension mismatch */
int readCouplingStrengths(string file)
{
	int phase = 1; // 1: reading early-phase values, 2: reading late-phase values
	int m = 0, n; // pre- and postsynaptic neuron
	double strength;
	ifstream f(file, ios::in); // the file handle
	string buf; // buffer to read one line

	if (!f.is_open())
		return 1;

	// Read early- and late-phase matrix
	while (getline (f, buf))
    	{
		istringstream iss(buf);

		if (!buf.empty())
		{
			n = 0;
			while(iss >> strength)
			{
				if (phase == 1)
					h[m][n] = strength;
				else
					z[m][n] = strength;

				n++;
			}
			m++;
		}
		else // blank line encountered
		{
			if (phase == 1) // now begins the second phase
			{
				if (m != pow2(Nl) || n != pow2(Nl)) // if dimensions do not match
				{
					f.close();
					return 0;
				}

				phase = 2;
				m = 0;
				n = 0;
			}
			else
				break;
		}
    	}

	f.close();

	if (m != pow2(Nl) || n != pow2(Nl)) // if dimensions do not match
	{
		return 0;
	}

	return 2;
}

/*** setStimulationEnd ***
 * Tells the Network instance the end of stimulation (even if not all stimuli are yet set) *
 * - int stim_end: the time step in which stimulation ends */
void setStimulationEnd(int stim_end)
{
	if (stim_end > stimulation_end)
		stimulation_end = stim_end;
}


/*** setSpikeStorageTime ***
 * Sets the number of timesteps for which spikes have to be kept *
 * - int storage_steps: the size of the storage timespan in timesteps */
void setSpikeStorageTime(int storage_steps)
{
	for (int m=0; m<N; m++)
	{
		neurons[m].setSpikeHistoryMemory(storage_steps);
	}

}

/*** resetLastSpikeIndex ***
 * Resets the last spike index of a neuron important to its calcium dynamics *
 * - int m: the neuron number */
void resetLastSpikeIndex(int m)
{
	last_Ca_spike_index[m] = 1;
}

/*** resetPlasticity ***
 * Depending on the arguments, undoes plastic changes that the network has undergone, *
 * resets calcium values or protein values *
 * - bool early_phase: resets all early-phase weights and calcium concentrations *
 * - bool late_phase: resets all late-phase weights *
 * - bool calcium: resets all calcium concentrations *
 * - bool proteins: resets all neuronal protein pools */
void resetPlasticity(bool early_phase, bool late_phase, bool calcium, bool proteins)
{
	for (int m=0; m<N; m++)
	{
		if (proteins)
			neurons[m].setProteinAmounts(0., 0., 0.);

		for (int n=0; n<N; n++) // reset synapses
		{
			if (early_phase)
			{
				if (conn[m][n])
					h[m][n] = h_0;
				else
					h[m][n] = 0.;
			}
			if (calcium)
				Ca[m][n] = 0.;
			if (late_phase)
				z[m][n] = 0.;
		}
	}
}

/*** reset ***
 * Resets the network and all neurons to initial state (but maintain connectivity) */
void reset()
{
	rg.seed(getClockSeed()); // set new seed by clock's epoch
	u_dist.reset(); // reset the uniform distribution for random numbers
	norm_dist.reset(); // reset the normal distribution for random numbers

	for (int m=0; m<N; m++)
	{
		neurons[m].reset();
		for (int n=0; n<N; n++) // reset synapses
		{
			if (conn[m][n])
				h[m][n] = h_0;
			else
				h[m][n] = 0.;
			Ca[m][n] = 0.;
			z[m][n] = 0.;
		}
		resetLastSpikeIndex(m);
		sum_h_diff[m] = 0.;
		sum_h_diff_p[m] = 0.;
		sum_h_diff_d[m] = 0.;
	}
	stimulation_end = 0;
	max_dev = 0.;
	tb_max_dev = 0;
#if PROTEIN_POOLS == POOLS_C || PROTEIN_POOLS == POOLS_PCD
	max_sum_diff = 0.;
	tb_max_sum_diff = 0;
#endif
#if PROTEIN_POOLS == POOLS_PD || PROTEIN_POOLS == POOLS_PCD
	max_sum_diff_p = 0.;
	tb_max_sum_diff_p = 0;
	max_sum_diff_d = 0.;
	tb_max_sum_diff_d = 0;
#endif

#ifdef TWO_NEURONS_ONE_SYNAPSE
	tag_glob = false;
	ps_glob = false;
#endif
}

/*** setCaConstants ***
 * Set constants for the calcium dynamics *
 * - double _theta_p: the potentiation threshold *
 * - double _theta_d: the potentiation threshold */
void setCaConstants(double _theta_p, double _theta_d, double _Ca_pre, double _Ca_post)
{
	theta_p = _theta_p;
	theta_d = _theta_d;
	Ca_pre = _Ca_pre;
	Ca_post = _Ca_post;
}

/*** setPSThresholds ***
 * Set thresholds for the onset of protein synthesis *
 * - double _theta_pro_P: the threshold for P synthesis (in units of h0) *
 * - double _theta_pro_C: the threshold for C synthesis (in units of h0) *
 * - double _theta_pro_D: the threshold for D synthesis (in units of h0) */
void setPSThresholds(double _theta_pro_P, double _theta_pro_C, double _theta_pro_D)
{
	theta_pro_p = _theta_pro_P*h_0;
	theta_pro_c = _theta_pro_C*h_0;
	theta_pro_d = _theta_pro_D*h_0;
}

/*** Constructor ***
 * Sets all parameters, creates neurons and synapses *
 * --> it is required to call setSynTimeConstant and setCouplingStrengths immediately *
 *     after calling this constructor! *
 * - double _dt: the length of one time step in s *
 * - int _Nl: the number of neurons in one line in excitatory population (row/column) *
 * - int _Nl_inh: the number of neurons in one line in inhibitory population (row/column) - line structure so that stimulation of inhib. *
                  population could be implemented more easily *
 * - double _p_c: connection probability *
 * - double _sigma_plasticity: standard deviation of the plasticity *
 * - double _z_max: the upper z bound */

Network(const double _dt, const int _Nl, const int _Nl_inh, double _p_c, double _sigma_plasticity, double _z_max) :
        dt(_dt), rg(getClockSeed()), u_dist(0.0,1.0), norm_dist(0.0,1.0), Nl(_Nl), Nl_inh(_Nl_inh), z_max(_z_max)
{
	N = pow2(Nl) + pow2(Nl_inh); // total number of neurons

	p_c = _p_c; // set connection probability

	t_syn_delay = 0.003; // from https://www.britannica.com/science/nervous-system/The-neuronal-membrane#ref606406, accessed 18-06-21
#if defined TWO_NEURONS_ONE_SYNAPSE && !defined TWO_NEURONS_ONE_SYNAPSE_ALT
	t_syn_delay = dt;
#endif

	t_syn_delay_steps = int(t_syn_delay/dt);

	// Biophysical parameters for stimulation, Ca dynamics and early phase
	t_Ca_delay = 0.0188; // from Graupner and Brunel (2012), hippocampal slices
	t_Ca_delay_steps = int(t_Ca_delay/dt);
	Ca_pre = 1.0; // from Graupner and Brunel (2012), hippocampal slices
	Ca_post = 0.2758; // from Graupner and Brunel (2012), hippocampal slices
	tau_Ca = 0.0488; // from Graupner and Brunel (2012), hippocampal slices
	tau_Ca_steps = int(tau_Ca/dt);
	tau_h = 688.4; // from Graupner and Brunel (2012), hippocampal slices
	gamma_p = 1645.6; // from Graupner and Brunel (2012), hippocampal slices
	gamma_d = 313.1; // from Graupner and Brunel (2012), hippocampal slices
	h_0 = 0.5*(gamma_p/(gamma_p+gamma_d)); // from Li et al. (2016)
	theta_p = 3.0; // from Li et al. (2016)
	theta_d = 1.2; // from Li et al. (2016)
	sigma_plasticity = _sigma_plasticity; // from Graupner and Brunel (2012) but corrected by 1/sqrt(1000)

	// Biophysical parameters for protein synthesis and late phase
	tau_pp = 1.0; // from Li et al. (2016)
	tau_pc = 1.0; // from Li et al. (2016)
	tau_pd = 1.0; // from Li et al. (2016)
	alpha_p = 1.0; // from Li et al. (2016)
	alpha_c = 1.0; // from Li et al. (2016)
	alpha_d = 1.0; // from Li et al. (2016)
	tau_z = 60.0; // from Li et al. (2016) - includes "gamma"
	theta_pro_p = 0.5*h_0; //
	theta_pro_c = 0.5*h_0; // from Li et al. (2016)
	theta_pro_d = 0.5*h_0; //
	theta_tag_p = 0.2*h_0; // from Li et al. (2016)
	theta_tag_d = 0.2*h_0; // from Li et al. (2016)

	// Create neurons and synapse matrices
	neurons = vector<Neuron> (N, Neuron(_dt));
	conn = new bool* [N];
	Ca = new double* [N];
	h = new double* [N];
	z = new double* [N];
	last_Ca_spike_index = new int [N];
	sum_h_diff = new double [N];
	sum_h_diff_p = new double [N];
	sum_h_diff_d = new double [N];

	for (int m=0; m<N; m++)
	{
		if (m < pow2(Nl)) // first Nl^2 neurons are excitatory
			neurons[m].setType(TYPE_EXC);
		else // remaining neurons are inhibitory
			neurons[m].setType(TYPE_INH);
		conn[m] = new bool [N];
		Ca[m] = new double [N];
		h[m] = new double [N];
		z[m] = new double [N];

		// create synaptic connections
		for (int n=0; n<N; n++)
		{
			conn[m][n] = false; // necessary for resetting the connections

			if (m != n) // if not on main diagonal (which should remain zero)
			{
				conn[m][n] = shallBeConnected(m, n); // use random generator depending on connection probability

			}
		}
	}

#ifdef TWO_NEURONS_ONE_SYNAPSE
	neurons[1].setPoisson(true);
	#ifdef PLASTICITY_OVER_FREQ
	neurons[0].setPoisson(true);
	#endif
#endif
}

/*** Destructor ***
 * Cleans up the allocated memory for arrays */
~Network()
{
	for(int i=0; i<N; i++)
	{
		delete[] conn[i];
		delete[] Ca[i];
		delete[] h[i];
		delete[] z[i];
	}

	delete[] conn;
	delete[] Ca;
	delete[] h;
	delete[] z;
	delete[] last_Ca_spike_index;
	delete[] sum_h_diff;
	delete[] sum_h_diff_p;
	delete[] sum_h_diff_d;
}

/* =============================================================================================================================== */
/* ==== Functions redirecting to corresponding functions in Neuron class ========================================================= */

/* Two versions are given for each function, one for consecutive and one for row/column numbering */

/*** getType ***
 * Returns the type of neuron (i|j) (inhbitory/excitatory) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the neuron type */
int getType(int i, int j) const
{
	return neurons[cNN(i,j)].getType();
}
int getType(int m) const
{
	return neurons[m].getType();
}

/*** getVoltage ***
 * Returns the membrane potential of neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the membrane potential in mV */
double getVoltage(int i, int j) const
{
	return neurons[cNN(i,j)].getVoltage();
}
double getVoltage(int m) const
{
	return neurons[m].getVoltage();
}

/*** getThreshold ***
 * Returns the value of the dynamic membrane threshold of neuron (i|j) *
 * - return: the membrane threshold in mV */
double getThreshold(int i, int j) const
{
	return neurons[cNN(i,j)].getThreshold();
}
double getThreshold(int m) const
{
	return neurons[m].getThreshold();
}

/*** getCurrent ***
 * Returns total current effecting neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous current in nA */
double getCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getCurrent();
}
double getCurrent(int m) const
{
	return neurons[m].getCurrent();
}

/*** getStimulusCurrent ***
 * Returns current evoked by external stimulation in neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous current stimulus in nA */
double getStimulusCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getStimulusCurrent();
}
double getStimulusCurrent(int m) const
{
	return neurons[m].getStimulusCurrent();
}

/*** getFluctCurrent ***
 * Returns fluctuating current evoked by external synapses in neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the instantaneous fluctuating current in nA */
double getFluctCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getFluctCurrent();
}
double getFluctCurrent(int m) const
{
	return neurons[m].getFluctCurrent();
}

/*** getConstCurrent ***
 * Returns the constant current elicited by the surrounding network (not this network!) in neuron (i|j) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the constant current in nA */
double getConstCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getConstCurrent();
}
double getConstCurrent(int m) const
{
	return neurons[m].getConstCurrent();
}

/*** getSigma ***
 * Returns the standard deviation of the white noise entering the external current *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the standard deviation in nA s^(1/2) */
double getSigma(int i, int j) const
{
	return neurons[cNN(i,j)].getSigma();
}
double getSigma(int m) const
{
	return neurons[m].getSigma();
}

/*** getSynapticCurrent ***
 * Returns the synaptic current that arrived in the previous time step *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the synaptic current in nA */
double getSynapticCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getSynapticCurrent();
}
double getSynapticCurrent(int m) const
{
	return neurons[m].getSynapticCurrent();
}

#if DENDR_SPIKES == ON
/*** getDendriticCurrent ***
 * Returns the current that dendritic spiking caused in the previous time step *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the synaptic current in nA */
double getDendriticCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getDendriticCurrent();
}
double getDendriticCurrent(int m) const
{
	return neurons[m].getDendriticCurrent();
}
#endif

#if COND_BASED_SYN == ON
/*** getExcSynapticCurrent ***
 * Returns the internal excitatory synaptic current that arrived in the previous time step *
 * - return: the excitatory synaptic current in nA */
double getExcSynapticCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getExcSynapticCurrent();
}
double getExcSynapticCurrent(int m) const
{
	return neurons[m].getExcSynapticCurrent();
}

/*** getInhSynapticCurrent ***
 * Returns the internal inhibitory synaptic current that arrived in the previous time step *
 * - return: the inhibitory synaptic current in nA */
double getInhSynapticCurrent(int i, int j) const
{
	return neurons[cNN(i,j)].getInhSynapticCurrent();
}
double getInhSynapticCurrent(int m) const
{
	return neurons[m].getInhSynapticCurrent();
}
#endif

/*** getActivity ***
 * Returns true if neuron (i|j) is spiking in this instant of duration dt *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: whether neuron is firing or not */
bool getActivity(int i, int j) const
{
	return neurons[cNN(i,j)].getActivity();
}
bool getActivity(int m) const
{
	return neurons[m].getActivity();
}

/*** spikeAt ***
 * Returns whether or not a spike has occurred at a given spike, begins searching *
 * from latest spike *
 * - int t_step: the timebin at which the spike should have occurred *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: true if a spike occurred, false if not */
bool spikeAt(int t_step, int i, int j) const
{
	return neurons[cNN(i,j)].spikeAt(t_step);
}
bool spikeAt(int t_step, int m) const
{
	return neurons[m].spikeAt(t_step);
}

/*** getSpikeTime ***
 * Returns the spike time for a given spike number (in temporal order, starting with 1) of neuron (i|j) *
 * ATTENTION: argument n should not exceed the result of getSpikeHistorySize() *
 * - int n: the number of the considered spike (in temporal order, starting with 1) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the spike time for the n-th spike (or -1 if there exists none) */
int getSpikeTime(int n, int i, int j) const
{
	return neurons[cNN(i,j)].getSpikeTime(n);
}
int getSpikeTime(int n, int m) const
{
	return neurons[m].getSpikeTime(n);
}

/*** removeSpikes ***
 * Removes a specified set of spikes from history, to save memory *
 * - int start: the number of the spike to start with (in temporal order, starting with 1)
 * - int end: the number of the spike to end with (in temporal order, starting with 1) *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
void removeSpikes(int start, int end, int i, int j)
{
	neurons[cNN(i,j)].removeSpikes(start, end);
}
void removeSpikes(int start, int end, int m)
{
	neurons[m].removeSpikes(start, end);
}


/*** getSpikeCount ***
 * Returns the number of spikes that have occurred since the last reset (including those that have been removed) of neuron (i|j) *
 * - return: the number of spikes */
int getSpikeCount(int i, int j) const
{
	return neurons[cNN(i,j)].getSpikeCount();
}
int getSpikeCount(int m) const
{
	return neurons[m].getSpikeCount();
}

/*** getSpikeHistorySize ***
 * Returns the current size of the spike history vector of neuron (i|j) *
 * - return: the size of the spike history vector */
int getSpikeHistorySize(int i, int j) const
{
	return neurons[cNN(i,j)].getSpikeHistorySize();
}
int getSpikeHistorySize(int m) const
{
	return neurons[m].getSpikeHistorySize();
}

/*** setCurrentStimulus ***
 * Sets a current stimulus for neuron (i|j)
 * - Stimulus& _cst: shape of one stimulus period *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located */
void setCurrentStimulus(Stimulus& _cst, int i, int j)
{
	neurons[cNN(i,j)].setCurrentStimulus(_cst);
}
void setCurrentStimulus(Stimulus& _cst, int m)
{
	neurons[m].setCurrentStimulus(_cst);
}

/*** getNumberIncoming ***
 * Returns the number of either inhibitory or excitatory incoming connections to this neuron *
 * from other neurons in the network *
 * - int type: the type of incoming connections (inh./exc.)
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the number of incoming connections */
int getNumberIncoming(int type, int i, int j) const
{
	return neurons[cNN(i,j)].getNumberIncoming(type);
}
int getNumberIncoming(int type, int m) const
{
	return neurons[m].getNumberIncoming(type);
}

/*** getNumberOutgoing ***
 * Returns the number of connections outgoing from this neuron to other *
 * neurons of a specific type *
 * - int type: the type of postsynaptic neurons (inh./exc.)
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: the number of outgoing connections */
int getNumberOutgoing(int type, int i, int j) const
{
	return neurons[cNN(i,j)].getNumberOutgoing(type);
}
int getNumberOutgoing(int type, int m) const
{
	return neurons[m].getNumberOutgoing(type);
}

/*** getPProteinAmount ***
 * Returns the LTP-related protein amount in a neuron *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: momentary LTP-related protein amount */
double getPProteinAmount(int i, int j) const
{
	return neurons[cNN(i,j)].getPProteinAmount();
}
double getPProteinAmount(int m) const
{
	return neurons[m].getPProteinAmount();
}

/*** getCProteinAmount ***
 * Returns the common protein amount in a neuron *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: momentary LTP-related protein amount */
double getCProteinAmount(int i, int j) const
{
	return neurons[cNN(i,j)].getCProteinAmount();
}
double getCProteinAmount(int m) const
{
	return neurons[m].getCProteinAmount();
}

/*** getDProteinAmount ***
 * Returns the LTD-related protein amount in a neuron *
 * - int i: the row where the neuron is located *
 * - int j: the column where the neuron is located *
 * - return: momentary LTD-related protein amount */
double getDProteinAmount(int i, int j) const
{
	return neurons[cNN(i,j)].getDProteinAmount();
}
double getDProteinAmount(int m) const
{
	return neurons[m].getDProteinAmount();
}


/*** saveNeuronParams ***
 * Saves all the neuron parameters (including the channel parameters) to a given file; *
 * all neurons have the same parameters, so the first one is taken */
void saveNeuronParams(ofstream *f) const
{
	neurons[0].saveNeuronParams(f);
}

/* =============================================================================================================================== */

};
