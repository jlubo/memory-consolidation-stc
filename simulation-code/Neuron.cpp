/*****************************************
 *** Model of a single LIF/MAT2 neuron ***
 *****************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***

 *** Uses the library boost/archive - ***
 *** Copyright 2002 Robert Ramey ***
 *** licensed under the Boost Software License v1.0 (http://www.boost.org/LICENSE_1_0.txt) ***

 *** Uses the library boost/serialization -  ***
 *** Copyright 2002 Robert Ramey, 2005 Matthias Troyer ***
 *** licensed under the Boost Software License v1.0 (http://www.boost.org/LICENSE_1_0.txt) ***/

#include <random>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/vector.hpp>
#include "Stimulus.cpp"

using namespace std;

/*** Neuron class ***
 * Represents one neuron */
class Neuron {	

friend class boost::serialization::access;

private:

/*** Computational parameters ***/
double dt; // s, one timestep for numerical simulation

/*** State variables ***/
double V; // mV, the current membrane potential at the soma
#if DENDR_SPIKES == ON
double I_dendr; // nA, the current evoked by dendritic spikes
double I_dendr_A; // nA, first contribution to the current evoked by dendritic spikes
double I_dendr_B; // nA, second contribution to the current evoked by dendritic spikes
double I_dendr_C; // nA, third contribution to the current evoked by dendritic spikes
double dendr_inp_integral; // nC, the integral over the charge deposited in the dendrite in the last 2 ms (or conductivity 
                           // in case of conductance-based synapses)
vector<double> dendr_inp_history; // nC, vector containing the PSC amplitudes of the last 2 ms
double dendr_int_window; // s, dendritic integration window
double refractory_dendr; // s, time span until refractory period for dendritic spikes is over
double t_ref_dendr; // s, absolute refractory period for dendritic spikes - has to be at least one timestep!
double dendr_spike_threshold; // nC, threshold that dendr_inp_integral has to cross for a dendritic spike to be evoked
double I_dendr_A_amp; // nC, amplitude of first current component of dendritic spikes
double I_dendr_B_amp; // nC, amplitude of second current component of dendritic spikes
double I_dendr_C_amp; // nC, amplitude of third current component of dendritic spikes
double tau_dendr_A; // s, decay time constant of first current component of dendritic spikes
double tau_dendr_B; // s, decay time constant of first current component of dendritic spikes
double tau_dendr_C; // s, decay time constant of first current component of dendritic spikes
#endif
#if NEURON_MODEL == MAT2
double ad_th; // mV, the adaptive voltage threshold
double exp1; // mV, the fast component of the adaptive voltage threshold
double exp2; // mV, the slow component of the adaptive voltage threshold
#endif
double p_P; // the protein amount for LTP in this neuron
double p_C; // the common protein amount for LTP and LTP in this neuron
double p_D; // the protein amount for LTD in this neuron
double I_stim; // nA, the externally applied stimulus current
double I_bg; // nA, the external background noise current (computed using Gaussian noise or an OU process with mean 0)
#if COND_BASED_SYN == ON
double I_int_exc; // nA, the synaptic input from excitatory network neurons affecting this neuron
double I_int_inh; // nA, the synaptic input from inhibitory network neurons affecting this neuron
double V_exc_syn_rev; // mV, the reversal potential for conductance-based, excitatory synapses
double V_inh_syn_rev; // mV, the reversal potential for conductance-based, inhibitory synapses
#endif
double I_int; // nA, the synaptic input from other network neurons affecting this neuron
double refractory; // s, time span until absolute refractory period is over
bool active; // specifies if neuron is currently spiking
int spike_count; // the total number of spikes that occurred since the last reset
vector<int> spike_history; // vector of all spike times (in units of timesteps) in the process since the last reset
int spike_history_reserve; // the maximum number of spikes
int inh_incoming; // number of incoming inhibitory connections in a network
int exc_incoming; // number of incoming excitatory connections in a network
int inh_outgoing; // number of outgoing connections to inhibitory neurons in a network
int exc_outgoing; // number of outgoing connections to excitatory neurons in a network
vector<int> outgoing; // vector of all outgoing connections to neurons in a network
Stimulus cst; // current stimulus for this neuron
minstd_rand0 rg; // default uniform generator for random numbers (seed is chosen in constructor)
normal_distribution<double> norm_dist; // normal distribution to obtain Gaussian white noise, constructed in Neuron class constructor
#if STIM_TYPE == POISSON_STIMULATION
bool poisson_neuron; // indicates if the neuron solely serves as Poisson spike generator
#endif

protected:

/*** Physical parameters ***/
double tau_mem; // s, the membrane time constant
double R_mem; // MΩ, resistance of the cell membrane
double V_rev; // mV, the reversal potential of the neuron
double V_reset; // mV, the reset potential of the neuron
double t_ref; // s, absolute refractory period - has to be at least one timestep!
#if NEURON_MODEL == LIF
double V_th; // mV, the threshold potential of the neuron
double V_spike; // mV, the height of an action potential
#elif NEURON_MODEL == MAT2
double ad_th_limit; // mV, the limit of the adaptive voltage threshold for time toward infinity
#endif
double tau_OU; // s, correlation time of the Ornstein-Uhlenbeck process
double sigma_WN; // nA s^1/2, standard deviation of the Gaussian white noise driving the OU process
#if SYNAPSE_MODEL == MONOEXP
double sigma_OU; // nA, standard deviation of the OU process
#endif
double I_0; // nA, the mean of the external synaptic inputs (i.e., of the OU process)
int type; // the type of this neuron (inhibitory/excitatory - TYPE_INH/TYPE_EXC)

/*** normalRandomNumber ***
 * Returns a random number drawn from a normal distribution with standard deviation 1 and mean 0 *
 * - return: the random number of type double (technically in units of sqrt(s)) */
inline double normalRandomNumber()
{
	double nd = norm_dist(rg);

	/* Check for mean and standard deviation *
	static int i = 0;
	static double mean =0.0;
	static double var =0.0;
	const int av_steps = 1000000;
	i++;
	mean += nd;
	var += pow2(nd); // compute variance based on an assumed mean of 0.0
	if (i == av_steps)
	{
		cout << endl << "mean: " << mean / double(av_steps) << endl;
		cout << "stddev: " << sqrt(var / double(av_steps)) << endl;
		mean = 0.0;
		var = 0.0;
		i=0;
	} */
	return nd;
}

public:

/*** saveNeuronParams ***
 * Saves all the neuron parameters to a given file */
void saveNeuronParams(ofstream *f) const
{
	*f << endl;
	*f << "Neuron parameters:" << endl;
	*f << "tau_mem = " << tau_mem << " s" << endl;
	*f << "R_mem = " << R_mem << " MΩ" << endl;
	*f << "V_rev = " << V_rev << " mV" << endl;
	*f << "t_ref = " << t_ref << " s" << endl;
#if NEURON_MODEL == LIF
	*f << "V_reset = " << V_reset << " mV" << endl;
	*f << "V_th = " << V_th << " mV" << endl;
	*f << "V_spike = " << V_spike << " mV" << endl;
#endif
#if COND_BASED_SYN == ON
	*f << "V_exc_syn_rev = " << V_exc_syn_rev << " mV" << endl;
	*f << "V_inh_syn_rev = " << V_inh_syn_rev << " mV" << endl;
#endif
	*f << "tau_OU = " << tau_OU << " s" << endl;
	*f << "sigma_WN = " << sigma_WN << " nA s^1/2" << endl;
	*f << "I_0 = " << I_0 << " nA" << endl;
#if COND_BASED_SYN == ON
	
#endif
#if DENDR_SPIKES == ON
	*f << "t_ref_dendr = " << t_ref_dendr << " s" << endl;
	*f << "dendr_spike_threshold = " << dendr_spike_threshold << " nC" << endl;
	*f << "dendr_int_window = " << dendr_int_window << " s" << endl;
	*f << "I_dendr_A_amp = " << I_dendr_A_amp << " nA" << endl;
	*f << "I_dendr_B_amp = " << I_dendr_B_amp << " nA" << endl;
	*f << "I_dendr_C_amp = " << I_dendr_C_amp << " nA" << endl;
	*f << "tau_dendr_A = " << tau_dendr_A << " s" << endl;
	*f << "tau_dendr_B = " << tau_dendr_B << " s" << endl;
	*f << "tau_dendr_C = " << tau_dendr_C << " s" << endl;
#endif
}

/*** serialize ***
 * Saves all state variables to a file using serialization from boost *
 * - ar: the archive stream *
 * - version: the archive version */
template<class Archive> void serialize(Archive &ar, const unsigned int version)
{
	ar & V;
#if DENDR_SPIKES == ON
	ar & I_dendr;
	ar & I_dendr_A;
	ar & I_dendr_B;
	ar & I_dendr_C;
	ar & dendr_inp_integral;
	ar & dendr_inp_history;
	ar & refractory_dendr;
#endif
#if NEURON_MODEL == MAT2
	ar & ad_th;
	ar & exp1;
	ar & exp2;
#endif
	ar & p_P;
	ar & p_C;
	ar & p_D;
	ar & I_stim;
	ar & I_bg;
#if COND_BASED_SYN == ON
	ar & I_int_exc;
	ar & I_int_inh;
	ar & V_exc_syn_rev;
	ar & V_inh_syn_rev;
#endif
	ar & I_int;
	ar & refractory;
	ar & active;
	ar & spike_count;
	ar & spike_history;
}

/*** getNumberIncoming ***
 * Returns the number of either inhibitory or excitatory incoming connections to this neuron *
 * from other neurons in a network *
 * - int type: the type of incoming connections (inh./exc.)
 * - return: the number of incoming connections */
int getNumberIncoming(int type) const
{
	if (type == TYPE_INH)
		return inh_incoming;
	else if (type == TYPE_EXC)
		return exc_incoming;
	else
		throw invalid_argument("Invalid neuron type.");
}

/*** getNumberOutgoing ***
 * Returns the number of connections outgoing from this neuron to other *
 * neurons (of a specific type) *
 * - int type [optional]: the type of the postsynaptic neuron (inh./exc.)
 * - return: the number of outgoing connections */
int getNumberOutgoing(int type) const
{
	if (type == TYPE_INH)
		return inh_outgoing;
	else if (type == TYPE_EXC)
		return exc_outgoing;
	else
		throw invalid_argument("Invalid neuron type.");
}
int getNumberOutgoing() const
{
	return outgoing.size();
}

/*** incNumberIncoming ***
 * Increases the number of incoming connections to this neuron (only to be used while *
 * a network is being built) *
 * - int type: the type of the presynaptic neuron (inh./exc.) */
void incNumberIncoming(int type)
{
	if (type == TYPE_INH)
		inh_incoming++;
	else if (type == TYPE_EXC)
		exc_incoming++;
}

/*** addOutgoingConnection ***
 * Adds an outgoing connection from this neuron to another one (only to be used while a network is being built) *
 * - index: the index of the neuron that receives the connection *
 * - int type: the type of the receiving neuron (inh./exc.) */
void addOutgoingConnection(int index, int type)
{
	outgoing.push_back(index);

	if (type == TYPE_INH)
		inh_outgoing++;
	else if (type == TYPE_EXC)
		exc_outgoing++;
}

/*** getOutgoingConnection ***
 * Returns the index of a neuron receiving an outgoing connection *
 * - arr_index: the index of the connection in the array *
 * - return: the index of the neuron that receives the connection */
int getOutgoingConnection(int arr_index) const
{
	return outgoing[arr_index];
}

/*** getVoltage ***
 * Returns the membrane potential of the neuron *
 * - return: the membrane potential in mV */
double getVoltage() const
{
	return V;
}

/*** getMembraneResistance ***
 * Returns the membrane resistance of the neuron *
 * - return: the membrane resistance in MΩ */
double getMembraneResistance() const
{
	return R_mem;
}

/*** getVoltageThreshold ***
 * Returns the value of the (possibly) dynamic membrane threshold of the neuron *
 * - return: the membrane threshold in mV */
double getVoltageThreshold() const
{
#if NEURON_MODEL == LIF
	return V_th;
#elif NEURON_MODEL == MAT2
	return ad_th;
#endif
}

/*** getCurrent ***
 * Returns total external current affecting the neuron *
 * - return: the instantaneous current in nA */
double getCurrent() const
{
	return I_stim+I_bg+I_int;
}

/*** getNetCurrent ***
 * Returns total current affecting the neuron, including external and leak currents *
 * - return: the instantaneous current in nA */
double getNetCurrent() const
{
	return I_stim+I_bg+I_int+(V_rev-V)/R_mem;
}

/*** getStimulusCurrent ***
 * Returns current evoked by external stimulation *
 * - return: the instantaneous current stimulus in nA */
double getStimulusCurrent() const
{
	return I_stim;
}

/*** getBGCurrent ***
 * Returns current external background current accounting for external inputs *
 * - return: the instantaneous external background current in nA */
double getBGCurrent() const
{
	return I_bg;
}

/*** getConstCurrent ***
 * Returns the mean of the external current *
 * - return: the constant current in nA */
double getConstCurrent() const
{
	return I_0;
}

/*** getSigma ***
 * Returns the standard deviation of the white noise entering the external current *
 * - return: the standard deviation in nA s^(1/2) */
double getSigma() const
{
	return sigma_WN;
}

/*** getSynapticCurrent ***
 * Returns the internal synaptic current that arrived in the previous timestep *
 * - return: the synaptic current in nA */
double getSynapticCurrent() const
{
	return I_int;
}

/*** setSynapticCurrent ***
 * Sets the synaptic input current from other neurons within the network to this neuron *
 * - _I_int: the synaptic current in nA */
void setSynapticCurrent(const double _I_int)
{
	I_int = _I_int;
}

/*** setSynapticPotential ***
 * Given the postsynaptic potential, sets the synaptic input current from other neurons *
 * within the network to this neuron *
 * - _V_PSP: the amplitude of the postsynaptic potential in mV */
void setSynapticPotential(const double _V_PSP)
{
	I_int = _V_PSP / R_mem;
}

#if COND_BASED_SYN == ON
/*** getExcSynapticCurrent ***
 * Returns the internal excitatory synaptic conductance of the previous timestep *
 * - return: the excitatory synaptic conductance in nS */
double getExcSynapticCurrent() const
{
	return I_int_exc;
}

/*** getInhSynapticCurrent ***
 * Returns the internal inhibitory synaptic conductance of the previous timestep *
 * - return: the inhibitory synaptic conductance in nS */
double getInhSynapticCurrent() const
{
	return I_int_inh;
}

/*** setExcSynapticCurrent ***
 * Sets the synaptic input current from excitatory neurons within the network to this neuron *
 * - _I_int_exc: the synaptic current/conductance in nA/nS */
void setExcSynapticCurrent(const double _I_int_exc)
{
	I_int_exc = _I_int_exc;
}

/*** setExcSynapticPotential ***
 * Given the postsynaptic potential, sets the synaptic input current from excitatory neurons *
 * within the network to this neuron *
 * - _V_PSP: the amplitude of the exc. postsynaptic potential in mV */
void setExcSynapticPotential(const double _V_PSP)
{
	I_int_exc = _V_PSP / R_mem;
}

/*** setInhSynapticCurrent ***
 * Sets the synaptic input current from inhibitory neurons within the network to this neuron *
 * - _I_int_inh: the synaptic current/conductance in nA/nS */
void setInhSynapticCurrent(const double _I_int_inh)
{
	I_int_inh = _I_int_inh;
}

/*** setInhSynapticPotential ***
 * Given the postsynaptic potential, sets the synaptic input current from inhibitory neurons *
 * within the network to this neuron *
 * - _V_PSP: the amplitude of the inh. postsynaptic potential in mV */
void setInhSynapticPotential(const double _I_int_inh)
{
	I_int_inh = _V_PSP / R_mem;
}
#endif

/*** increaseExcSynapticCurrent ***
 * Adds a specified contribution to the excitatory synaptic input current/conductance *
 * - contr: the contribution to be added to the synaptic current/conductance in nA/nS */
void increaseExcSynapticCurrent(const double contr)
{
#if COND_BASED_SYN == ON
	I_int_exc += contr;
#else
	I_int += contr;
#endif
}

/*** increaseExcSynapticPotential ***
 * Given the postsynaptic potential, adds a contribution to the excitatory synaptic input current/conductance *
 * - _V_PSP: the amplitude of the exc. postsynaptic potential in mV */
void increaseExcSynapticPotential(const double _V_PSP)
{
	double contr = _V_PSP / R_mem;
#if COND_BASED_SYN == ON
	I_int_exc += contr;
#else
	I_int += contr;
#endif
}

/*** increaseInhSynapticCurrent ***
 * Adds a specified contribution to the inhibitory synaptic input current/conductance *
 * - contr: the contribution to be added to the synaptic current/conductance in nA/nS */
void increaseInhSynapticCurrent(const double contr)
{
#if COND_BASED_SYN == ON
	I_int_inh += contr;
#else
	I_int -= contr;
#endif
}

/*** increaseInhSynapticPotential ***
 * Given the postsynaptic potential, adds a contribution to the inhibitory synaptic input current/conductance *
 * - _V_PSP: the amplitude of the inh. postsynaptic potential in mV */
void increaseInhSynapticPotential(const double _V_PSP)
{
	double contr = _V_PSP / R_mem;
#if COND_BASED_SYN == ON
	I_int_inh += contr;
#else
	I_int -= contr;
#endif
}

#if DENDR_SPIKES == ON
/*** updateDendriteInput ***
 * Adds the the amplitude of a PSC current to the history and to the integral of the deposited charge *
 * - psc_amplitude: the PSC amplitde in nA/nS/nC */
void updateDendriteInput(const double psc_amplitude)
{
	dendr_inp_history[dendr_inp_history.size()-1] += psc_amplitude;
	dendr_inp_integral += psc_amplitude;
}

/*** getDendriticCurrent ***
 * Returns the current that dendritic spiking caused in the previous timestep *
 * - return: the synaptic current in nA */
double getDendriticCurrent() const
{
	return I_dendr;
}

/*** compDendriticCurrent ***
 * Computes the dendritic current at time t as in Jahnke et al., 2015 *
 * - t: the time counted from the beginning of a dendritic spike *
void compDendriticCurrent(double t)
{    
    I_dendr = -55.*exp(-t/0.0002) + 64.*exp(-t/0.0003) - 9.*exp(-t/0.0007);
}*/
#endif

/*** getActivity ***
 * Returns true if the neuron is spiking in this instant of duration dt *
 * ATTENTION: can cause complications when used in a population for it does *
 * not include information about the exact spike time *
 * - return: whether neuron is firing or not */
bool getActivity() const
{
	if (active)
		return true;
	else
		return false;
}

/*** getSpikeTime ***
 * Returns the spike time for a given spike number (ATTENTION: argument n should not exceed the result of getSpikeHistorySize() -
 * this is not checked here for performance reasons) *
 * - int n: the number of the spike (in temporal order, starting with 1)
 * - return: the spike time in units of time bins for the n-th spike (or -1 if it does not exist) */
int getSpikeTime(int n) const
{
	return spike_history[n-1];
}

/*** spikeAt ***
 * Returns whether or not a spike has occurred at a given timestep, begins searching *
 * from latest spike *
 * - int t_step: the time bin at which the spike should have occurred
 * - return: true if a spike occurred, false if not */
bool spikeAt(int t_step) const
{
	for(int i=spike_history.size()-1; i>=0; i--)
	{
		if (spike_history[i] == t_step) // value that is looked for
			return true;
		else if (spike_history[i] < t_step) // already below value that is looked for
			return false;
	}
	return false;
}

/*** spikesInInterval ***
 * Counts the number of spikes in a given interval *
 * - int tb_start: the time bin at which the interval begins *
 * - int tb_end: one time bin after the interval ends *
 * - return: the number of spikes in the given interval */
int spikesInInterval(int tb_start, int tb_end) const
{
	int count = 0;

	for(int i=spike_history.size()-1; i>=0; i--)
	{
		if (spike_history[i] >= tb_start && spike_history[i] < tb_end)
			count++;
	}

	return count;
}

/*** hideSpikeTime ***
 * "Hides" the spike time for a given spike number, which means that *
 * it will be stored as a negative number, such that it will not anymore *
 * be taken into consideration for synaptic transmission *
 * (WAS ONLY NECESSARY FOR AN APPROXIMATION THAT IS NOT USED ANYMORE) *
 * - int n: the number of the spike (in temporal order, starting with 1)
 * - return: the spike time for the n-th spike (or -1.0 if it does not exist) */
void hideSpikeTime(int n)
{
	if (n <= spike_history.size() && n >= 1)
		spike_history[n-1] *= -1;
}

/*** removeSpikes ***
 * Removes a specified set of spikes from history, to save memory *
 * - int start: the number of the spike to start with (in temporal order, starting with 1)
 * - int end: the number of the spike to end with (in temporal order, starting with 1) */
void removeSpikes(int start, int end)
{
	spike_history.erase(spike_history.begin()+start-1,spike_history.begin()+end);
}

/*** getSpikeCount ***
 * Returns the number of spikes that have occurred since the last reset (including those that have been removed) *
 * - return: the number of spikes */
int getSpikeCount() const
{
	return spike_count;
}

/*** getSpikeHistorySize ***
 * Returns the current size of the spike history vector *
 * - return: the size of the spike history vector */
int getSpikeHistorySize() const
{
	return spike_history.size();
}

/*** processTimeStep ***
 * Processes one timestep (of duration delta_t) for the neuron * 
 * - int tb_step: timestep at which to evaluate stimulus (< 0 before stimulus onset) *
 * - int tb_init: initial timestep for simple decay process (should be positive only in decaying state!) */
void processTimeStep(int tb_step, int tb_init)
{
	double delta_t; // duration of the timestep in seconds, either dt or tb_step-tb_init

	if (tb_init < 0)
		delta_t = dt;
	else
		delta_t = (tb_step - tb_init) * dt;

#if SYNAPSE_MODEL == DELTA
	I_bg = normalRandomNumber() * sqrt(1/delta_t) * sigma_WN + I_0;
#elif SYNAPSE_MODEL == MONOEXP
	I_bg = (I_bg-I_0) * exp(-delta_t/tau_OU) + normalRandomNumber() * sqrt(1. - exp(-2.*delta_t/tau_OU)) * sigma_OU + I_0; // compute external synaptic input in nA
#endif

#if COND_BASED_SYN == ON
	//I_int = (I_int_exc * V_exc_syn_rev + I_int_inh * V_inh_syn_rev - (I_int_exc + I_int_inh) * V) / 1000.; // divide by 1000 to get from nS*mV=pA to nA
	I_int = (I_int_exc * (V_exc_syn_rev - V) + I_int_inh * (V_inh_syn_rev - V)) / 1000.; // divide by 1000 to get from nS*mV=pA to nA
#endif

#if STIM_TYPE == POISSON_STIMULATION
	if (poisson_neuron == false)
	{
#endif

#if NEURON_MODEL == MAT2

		// MAT(2) neuron
		V = V * exp(-delta_t/tau_mem) + R_mem*(I_int + I_bg) * (1. - exp(-delta_t/tau_mem)); // compute mem. pot. in mV (analytical solution)

		exp1 = exp1 * exp(-delta_t/0.01); // fast threshold relaxation
		exp2 = exp2 * exp(-delta_t/0.2); // slow threshold relaxation

		if (active)
		{

			exp1 = exp1 + 0.015; // add new spike with full contribution alpha_1
			exp2 = exp2 + 0.003; // add new spike with full contribution alpha_2
			
			active = false;
		}
		
		ad_th = ad_th_limit + exp1 + exp2; // update adaptive threshold

#elif NEURON_MODEL == LIF

#if DENDR_SPIKES == ON

		// exponential decay of dendritic spikes
		I_dendr_A *= exp(- delta_t / tau_dendr_A);
		I_dendr_B *= exp(- delta_t / tau_dendr_B);
		I_dendr_C *= exp(- delta_t / tau_dendr_C);

		if (refractory_dendr > EPSILON) // if in refractory period for dendritic spikes
		{
			refractory_dendr -= delta_t;
		}
		else
		{
			if (dendr_inp_integral > dendr_spike_threshold) // threshold has been crossed
			{
				// dendrite spike contributions do not have to be added up because possible remaining contributions should have decayed to zero
				I_dendr_A = -55.; //I_dendr_A -= 55.;
				I_dendr_B = 64.; //I_dendr_B += 64.;
				I_dendr_C = -9.; //I_dendr_C -= 9.;

				refractory_dendr = t_ref_dendr;
			}
		}

		//compDendriticCurrent(tb_step*delta_t - I_dendr_A);
		I_dendr = I_dendr_A + I_dendr_B + I_dendr_C;
			
		dendr_inp_integral -= dendr_inp_history[0]; // remove oldest contributions from integral
		dendr_inp_history.erase(dendr_inp_history.begin()); // remove oldest contributions from history
		dendr_inp_history.push_back(0.); // add slot for new contributions

#endif

		// LIF
		//V += delta_t/tau_mem * (- V + V_rev + R_mem*(I_bg + I_int)); // compute mem. pot. in mV (Euler method)
		V = V * exp(-delta_t/tau_mem) + (V_rev + R_mem*(  I_bg 
			                                        + I_int
#if DENDR_SPIKES == ON
		                                                + I_dendr
#endif
		                                                         )) * (1. - exp(-delta_t/tau_mem)); // compute mem. pot. in mV (analytical solution)
#endif

#if STIM_TYPE == POISSON_STIMULATION
	} // poisson_neuron == false

	else
	{
		active = false;
		V = V_reset;
	} // poisson_neuron == true
#endif
	if (cst.isSet() && tb_init < 0 && abs(I_stim = cst.get(tb_step)) > EPSILON) // stimulation; get stimulus current in nA
	{
#if STIM_TYPE == POISSON_STIMULATION || defined TWO_NEURONS_ONE_SYNAPSE_BASIC_EARLY
	#if NEURON_MODEL == MAT2
		V = ad_th + EPSILON; // definite spiking (the magnitude of I_stim is not important as long as it is greater than zero)
	#elif NEURON_MODEL == LIF
		V = V_th + EPSILON;
	#endif
#else
		V += R_mem * I_stim * (1. - exp(-delta_t/tau_mem));
#endif
	}

	if (refractory > EPSILON // if in refractory period
#if STIM_TYPE == POISSON_STIMULATION
	   && poisson_neuron == false
#endif
	   )
	{
#if NEURON_MODEL == LIF
		active = false;
		
		V = V_reset;
#endif
		refractory -= delta_t;
	}
	else
	{
#if NEURON_MODEL == LIF
		if (V >= V_th) // threshold crossing
		{
			V = V_reset;
#elif NEURON_MODEL == MAT2
		if (V >= ad_th) // threshold crossing
		{
#endif
			refractory = t_ref;
			if (tb_step >= 0) // count only those spikes that occur after the stimulus onset
				spike_history.push_back(tb_step);

			active = true;
			spike_count++;
		}
	}

}

/*** setProteinAmounts ***
 * Sets the protein amounts in the neuron *
 * - double _p_P: momentary protein amount for potentiation *
 * - double _p_C: momentary common protein amount *
 * - double _p_D: momentary protein amount for depression */
void setProteinAmounts(const double _p_P, const double _p_C, const double _p_D)
{
	p_P = _p_P;
	p_C = _p_C;
	p_D = _p_D;
}

/*** getPProteinAmount ***
 * Returns the protein amount for potentiation in the neuron *
 * - return: momentary protein amount */
double getPProteinAmount() const
{
	return p_P;
}

/*** getCProteinAmount ***
 * Returns the common protein amount in the neuron *
 * - return: momentary protein amount */
double getCProteinAmount() const
{
	return p_C;
}

/*** getDProteinAmount ***
 * Returns the protein amount for depression in the neuron *
 * - return: momentary protein amount */
double getDProteinAmount() const
{
	return p_D;
}

/*** setCurrentStimulus ***
 * Sets a current stimulus for the neuron *
 * - Stimulus& _cst: shape of one stimulus period */
void setCurrentStimulus(const Stimulus& _cst)
{
	cst = _cst;
}

/*** isStimulusSet ***
 * Returns true if a stimulus has been set *
 * - return: the value of isSet() in Stimulus class */
bool isStimulusSet() const
{
	return cst.isSet();
}

/*** multiplyCurrentStimulus ***
 * Multiplies the set current stimulus by a real number *
 * - double r: number to multiply */
void multiplyCurrentStimulus(double r)
{
	cst.multiplyBy(r);
}

/*** setConstCurrent ***
 * Sets the constant current (mean of Gaussian white noise or OU process) to a newly defined value *
 * - double _I_0: constant current in nA */
void setConstCurrent(double _I_0)
{
	I_0 = _I_0;
}

/*** setSigma ***
 * Sets the standard deviation of the external input current (i.e., Gaussian white noise or the OU process)*
 * - double _sigma: standard deviation in nA s^1/2 */
void setSigma(double _sigma)
{
	sigma_WN = _sigma;

#if SYNAPSE_MODEL == MONOEXP
	sigma_OU = sigma_WN / sqrt(2.*tau_OU); // sigma of the OU process; required unit for tau_OU is [s]
#endif
}

/*** setTauOU ***
 * Sets the time constant of the Ornstein-Uhlenbeck process *
 * (synaptic time constant of assumed input synapses) *
 * - double _tau_OU: the synaptic time constant in s*/
void setTauOU(double _tau_OU)
{
	tau_OU = _tau_OU;
}

/*** setType ***
 * Sets the type of this neuron (inhibitory/excitatory) *
 * - int _type: the neuron type */
void setType(int _type)
{
	type = _type;
}

/*** setSpikeHistoryMemory ***
 * Sets the RAM size that shall be reserved for the spike history *
 * (WARNING: choosing this too large can cause the simulation to be killed by the OS!) *
 * - int storage_steps: the size of the storage timespan in timesteps *
 * - return: the reserved size of the spike history vector */
int setSpikeHistoryMemory(int storage_steps)
{
	spike_history_reserve = int(round(storage_steps*(dt/t_ref)));
	return spike_history_reserve;
}

/*** getType ***
 * Returns the type of this neuron (inhbitory/excitatory) *
 * - return: the neuron type */
int getType() const
{
	return type;
}

#if STIM_TYPE == POISSON_STIMULATION
/*** setPoisson ***
 * Allows for making this neuron a Poisson neuron (used to generate Poissonian spikes only, without own voltage dynamics) *
 * - bool _poisson_neuron: indicates if the neuron shall be Poissonian */
void setPoisson(bool _poisson_neuron)
{
	poisson_neuron = _poisson_neuron;
}
#endif

/*** resetConnections ***
 * Resets the connections *
 * - only_exc: if true, only counters for excitatory connections are reset */
void resetConnections(bool only_exc = false)
{
	exc_outgoing = 0;
	exc_incoming = 0;

	if (!only_exc)
	{
		inh_outgoing = 0;
		inh_incoming = 0;
	}

	vector<int>().swap(outgoing);
}

/*** reset ***
 * Resets neuron to initial state */
void reset()
{
	vector<int>().swap(spike_history); // additionally to clearing the vector, reallocates it
	spike_history.reserve(spike_history_reserve); // pre-allocates space for the maximum number of spikes

	V = V_rev;

#if DENDR_SPIKES == ON
	I_dendr = 0.;
	I_dendr_A = 0.;//-1000.;
	I_dendr_B = 0.;
	I_dendr_C = 0.;
	dendr_inp_history.assign(int(round(dendr_int_window/dt)), 0.);
	dendr_inp_integral = 0.;
	refractory_dendr = 0.;
#endif
#if NEURON_MODEL == MAT2
	ad_th = ad_th_limit;
	exp1 = 0.0;
	exp2 = 0.0;
#endif
	refractory = 0.0; // neuron ready to fire
	p_P = 0.0;
	p_C = 0.0;
	p_D = 0.0;
	I_stim = 0.0;
	I_bg = 0.0;
#if COND_BASED_SYN == ON
	I_int_exc = 0.0;
	I_int_inh = 0.0;
#endif
	I_int = 0.0;

	active = false;
	spike_count = 0;
	
	rg.seed(getClockSeed()); // set new seed by clock's epoch
	norm_dist.reset(); // reset the normal distribution for random numbers
}


/*** Constructor ***
 * Sets all parameters on experimentally determined values *
 * _dt: size of a timestep in seconds */
Neuron(const double _dt) : 
	dt(_dt), rg(getClockSeed()), norm_dist(0.0,1.0)
{
	tau_mem = 0.010; // from Yamauchi et al., 2011

#if NEURON_MODEL == MAT2

	t_ref = 0.002; // from Kobayashi et al., 2009

#if defined TWO_NEURONS_ONE_SYNAPSE_LI2016
	R_mem = 0.0001; // from Li et al., 2016 (mind that the quantity R here is a different one than there!)
#else
	R_mem = 0.01;
	//R_mem = 50.0; // from Kobayashi et al., 2009 / Yamauchi et al., 2009
#endif

	V_rev = 0.0;
	V_reset = 0.0; // estimated to model afterhypolarization in combination with t_ref (see also Dayan & Abbott, p. 4)
	ad_th_limit = 0.005; // from Li et al., 2016

#elif NEURON_MODEL == LIF

	t_ref = 0.002; // estimated to model afterhypolarization in combination with V_reset
	R_mem = 10.0; // from Dayan & Abbott, fig. 5.
	V_th = -55.0; // from Dayan & Abbott, p. 162
	V_spike = 35.0; // estimated from sketch in Bachelor's thesis
	V_rev = -65.0; // from Dayan & Abbott, fig. 5.
	V_reset = -70.0; // estimated to model afterhypolarization in combination with t_ref (see also Dayan & Abbott, p. 4)

#endif
	sigma_WN = 0.;
#if SYNAPSE_MODEL == MONOEXP
	sigma_OU = 0.;
#endif
	I_0 = 0.;
	tau_OU = 0.005; // estimated

#if STIM_TYPE == POISSON_STIMULATION
	poisson_neuron = false;
#endif

#if COND_BASED_SYN == ON
	V_exc_syn_rev = 0.; // Moritz' Master's thesis
	V_inh_syn_rev = -70.; // Moritz' Master's thesis
#endif
#if DENDR_SPIKES == ON
	dendr_spike_threshold = 6.5; // estimated (Jahnke et al., 2015: 8.65)
	dendr_int_window = 0.002;
	t_ref_dendr = 0.005;
	I_dendr_A_amp = -55.;
	I_dendr_B_amp = 64.; 
	I_dendr_C_amp = -9.;
	tau_dendr_A = 0.0002;
	tau_dendr_B = 0.0003;
	tau_dendr_C = 0.0007;
#endif

	spike_history_reserve = 100;
	reset();
	resetConnections();
}

/*** Destructor ***
 * Frees the allocated memory */
~Neuron()
{

}


};


















