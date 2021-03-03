/***************************************************************
 *** Main function, instantiates the NetworkSimulation class ***
 ***************************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

#include "NetworkSimulation.cpp"

extern uint8_t data[] asm("_binary_code_zip_start"); // start of file attachment "code.zip"
extern uint8_t data_end[] asm("_binary_code_zip_end"); // end of file attachment "code.zip"
extern uint8_t data2[] asm("_binary_plotFunctions_py_start"); // start of file attachment "plotFunctions.py"
extern uint8_t data2_end[] asm("_binary_plotFunctions_py_end"); // end of file attachment "plotFunctions.py"

int Nl = 40; // number of neurons in one row of the excitatory population (total number: Nl^2)
int Nl_inh = 20; // number of neurons in one row of the inhibitory population (total number: Nl_inh^2)

double dt = 0.0002; // s, duration of one time step 
double t_max = 28810.0; // s, total duration of the simulation
double t_wfr = 0.5; // s, size of the sliding time window for firing rate computation

double w_ee = 1.0; // coupling strength for excitatory inputs to excitatory neurons as multiple of h_0 (thus equals h_0) 
double w_ei = 2.0; // coupling strength for excitatory inputs to inhibitory neurons as multiple of h_0
double w_ie = 4.0; // coupling strength for inhibitory inputs to excitatory neurons as multiple of h_0
double w_ii = 4.0; // coupling strength for inhibitory inputs to inhibitory neurons as multiple of h_0
double p_c = 0.1; // probability of connection between two neurons
double tau_syn = 0.005; // s, synaptic time constant

double I_0 = 0.15; // nA, mean current applied to the individual neurons (e.g., mean of the OU process)
double sigma_WN = 0.05; // nA s^1/2, standard deviation of the input current (e.g., std. dev. of the Gaussian white noise driving the OU process)

double theta_p = 3.0; // the calcium threshold for early-phase potentiation
double theta_d = 1.2; // the calcium threshold for early-phase depression
double Ca_pre = 0.6; // the postsynaptic calcium contribution evoked by a presynaptic spike (Li et al., 2016, not adjusted: 1.0)
double Ca_post = 0.16548; // the postsynaptic calcium contribution evoked by a postsynaptic spike (Li et al., 2016, not adjusted: 0.2758)
double sigma_plasticity = 9.1844 / sqrt(1000); // nA s, standard deviation of the Gaussian white noise contributing to early-phase plasticity

double theta_pro_P = 0.5; // protein synthesis threshold for pool P in units of h_0
double theta_pro_C = 0.5; // protein synthesis threshold for pool C in units of h_0
double theta_pro_D = 0.5; // protein synthesis threshold for pool D in units of h_0
double z_max = 1.; // the upper late-phase bound

double recall_fraction = 0.5; // the fraction of neurons in the cell assembly that is stimulated to trigger recall
int N_stim = 25; // number of hypothetical synapses used for stimulation

bool ff_enabled = true; // specifies if fast-forward mode (no computation of spiking dynamics) is allowed
int output_period = 10; // in general, every "output_period-th" time step, data will be recorded for plotting
int net_output_period = -1; //int(round(50.0/dt)); // every "net_output_period-th" time step, a network plot will be created;
                            // NOTE: should be chosen to be a multiple of the FF mode timestep delta_t;
                            // ATTENTION: should not be chosen too small because it causes removal of spikes


/*** main ***
 * - argc: the number of arguments *
 * - argv: array containing the arguments */
int main(int argc, char** argv) 
{
	char* cd = get_current_dir_name();
	string path = string(cd) + string("/") + dateStr("", true); // path to working directory
	//string path = string("./") + dateStr("", true); // path to working directory
	string parfile = ""; // a text file specifying the parameters
	string prot_learn = ""; // the used stimulus protocol for training
	string prot_recall = ""; // the used stimulus protocol for recall
	string purpose = ""; // short description of the purpose of this simulation
	bool purpose_set = false; // specifies whether the purpose variable has been set by argument (to enable setting no purpose)

	free(cd); // get_current_dir_name() uses malloc() for memory allocation

	// Process arguments
	for(int i=1; i<argc; i++)
	{
		char* pt;
		double read = 0.0;
		// specifying general settings
		if (strstr(argv[i], "-learn=") == argv[i])
		{
			pt = strstr(argv[i], "=") + 1;
			prot_learn = string(pt);
			continue;
		}
		else if (strstr(argv[i], "-recall=") == argv[i])
		{
			pt = strstr(argv[i], "=") + 1;
			prot_recall = string(pt);
			continue;
		}
		else if (strstr(argv[i], "-purpose=") == argv[i])
		{
			pt = strstr(argv[i], "=") + 1;
			purpose = string(pt);
			purpose_set = true;
			continue;
		}
		else if (strstr(argv[i], "-nopurpose") == argv[i])
		{
			purpose_set = true;
			continue;
		}
		else if (strstr(argv[i], "-noff") == argv[i])
		{
			ff_enabled = false;
			continue;
		}
		else if( (pt = strstr(argv[i], "=")) != NULL )
		{
			pt++;
			read = atof(pt);
		}
		
		// specifying numeric variables
		if (strstr(argv[i], "-p_c=") == argv[i] || strstr(argv[i], "-pc=") == argv[i])
			p_c = read;
		else if (strstr(argv[i], "-tau_syn=") == argv[i])
			tau_syn = read;
		else if (strstr(argv[i], "-w_ee=") == argv[i])
			w_ee = read;
		else if (strstr(argv[i], "-w_ei=") == argv[i])
			w_ei = read;
		else if (strstr(argv[i], "-w_ie=") == argv[i])
			w_ie = read;
		else if (strstr(argv[i], "-w_ii=") == argv[i])
			w_ii = read;
		else if (strstr(argv[i], "-z_max=") == argv[i] || strstr(argv[i], "-zmax=") == argv[i])
			z_max = read;
		else if (strstr(argv[i], "-theta_pro_P=") == argv[i])
			theta_pro_P = read;
		else if (strstr(argv[i], "-theta_pro_C=") == argv[i] || strstr(argv[i], "-theta_pro=") == argv[i])
			theta_pro_C = read;
		else if (strstr(argv[i], "-theta_pro_D=") == argv[i])
			theta_pro_D = read;
		else if (strstr(argv[i], "-I_0=") == argv[i] || strstr(argv[i], "-I_const=") == argv[i])
			I_0 = read;
		else if (strstr(argv[i], "-Nl_exc=") == argv[i] || strstr(argv[i], "-Nl=") == argv[i])
			Nl = int(read);
		else if (strstr(argv[i], "-Nl_inh=") == argv[i])
			Nl_inh = int(read);
		else if (strstr(argv[i], "-N_stim=") == argv[i])
			N_stim = int(read);
		else if (strstr(argv[i], "-t_max=") == argv[i] || strstr(argv[i], "-tmax=") == argv[i])
			t_max = read;
		else if (strstr(argv[i], "-t_wfr=") == argv[i])
			t_wfr = read;
		else if (strstr(argv[i], "-sigma_WN=") == argv[i] || strstr(argv[i], "-sigma_wn") == argv[i])
			sigma_WN = read;
		else if (strstr(argv[i], "-sigma_plasticity=") == argv[i])
			sigma_plasticity = read;
		else if (strstr(argv[i], "-output_period=") == argv[i])
			output_period = int(read);
		else if (strstr(argv[i], "-net_output_period=") == argv[i])
			net_output_period = int(read);
		else if (strstr(argv[i], "-theta_p=") == argv[i])
			theta_p = read;
		else if (strstr(argv[i], "-theta_d=") == argv[i])
			theta_d = read;
		else if (strstr(argv[i], "-Ca_pre=") == argv[i])
			Ca_pre = read;
		else if (strstr(argv[i], "-Ca_post=") == argv[i])
			Ca_post = read;
		else if (strstr(argv[i], "-r=") == argv[i])
			recall_fraction = read;

		// specifying parameter file (all parameters passed previously will be overwritten)
		//else if (strstr(argv[i], "-F") == argv[i] && (i+1) < argc && strstr(argv[i+1], "-") != argv[i])
		//{
		//	parfile = string(argv[++i]);

			// TODO read parfile
		//	continue;
		//}
		else
		{
			cout << "Unrecognized command line option: " << argv[i] << endl;
			return -1;
		}
	}

#ifdef TWO_NEURONS_ONE_SYNAPSE
	#ifndef PLASTICITY_OVER_FREQ
	cout << "TWO_NEURONS_ONE_SYNAPSE" << endl;
	t_max = 28800.0;
	output_period = 100;
	#else
	cout << "PLASTICITY_OVER_FREQ" << endl;
	#endif
	Nl = 2;
	Nl_inh = 0;
	I_0 = 0.;
	sigma_WN = 0.;
	N_stim = 1;
#elif defined SEEK_I_0
	purpose = "seek";
	purpose_set = true;
	output_period = 1000;
	net_output_period = int(t_max/dt) - 1;
#endif

	// Demand the purpose, if empty
	if (!purpose_set)
	{
		cout << "Purpose: ";
		getline(cin, purpose);
	}
	
	// Create working directory and files
	if (!prot_learn.empty())
		path += string("_") + prot_learn;
	if (!purpose.empty())
		path += string(" ") + purpose;
	system(concat("mkdir -p \"", path + "\"").c_str()); // create working directory
	if (chdir(path.c_str()) == -1) { // try to change directory
		showChDirErrMessage();
		return -1;
	}
	extractFileFromBinary("code.zip", data, data_end); // extracts code.zip out of the binary
	extractFileFromBinary("plotFunctions.py", data2, data2_end); // extracts plotFunctions.py out of the binary
	writeRunFile(string("run"), argc, argv); // writes the command line that was used to run this program in a file

	// Create NetworkSimulation object and set computational parameters
	NetworkSimulation sn = NetworkSimulation(Nl, Nl_inh, dt, t_max, p_c, sigma_plasticity, z_max, t_wfr, ff_enabled);

	// Seeking I_0
#ifdef SEEK_I_0 //----------------------------------------------------------------------

	ofstream seekICData(dateStr("_I_0_nu.txt"));
	ofstream seekICResult(dateStr("_I_0.txt"));
	double seekICVar; // firing rate value determined in current simulation
	double seekICVar_old = -1.0; // firing rate value determined in previous simulation, initially has to be set to -1.0
	
	double I_0_start = I_0; // starting with the adjusted I_0
	double I_0_step = 0.001; // current increment step size

	sn.setSeekICVar(&seekICVar); // pass address of seekICVar to NetworkSimulation object

#endif //-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-
	
	sn.setParams(I_0, sigma_WN, tau_syn, w_ee, w_ei, w_ie, w_ii,
	             prot_learn, prot_recall, output_period, net_output_period, N_stim, theta_p, theta_d, Ca_pre, Ca_post, theta_pro_P, theta_pro_C, theta_pro_D, recall_fraction);

#ifndef SEEK_I_0

	sn.simulate(path, true, purpose); // run the actual simulation

#else //------------------------------------------------------------------------------------

	while(true) // is cancelled by "break" once target firing rate has been found
	{
		sn.setParams(I_0, sigma_WN, tau_syn, w_ee, w_ei, w_ie, w_ii,
	             prot_learn, prot_recall, output_period, net_output_period, N_stim, theta_p, theta_d, Ca_pre, Ca_post, theta_pro_P, theta_pro_C, theta_pro_D, recall_fraction); // TODO: syntactic redundance

		sn.simulate(path, (I_0 == I_0_start) ? true : false, purpose); // run the actual simulation

		// after the simulation, seekICVar contains the mean firing rate
		seekICData << fixed << I_0 << "\t\t\t";
		if (seekICVar > 0) // zeros should not be used
			seekICData << seekICVar << endl;
		else
			seekICData << "nan" << endl;
	
		// if specified frequency has been crossed, compute I_0 using linear interpolation
		if ( (seekICVar_old > SEEK_I_0 && seekICVar <= SEEK_I_0) || // crossing from above SEEK_I_0
			  (seekICVar_old < SEEK_I_0 && seekICVar >= SEEK_I_0 && seekICVar_old >= 0.0) ) // crossing from below SEEK_I_0
		{
			double m = (seekICVar_old-seekICVar)/(I_0+I_0_step-I_0);
			double b = seekICVar_old - m*(I_0+I_0_step);
			seekICResult << fixed << (SEEK_I_0 - b) / m << endl;
			break;
		}

		// initial value of I_0 caused a firing rate lower than target firing rate -> increase I_0 from now on
		if (seekICVar < SEEK_I_0 && seekICVar_old < 0.0) // only if firing rate was below SEEK_I_0 and if this was the initial step
			I_0_step *= (-1.0);

		seekICVar_old = seekICVar;
		I_0 -= I_0_step;
	}
	seekICVar_old = -1.0; // reset seekICVar_old
	I_0_step = (I_0_step < 0 ? I_0_step * (-1.0) : I_0_step); // reset I_0_step
	//I_0 = I_0_start; // reset I_0 // do not reset - for next round should start with previous value
	seekICData << endl;

	seekICData.close();
	seekICResult.close();

#endif //-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-

	if (chdir(path.c_str()) == -1) { // try to re-change directory
		showChDirErrMessage();
		return -1;
	}

	//system("sudo shutdown -P now");

#ifdef TWO_NEURONS_ONE_SYNAPSE
	// save which was the kind of plasticity evoked by the given stimulus
	int plast_type = sn.getPlasticityType(); 
	double max_dev = sn.getMaxDev(); 
	
	chdir("..");
	ofstream fpt("plasticity_type.txt", ofstream::out | ofstream::app);

	fpt << prot_learn << "\t\t" << prot_recall << "\t\t" << max_dev << "\t\t" << plast_type << "\t\t# " << dateStr() << " ";

	if (plast_type == 0)
		fpt << "ELTP";
	else if (plast_type == 1)
		fpt << "ELTP-T";
	else if (plast_type == 2)
		fpt << "LLTP";
	else if (plast_type == 3)
		fpt << "ELTD";
	else if (plast_type == 4)
		fpt << "ELTD-T";
	else if (plast_type == 5)
		fpt << "LLTD";
	fpt << endl;
	
	fpt.close();
#endif

	return 0;
}
