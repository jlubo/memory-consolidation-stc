/*****************************************************************************************
 ***     Pre-defined stimulus protocols for learning/consolidating cell assemblies     ***
 *****************************************************************************************/

/*** Copyright 2017-2021 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

/*** stimFunc ***
 * Function for either adding a fixed rectangular pulse, Gaussian stimulation, Ornstein-Uhlenbeck stimulation *
 * or Poissonian stimulation *
 * - st: pointer to the Stimulus object to be modified *
 * - frequency: the stimulation frequency *
 * - w_stim: coupling strength between input layer and receiving layer *
 * - N_stim: number of neurons in the input layer *
 * - tau_syn: the synaptic time constant *
 * - index: index of the interval the stimulation is added to */
void stimFunc(Stimulus* st, double frequency, double w_stim, int N_stim, double tau_syn, int index)
{
#if STIM_TYPE == DET_STIMULATION
	st->setDetStimulation(round(1/(shape_freq*dt)), index); // set period of deterministic stimulation
	st->addRectPulse(w_stim*N_stim, 0, 1, index); // add fixed rectangular pulse
#elif STIM_TYPE == GAUSS_STIMULATION
	st->setGaussStimulation(w_stim, N_stim, frequency, index); // use Gaussian white-noise stimulation
#elif STIM_TYPE == OU_STIMULATION
	st->setOUStimulation(w_stim, N_stim, frequency, tau_syn, index); // use stimulation by an Ornstein-Uhlenbeck process
#elif STIM_TYPE == POISSON_STIMULATION
	st->setPoissonStimulation(w_stim, N_stim, frequency, index); // use Poisson-distributed pulses
#endif

}


/*** createStimulusFromProtocols ***
 * Creates a Stimulus object according to specified stimulation protocols *
 * - prot_learn: string specifying the learning protocol that shall be used *
 * - prot_recall: string specifying the learning protocol that shall be used *
 * - dt: duration of one time step in seconds *
 * - w_stim: coupling strength between input layer and receiving layer *
 * - N_stim: number of neurons in the input layer *
 * - tau_syn: the synaptic time constant *
 * - return: Stimulus object */
Stimulus createStimulusFromProtocols(string prot_learn, string prot_recall, double dt, double w_stim, int N_stim, double tau_syn)
{
	Stimulus st(dt); // new Stimulus object
	double frequency; // stimulation frequency
	char* pt; // pointer for processing the protocol string
	int index; // interval index

	// learning protocol:
	///////////////////////

	pt = (char*) prot_learn.c_str();
	if (!prot_learn.compare("STET"))
	{
		frequency = 100.0; // Hz
		index = st.addStimulationInterval(int(round(3600.0/dt)), int(round(3601.0/dt))); // add stimulation interval
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(4200.0/dt)), int(round(4201.0/dt))); // add stimulation interval
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(4800.0/dt)), int(round(4801.0/dt))); // add stimulation interval
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("WTET"))
	{
		frequency = 100.0; // Hz
		int index = st.addStimulationInterval(int(round(3600.0/dt)), int(round(3600.2/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("SLFS"))
	{
		frequency = 20.0; // Hz
		for (int i=0;i<900;i++)
		{
			index = st.addStimulationInterval(int(round((3600.0 + i*23./20.)/dt)), int(round((3600.0 + i*23./20. + 3./20.)/dt))); // add start and end times of stimulation
			stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to interval
		}
	}
#if STIM_TYPE == POISSON_STIMULATION
	else if (!prot_learn.compare("WLFS"))
	{
		frequency = 1.0; // Hz
		index = st.addStimulationInterval(int(round(3600.0/dt)), int(round(4500.0/dt))); // add start and end time of stimulation
		st.setPoissonStimulation(w_stim, N_stim, frequency, index); // Poisson-distributed stimulation must be used
	}
#endif
	else if (strstr(pt, "F") == pt && strstr(pt, "D") > pt + 1) // "generic" protocol
	{
		double duration, at;
		char* pt2 = strstr(pt, "D");
		char* pt3 = strstr(pt, "at");

		pt2[0] = '\0'; // terminate first string for frequency (replacing 'D' by '\0')
		pt++; // skip 'F'
		frequency = atof(pt);

		if (pt3 > 0) // if time of occurrence is specified
			pt3[0] = '\0'; // terminate string for duration (replacing 'a' by '\0')
		pt2++; // skip '\0'
		duration = atof(pt2);

		if (pt3 > 0) // if time of occurrence is specified
		{
			pt3 += 2; // skip '\0' and 't'
			at = atof(pt3);

			index = st.addStimulationInterval(int(round(at/dt)), int(round((at + duration/10)/dt))); // add stimulation interval for actual recall
			stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		}
		else // "standard" learning at t = 100.0 s
		{
			index = st.addStimulationInterval(int(round(100.0/dt)), int(round((100.0 + duration/10)/dt))); // add stimulation interval (by specifying start and end time)
			stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		}
	}
	else if (!prot_learn.compare("TEST"))
	{
		frequency = 100.0; // Hz
		//index = st.addStimulationInterval(int(round(2.0/dt)), int(round(2.5/dt))); // add start and end time of stimulation
		index = st.addStimulationInterval(int(round(1.0/dt)), int(round(1.04/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(3.0/dt)), int(round(3.04/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(5.0/dt)), int(round(5.04/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("TRIPLET"))
	{
		frequency = 100.0; // Hz
		//index = st.addStimulationInterval(int(round(2.0/dt)), int(round(2.5/dt))); // add start and end time of stimulation
		index = st.addStimulationInterval(int(round(10.0/dt)), int(round(10.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(10.5/dt)), int(round(10.6/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(11.0/dt)), int(round(11.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		//index = st.addStimulationInterval(int(round(11.5/dt)), int(round(11.6/dt))); // add start and end time of stimulation
		//stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (strstr(pt, "TRIPLET") == pt && strstr(pt, "f") > pt) // TRIPLET with variable frequency
	{
		char* pt2 = strstr(pt, "f");
		pt2[0] = '\0'; // terminate first string for frequency (replacing 'f' by '\0')
		pt2++; // skip 'f'
		frequency = atof(pt2);

		index = st.addStimulationInterval(int(round(10.0/dt)), int(round(10.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(10.5/dt)), int(round(10.6/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(11.0/dt)), int(round(11.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval

	}
	else if (!prot_learn.compare("TESTA"))
	{
		frequency = 50.0; // Hz
		//index = st.addStimulationInterval(int(round(2.0/dt)), int(round(2.5/dt))); // add start and end time of stimulation
		index = st.addStimulationInterval(int(round(1.0/dt)), int(round(1.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(1.5/dt)), int(round(1.6/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(2.0/dt)), int(round(2.1/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		index = st.addStimulationInterval(int(round(2.5/dt)), int(round(2.6/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		//index = st.addStimulationInterval(int(round(3.0/dt)), int(round(3.1/dt))); // add start and end time of stimulation
		//stimFunc(&st, 2*frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("TEST10"))
	{
		frequency = 100.0; // Hz
		index = st.addStimulationInterval(int(round(1.0/dt)), int(round(2.0/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("TEST3"))
	{
		frequency = 100.0; // Hz
		index = st.addStimulationInterval(int(round(1.0/dt)), int(round(1.3/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_learn.compare("TESTSAV"))
	{
		frequency = 100.0; // Hz
		index = st.addStimulationInterval(int(round(0.1/dt)), int(round(0.2/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}

	// recall protocol:
	/////////////////////

	pt = (char*) prot_recall.c_str();
	if (strstr(pt, "F") == pt && strstr(pt, "D") > pt + 1) // "generic" protocol
	{
		double duration, at;
		//double strength = -1.; // TODO can be implemented as a parameter
		char* pt2 = strstr(pt, "D");
		char* pt3 = strstr(pt, "at");
		int index2, index3;

		pt2[0] = '\0'; // terminate first string for frequency (replacing 'D' by '\0')
		pt++; // skip 'F'
		frequency = atof(pt);

		if (pt3 > 0) // if time of occurrence is specified
			pt3[0] = '\0'; // terminate string for duration (replacing 'a' by '\0')
		pt2++; // skip '\0'
		duration = atof(pt2);

		if (pt3 > 0) // if time of occurrence is specified
		{
			pt3 += 2; // skip '\0' and 't'
			at = atof(pt3);

			index = st.addStimulationInterval(int(round(at/dt)), int(round((at + duration/10)/dt))); // add stimulation interval for actual recall
			stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		}
		else // "standard" recall at t = 150.0 s
		{
			index = st.addStimulationInterval(int(round(150.0/dt)), int(round((150.0 + duration/10)/dt))); // add stimulation interval for actual recall
			stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		}

		//if (strength >= 0.)
		//{
		//	st.multiplyBy(strength, index); // alter strength of stimulation
		//}
	}
	else if (!prot_recall.compare("TESTA"))
	{
		frequency = 100.0; // Hz
		//index = st.addStimulationInterval(int(round(1.0/dt)), int(round(1.2/dt))); // add start and end time of stimulation
		//stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		//index = st.addStimulationInterval(int(round(3.0/dt)), int(round(3.2/dt))); // add start and end time of stimulation
		index = st.addStimulationInterval(int(round(3.0/dt)), int(round(3.005/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}
	else if (!prot_recall.compare("TESTSAV"))
	{
		frequency = 100.0; // Hz
		//index = st.addStimulationInterval(int(round(1.0/dt)), int(round(1.2/dt))); // add start and end time of stimulation
		//stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
		//index = st.addStimulationInterval(int(round(3.0/dt)), int(round(3.2/dt))); // add start and end time of stimulation
		index = st.addStimulationInterval(int(round(0.2/dt)), int(round(0.205/dt))); // add start and end time of stimulation
		stimFunc(&st, frequency, w_stim, N_stim, tau_syn, index); // actually add stimulation to the interval
	}

#if STIM_PREPROC == ON
	st.preProcess();
#endif

	return st;
}

/*** createOscillStimulus ***
 * Creates a Stimulus object with sinusoidal oscillating input during the whole simulation *
 * - dt: duration of one time step in seconds *
 * - tb_max: number of time steps for the whole simulation *
 * - period: period for the sine-shaped input current *
 * - mean: mean value for the sine-shaped input current *
 * - amplitude: amplitude for the sine-shaped input current *
 * - return: new Stimulus object with  oscillation */
Stimulus createOscillStimulus(double dt, int tb_max, int period, double mean, double amplitude)
{
	Stimulus st(dt); // new Stimulus object
	int index; // interval index

	index = st.addStimulationInterval(0, tb_max); // add stimulation interval

	st.setDetStimulation(period, index); // set deterministic stimulation
	st.addRectPulse(mean, 0, period, index); // add fixed rectangular pulse (serves as mean for the sine oscillation)
	st.addSinePulse(amplitude, 0, period, index); // add sine pulse with defined period and amplitude

#if STIM_PREPROC == ON
	st.preProcess();
#endif

	return st;
}
