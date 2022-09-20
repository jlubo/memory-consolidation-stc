/************************************************************************************************
 ***  Class that allows to define stimulation intervals with different shape. The shape is    ***
 ***     periodically repeated and can add to non-deterministic contributions of certain      ***
 ***    frequency and magnitude. For non-deterministic stimulation, Poissonian, Gaussian,     ***
 ***                   or Ornstein-Uhlenbeck-type noise can be used.                          ***
 ************************************************************************************************/

/*** Copyright 2017-2022 Jannik Luboeinski ***
 *** licensed under Apache-2.0 (http://www.apache.org/licenses/LICENSE-2.0) ***/

class Stimulus
{

#if (STIM_TYPE != POISSON_STIMULATION && STIM_TYPE != OU_STIMULATION && STIM_TYPE != GAUSS_STIMULATION && STIM_TYPE != DET_STIMULATION)
	#error "Unsupported option for STIM_TYPE."
#endif

	/*** Interval class ***
	 * An object of this class represents a temporal interval and specifies the stimulation that occurs in that interval */
	class Interval
	{
		private:

		int n; // number of timesteps for the shape
		double* shape; // deterministic stimulus magnitude values for all timesteps
#if STIM_TYPE != DET_STIMULATION
		minstd_rand0 rg; // default uniform generator for random numbers
#if STIM_TYPE == POISSON_STIMULATION
		double expc_spikes; // expected Poisson spikes per timestep (Poisson spike occurrence frequency times the duration of one timestep)
		double poisson_contrib; // magnitude of contribution of one Poisson spike
		int N_P; // number of Poisson neurons
		//double* prob_dist; // probability distribution for firing of 1,2,3...N_P Poisson neurons in one timestep
		uniform_real_distribution<double> u_dist; // uniform distribution
#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION
		normal_distribution<double> n_dist; // normal distribution to obtain Gaussian white noise, constructed in Neuron class constructor
		double mean_stim; // mean of the Gaussian white noise or OU process used for stimulation
		double sigma_stim; // standard deviation of the Gaussian white noise or OU process used for stimulation ("discrete" standard deviation, contains 1/sqrt(dt) already)
		double expOU; // exponential decay factor for one timestep of OU process
		double noisePrefactorOU; // pre-factor for one timestep for white noise in OU formula
		double stim_prev; // value of the OU process in previous timestep
#endif
#endif
		public:

		const int start; // timestep at which the stimulus begins
		const int end; // timestep at which the stimulus begins

		/*** getShapeLength ***
		 * Returns the length of one period of a deterministic Stimulus shape as number of time bins *
		 * - return: period length */
		int getShapeLength() const
		{
			return n;
		}

		/*** addRectPulse ***
		 * Adds a rectangle pulse of a given length and at a given time *
		 * to the deterministic shape *
		 * - magnitude: height of the rectangle pulse
		 * - offset: time in steps (within a period) at which rectangle pulse shall start 
		 * - len: time length of the rectangle pulse in steps */
		void addRectPulse(double magnitude, int offset, int len)
		{
			if (offset >= n || offset < 0)
				throw invalid_argument("Invalid offset value.");

			for (int i=0; i<(len > n-offset ? n-offset : len); i++)
			{
				shape[offset+i] += magnitude;
			}
		}

		/*** addSinePulse ***
		 * Adds a sine pulse (a whole period) of a given period and at a given time *
		 * to the deterministic shape *
		 * - magnitude: amplitude of sine function
		 * - offset: time in steps (within a period) at which sinus pulse shall start 
		 * - period: time length of one sinus period in steps */
		void addSinePulse(double magnitude, int offset, int period)
		{
			if (offset >= n || offset < 0)
				throw invalid_argument(string("Invalid offset value (must be between 0 and ") + to_string(n) + string(")."));

			for (int i=0; i<(period > n-offset ? n-offset : period); i++)
			{
				shape[offset+i] += magnitude * sin( (2.0*M_PI/period) * i );
			}
		}

		/*** multiplyBy ***
		 * Multiplies the deterministic shape and possibly the Poisson contribution or the Gaussian mean by a real number *
		 * - r: number to multiply with */
		void multiplyBy(double r)
		{
			for(int i=0; i<n; i++)
			{
				shape[i] *= r;
			}
#if STIM_TYPE == POISSON_STIMULATION
			poisson_contrib *= r;
#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION
			mean_stim *= r;
			sigma_stim *= r;
#endif		
		}
	

#if STIM_TYPE == POISSON_STIMULATION
		/*** setPoissonStimulation ***
		 * Sets the strength, number of neurons and frequency for spike generation through a homogeneous Poisson process *
		 * (Poisson dynamics can be switched off by suitable parameters) *
		 * - strength: the strength of contribution by one Poisson spike *
		 * - N: number of Poisson neurons used *
		 * - nd_freq: the frequency for Poisson stimulation in Hz *
		 * - dt: duration of one timestep in seconds */
		void setPoissonStimulation(double strength, int N, double nd_freq, double dt)
		{
			if (strength > EPSILON && N > 0 && nd_freq > EPSILON)
			{
				poisson_contrib = strength;				
				expc_spikes = nd_freq * dt;
				N_P = N;			

				// compute probability distribution for firing of 1,2,3...N_P Poisson neurons in one timestep
				/*freeProbDist();
				prob_dist = new double[N_P];

				double cum_rem = 0.;
				for (int i=N_P-1; i>=0; i--)
				{
					prob_dist[i] = pow(fdt, i+1) - cum_rem; // probability that exactly i+1 Poisson neurons fire
					cum_rem += prob_dist[i];
				}*/
			}
			else // switch off dynamics
			{
				N_P = 0;
			}
		}

		/*** getPoissonContribution ***
		 * Returns the magnitude of contribution of one Poisson spike *
		 * - return: contribution of one Poisson spike */
		double getPoissonContribution() const
		{
			return poisson_contrib;
		}

		/*** getSpikeExpectance ***
		 * Returns the expected value of spikes per timestep generated by the homogeneous Poisson process *
		 * - return: expected value of spikes per timestep */
		double getSpikeExpectance() const
		{
			return expc_spikes;
		}

		/*** getPoissonNeuronNumber ***
		 * Returns the number of Poisson neurons *
		 * - return: number of Poisson neurons */
		int getPoissonNeuronNumber() const
		{
			return N_P;
		}


#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION

		/*** getMean ***
		 * Returns the mean stimulus magnitude *
		 * - return: mean of Gaussian white noise or OU process */
		double getMean() const
		{
			return mean_stim;
		}

		/*** getSigma ***
		 * Returns the standard deviation of the stimulus *
		 * - return: standard deviation of Gaussian white noise or OU process */
		double getSigma() const
		{
			return sigma_stim;
		}

		
		/*** getExpOU ***
		 * Returns the exponential decay factor for OU formula *
		 * - return: exponential decay factor for OU formula */
		double getExpOU() const
		{
			return expOU;
		}

		/*** getNoisePrefactorOU ***
		 * Returns the pre-factor for white noise in the OU formula *
		 * - return: pre-factor for white noise in OU formula */
		double getNoisePrefactorOU() const
		{
			return noisePrefactorOU;
		}

		/*** setGaussStimulation ***
		 * Sets the mean and standard deviation for Gaussian white noise simulating inputs *
		 * from a population of neurons using the diffusion approximation *
		 * (dynamics can be switched off by suitable parameters) *
		 * - strength: the strength of stimulation through one hypothetical input synapse *
		 * - N: number of neurons simulated by the Gaussian process *
		 * - nd_freq: the frequency for Gauss stimulation in Hz *
		 * - dt: duration of one timestep in seconds */
		void setGaussStimulation(double strength, int N, double nd_freq, double dt)
		{
			if (strength > EPSILON && N > 0 && nd_freq > EPSILON)
			{
				mean_stim = strength * N * nd_freq; // mean
				sigma_stim = strength * sqrt(N * nd_freq / dt); // standard deviation
			}
			else // switch off dynamics
			{
				mean_stim = 0.;
				sigma_stim = 0.;
			}
		}

		/*** setOUStimulation ***
		 * Sets the mean and standard deviation for OU process simulating inputs *
		 * from a population of neurons *
		 * (dynamics can be switched off by suitable parameters) *
		 * - strength: the strength of stimulation through one hypothetical input synapse *
		 * - N: number of neurons simulated by the OU process *
		 * - nd_freq: the frequency for Gauss stimulation in Hz *
		 * - tau_OU: the correlation time of the Ornstein-Uhlenbeck process in s *
		 * - dt: duration of one timestep in seconds */
		void setOUStimulation(double strength, int N, double nd_freq, double tau_OU, double dt)
		{
			if (strength > EPSILON && N > 0 && nd_freq > EPSILON)
			{
				mean_stim = strength * N * nd_freq; // mean
				sigma_stim = strength * sqrt(N * nd_freq) / sqrt(2.*tau_OU); // standard deviation (here no division by sqrt(dt) because of OU formula)
				expOU = exp(-dt/tau_OU); // exponential decay factor for OU formula
				noisePrefactorOU = sigma_stim * sqrt(1. - exp(-2.*dt/tau_OU)); // pre-factor for white noise in OU formula
				stim_prev = mean_stim; // set start value on mean
			}
			else // switch off dynamics
			{
				mean_stim = 0.;
				sigma_stim = 0.;
				expOU = 0.;
				noisePrefactorOU = 0.;
				stim_prev = 0.;
			}
		}

#endif

		/*** setDetStimulation ***
		 * Sets the period for deterministic stimulation shape *
		 * _n: the number of timesteps for one shape period */
		void setDetStimulation(int _n)
		{
			freeShape(); // free if there was an old shape
				
			n = _n;

			if (n > 0) // if there is a new shape
			{
				shape = new double[n];
				for (int i = 0; i < n; i++)
				{
					shape[i] = 0.0;
				}
			}
			
			
		}
			
		/*** isShapeDefined ***
		 * Tests whether or not a deterministic shape of the Interval object is defined *
		 * - return: true if a shape is defined, false if not */
		inline bool isShapeDefined() const
		{
			if (n > 0)
				return true;
			else
				return false;
		}

		/*** readStimulusShape ***
		 * Returns the value of the stimulus shape (the stimulus period) at a certain time bin *
		 * - i: time bin in stimulus period
		 * - return: stimulus shape at given time */
		double readStimulusShape(int i) const
		{
			if (isShapeDefined())
				return shape[i % n];
			else
				return 0.;
		}

		/*** freeShape ***
		 * Frees the memory reserved for the shape array, in case there is reserved memory */
		void freeShape()
		{
			if (isShapeDefined())
			{
				delete[] shape;
				n = 0;
			}
		}

		/*** get ***
		 * Returns the current magnitude of the stimulation in the interval that may consist of a deterministic *
		 * contribution and of a non-deterministic contribution *
		 * - t_step: the current timestep *
		 * - return: the magnitude of stimulation */
		double get(int t_step)
		{
			double ret;

			ret = readStimulusShape(t_step-start); // deterministic contribution

#if STIM_TYPE == POISSON_STIMULATION

			for (int i=0; i<N_P; i++) // TODO: could be optimized using prob_dist
			{
				if (u_dist(rg) < expc_spikes) // Poisson spiking 
					ret += poisson_contrib; // contribution by homogeneous Poisson process
			}

#elif STIM_TYPE == GAUSS_STIMULATION

			ret += mean_stim + sigma_stim * n_dist(rg); // Gaussian white noise contribution

#elif STIM_TYPE == OU_STIMULATION

			stim_prev = (stim_prev - mean_stim) * expOU + noisePrefactorOU * n_dist(rg) + mean_stim;
			ret += stim_prev; // OU process contribution

#endif
			return ret;
		}

		/*** copy ***
		 * Copies attributes from another Interval object into this (used by copy constuctor and assignment operator) *
		 * - org: reference to an Interval object */
	  	void copy(const Interval& org)
		{
			// deterministic stimulation
			n = org.getShapeLength();
			if (n > 0) // if there is a shape
			{
				shape = new double[n];
				for (int i = 0; i < n; i++)
				{
					shape[i] = org.readStimulusShape(i);
				}
			}

			// non-deterministic stimulation
#if STIM_TYPE == POISSON_STIMULATION
			poisson_contrib = org.getPoissonContribution();
			expc_spikes = org.getSpikeExpectance();
			N_P = org.getPoissonNeuronNumber();
			u_dist = uniform_real_distribution<double>(0.0,1.0);
#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION
			mean_stim = org.getMean();
			sigma_stim = org.getSigma();
			expOU = org.getExpOU();
			noisePrefactorOU = org.getNoisePrefactorOU();
			stim_prev = org.getMean();
#endif
		}


		/*** Assignment operator ***
		 * Reads attributes from an Interval object and copies them into an existing Interval object; *
		 * the variables "start" and "end" have to be the same in the original object and the object to assign to
		 * - org: reference to an Interval object */
		Interval& operator=(const Interval& org)
		{
			if (this == &org) 
				return *this;

			if (start != org.start || end != org.end)
				throw length_error("The Interval object being assigned has differing \"start\" or \"end\" value.");
			
			freeShape();
			copy(org);
		
			return *this;
		}

		/*** Copy constructor ***
		 * Reads attributes from an Interval object and copies them into a new Interval object *
		 * - org: reference to an Interval object */
	  	Interval(const Interval& org) : start(org.start), end(org.end)
		{
			copy(org);
#if STIM_TYPE != DET_STIMULATION
			rg = minstd_rand0(getClockSeed());
#endif
		}

		/*** Default constructor ***
		 * _start: the timestep with which the interval starts *
		 * _end: the timestep with which the interval ends */
		Interval(int _start, int _end) : start(_start), end(_end)
		{			
			n = 0;
#if STIM_TYPE == POISSON_STIMULATION
			expc_spikes = 0.;
			poisson_contrib = 0.;
			N_P = 0;
			u_dist = uniform_real_distribution<double>(0.0,1.0);
#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION
			mean_stim = 0.;
			sigma_stim = 0.;
			expOU = 0.;
			noisePrefactorOU = 0.;
			stim_prev = 0.;
#endif
#if STIM_TYPE != DET_STIMULATION
			rg = minstd_rand0(getClockSeed());
#endif
		}

		/*** ~Interval ***
	 	 * Destructor */
		~Interval()
		{
			freeShape();
		}
	};

	vector<Interval> intervals; // vector of begin/end timesteps for stimulation intervals
#if STIM_PREPROC == ON
	vector<bool> stim_flags; // pre-processed vector specifying if timesteps exhibit stimulation or not (should be stored in bits, not bytes!)
#endif
	int stimulation_start; // timestep at which all stimulation begins
	int stimulation_end; // timestep at which all stimulation ends
	double dt; // duration of one timestep in s
	char* ppdata; // pre-processed stimulus data (indices for each timestep that indicate if there is stimulation and what kind of stimulation it is)

	public:
#if STIM_PREPROC == ON
	/*** preProcess ***
	 * Pre-processes the stimulus such that runtime is reduced (at the cost of memory consumption) */
	void preProcess() 
	{
		if (stimulation_end < stimulation_start) // if there is no stimulation
			return;

		cout << "Stimulus pre-processing will use " << dtos((stimulation_end-stimulation_start) / pow2(1024.) / 8, 2) << " MB." << endl;

		for (int t_step=stimulation_start; t_step<stimulation_end; t_step++) // loop over timesteps
		{
			stim_flags[t_step-stimulation_start] = false;

			for (int i=0; i<intervals.size(); i++) // loop over stimulus intervals
			{
				if (t_step >= intervals[i].start && t_step < intervals[i].end)
					stim_flags[t_step-stimulation_start] = true;
			}
		}
	}
#endif

	/*** get ***
	 * If there is a stimulus defined for the given time, returns the stimulus magnitude at this time *
	 * (stimulation does only take place in defined stimulus intervals, the start of an interval marks *
	 * the start of a stimulus shape period) *
	 * - t_step: timestep at which to evaluate stimulus
	 * - return: stimulus at given time */
	double get(int t_step)
	{	
		if (t_step >= stimulation_start && t_step < stimulation_end)
		{
#if STIM_PREPROC == ON
			if (stim_flags[t_step-stimulation_start] == 1)
			{
#endif
				for (int i=0; i<intervals.size(); i++) // loop over stimulus intervals
				{
					if (t_step >= intervals[i].start && t_step < intervals[i].end)
						return intervals[i].get(t_step); // get stimulus for this interval and timestep
				}
#if STIM_PREPROC == ON
			}
#endif
		}
		return 0.; // timestep does not lie within any interval / does not contain stimulation
	}


	/*** stimExists ***
	 * If there is a stimulus defined for the given time, returns true *
	 * - t_step: timestep at which to evaluate stimulus
	 * - return: does stimulus exist at given time */
	bool stimExists(int t_step)
	{	
		if (t_step >= stimulation_start && t_step < stimulation_end)
		{
#if STIM_PREPROC == ON
			if (stim_flags[t_step-stimulation_start] != 0) // read out flag
				return true;
#else
			for (int i=0; i<intervals.size(); i++) // loop over stimulus intervals
			{
				if (t_step >= intervals[i].start && t_step < intervals[i].end)
					return true;
			}
#endif
		}
		return false; // timestep does not lie within any interval / does not contain stimulation
	}



	/*** addStimulationInterval ***
	 * Adds stimulation for an interval, given in absolute time; stimulation does only take place *
	 * within intervals that are added using this method *
	 * - t_step_start: the beginning of the stimulation interval
	 * - t_step_end: the end of the stimulation interval *
	 * - return: interval index (1-based) if it was added successfully; 0 if not */
	int addStimulationInterval(int t_step_start, int t_step_end)
	{
		if (t_step_start < 0 || t_step_end < 0 || (t_step_end <= t_step_start))
			throw invalid_argument("Invalid interval boundaries.");

		for (int i=0; i<intervals.size(); i++) // check if this interval overlaps with a previously defined one
		{
			if ( (t_step_start >= intervals[i].start && t_step_start <= intervals[i].end) // start lies within another interval
			  || (t_step_end >= intervals[i].start && t_step_end <= intervals[i].end) // end lies within another interval
			  || (intervals[i].start >= t_step_start && intervals[i].end <= t_step_end) ) // another interval lies between start and end
			{
				throw invalid_argument("Overlapping stimulation intervals");
			}
		}
		
		intervals.push_back(Interval(t_step_start, t_step_end));

		if (stimulation_start > t_step_start)
			stimulation_start = t_step_start; // update start time of stimulation
		if (stimulation_end < t_step_end)
			stimulation_end = t_step_end; // update end time of stimulation

#if STIM_PREPROC == ON
		stim_flags.assign(stimulation_end-stimulation_start, true); // initially, assume that all timesteps within the interval contain stimulation
#endif

		return intervals.size(); // index
	}

	
	/*** testIntervalIndex ***
	 * Tests the validity of an interval index; if its is invalid, throws an exception *
	 * - index: an interval */
	void testIntervalIndex(int index)
	{
		if (index <= 0 && index > intervals.size())
			throw invalid_argument("Invalid interval index.");
	}
	
	/*** setDetStimulation ***
	 * Sets the period for deterministic stimulation shape of a stimulus interval *
	 * (or deactivates deterministic stimulation) *
	 * - shape_freq [overloaded]: the repetition frequency for the shape
	 * - n_steps [overloaded]: the period duration for the period of the shape (in timesteps)
	 * - index: index of the interval the operation is applied to */
	void setDetStimulation(double shape_freq, int index)
	{
		testIntervalIndex(index);
		intervals[index-1].setDetStimulation(round(1/(shape_freq*dt)));	
	}
	void setDetStimulation(int n_steps, int index)
	{
		testIntervalIndex(index);
		intervals[index-1].setDetStimulation(n_steps);	
	}

	/*** addRectPulse ***
	 * Adds a rectangle pulse of a given length and at a given time *
	 * to the deterministic shape of a stimulus interval *
	 * - magnitude: height of the rectangle pulse
	 * - offset: time in steps (within a period) at which rectangle pulse shall start 
	 * - len: time length of the rectangle pulse in steps
	 * - index: index of the interval the operation is applied to */
	void addRectPulse(double magnitude, int offset, int len, int index)
	{
		testIntervalIndex(index);

		intervals[index-1].addRectPulse(magnitude, offset, len);
	}

	/*** addSinePulse ***
	 * Adds a sine pulse (a whole period) of a given period and at a given time *
	 * to the deterministic shape of a stimulus interval *
	 * - magnitude: amplitude of sine function
	 * - offset: time in steps (within a period) at which sinus pulse shall start 
	 * - period: time length of one sinus period in steps
	 * - index: index of the interval the operation is applied to */
	void addSinePulse(double magnitude, int offset, int period, int index)
	{
		testIntervalIndex(index);

		intervals[index-1].addSinePulse(magnitude, offset, period);
	}

	/*** multiplyBy ***
	 * Multiplies the deterministic shape and possibly the Poisson contribution or the Gaussian mean by a real number *
	 * (either of one specified interval or of all intervals) *
	 * - r: number to multiply with *
	 * - index [optional]: index of the interval the operation is applied to */
	void multiplyBy(double r, int index)
	{
		testIntervalIndex(index);

		intervals[index-1].multiplyBy(r);
	}
	void multiplyBy(double r)
	{
		for (int i=0; i<intervals.size(); i++) // loop over stimulus intervals
		{
			intervals[i].multiplyBy(r);
		}
	}

#if STIM_TYPE == POISSON_STIMULATION

	/*** setPoissonStimulation ***
	 * Sets the strength, number of neurons and frequency for spike generation through a homogeneous Poisson process *
	 * (Poisson dynamics can be switched off by suitable parameters) *
	 * - strength: the strength of contribution by one Poisson spike *
	 * - N: number of Poisson neurons used *
	 * - nd_freq: the frequency for Poisson stimulation in Hz *
	 * - index: index of the interval the operation is applied to */
	void setPoissonStimulation(double strength, int N, double nd_freq, int index)
	{
		testIntervalIndex(index);
		
		intervals[index-1].setPoissonStimulation(strength, N, nd_freq, dt);
	}


#elif STIM_TYPE == GAUSS_STIMULATION

	/*** setGaussStimulation ***
	 * Sets the mean and standard deviation for Gaussian white noise simulating inputs *
	 * from a population of neurons using the diffusion approximation *
	 * (dynamics can be switched off by suitable parameters) *
	 * - strength: the strength of stimulation through one hypothetical input synapse *
	 * - N: number of neurons simulated by the Gaussian process *
	 * - nd_freq: the frequency for Gauss stimulation in Hz *
	 * - index: index of the interval the operation is applied to */
	void setGaussStimulation(double strength, int N, double nd_freq, int index)
	{
		testIntervalIndex(index);

		intervals[index-1].setGaussStimulation(strength, N, nd_freq, dt);
	}

#elif STIM_TYPE == OU_STIMULATION

	/*** setOUStimulation ***
	 * Sets the mean and standard deviation for an OU process simulating inputs *
	 * from a population of neurons using the diffusion approximation *
	 * (dynamics can be switched off by suitable parameters) *
	 * - strength: the strength of stimulation through one hypothetical input synapse *
	 * - N: number of neurons simulated by the Gaussian process *
	 * - nd_freq: the frequency for Gauss stimulation in Hz *
	 * - tau_OU: the correlation time of the Ornstein-Uhlenbeck process in s *
	 * - index: index of the interval the operation is applied to */
	void setOUStimulation(double strength, int N, double nd_freq, double tau_OU, int index)
	{
		testIntervalIndex(index);

		intervals[index-1].setOUStimulation(strength, N, nd_freq, tau_OU, dt);
	}

#endif

	/*** plotAll ***
	 * Plots all stimulation intervals using a GNUplot script *
	 * - name: the file name without suffix */
	void plotAll(string name) const
	{
		string txt, gpl, pdf;
		ofstream f;
      
		// prepare file names
		txt = name + string(".txt");
		gpl = name + string(".gpl");
		pdf = name + string(".pdf");

		// open data file
		f.open(txt);
		if (!f.is_open())
			throw runtime_error(string("File ") + string(txt) + string(" could not be opened."));

		cout << "Plotting Stimulus object labeled '" << name << "'" << endl;
		cout << "\tstimulation_start = " << stimulation_start*dt << " s, stimulation_end = " << stimulation_end*dt << " s" << endl;
		cout << "\tdt = " << dt << " s" << endl;

		// write data file
		f << fixed << (stimulation_start-100)*dt << "\t\t\t" // some timesteps before stimulation
		  << 0. << "\t\t\t" << false << "\t\t\t" << NAN << "\t\t\t" << NAN << "\t\t\t"
#if STIM_PREPROC == ON
		  << 0
#endif
		  << endl 
		  << fixed << (stimulation_start-1)*dt << "\t\t\t"
		  << 0. << "\t\t\t" << false << "\t\t\t" << NAN << "\t\t\t" << NAN << "\t\t\t"
#if STIM_PREPROC == ON
		  << 0
#endif
		  << endl;
		for (int t_step=stimulation_start; t_step<stimulation_end; t_step++) // loop over timesteps
		{
			double shape = 0.;
			bool poisson = false;
			double s_mean = NAN;
			double s_sigma = NAN;

			for (int i=0; i<intervals.size(); i++) // loop over stimulus intervals
			{
				if (t_step >= intervals[i].start && t_step < intervals[i].end)
				{
					if (intervals[i].isShapeDefined())
						shape = intervals[i].readStimulusShape(t_step-intervals[i].start);
#if STIM_TYPE == POISSON_STIMULATION
					poisson = (intervals[i].getPoissonNeuronNumber() > 0);
#elif STIM_TYPE == GAUSS_STIMULATION || STIM_TYPE == OU_STIMULATION
					if (intervals[i].getMean() > EPSILON || intervals[i].getSigma() > EPSILON)
					{
						s_mean = intervals[i].getMean();
						s_sigma = intervals[i].getSigma(); // sigma for white noise is scaled/distorted by 1/sqrt(dt)
					}
#endif
				}
			}

			f << fixed << t_step*dt << "\t\t\t" 
			  << shape << "\t\t\t" << poisson << "\t\t\t" << s_mean << "\t\t\t" << s_sigma << "\t\t\t"
#if STIM_PREPROC == ON
			  << stim_flags[t_step-stimulation_start]
#endif
			  << endl;
		}
		f << fixed << (stimulation_end+1)*dt << "\t\t\t" // some timesteps after stimulation
		  << 0. << "\t\t\t" << false << "\t\t\t" << NAN << "\t\t\t" << NAN << "\t\t\t"
#if STIM_PREPROC == ON
		  << 0
#endif
		  << endl
		  << fixed << (stimulation_end+100)*dt << "\t\t\t"
		  << 0. << "\t\t\t" << false << "\t\t\t" << NAN << "\t\t\t" << NAN << "\t\t\t"
#if STIM_PREPROC == ON
		  << 0
#endif
		  << endl;
		f.close();

		// open script file
		f.open(gpl);
		if (!f.is_open())
			throw runtime_error(string("File ") + string(gpl) + string(" could not be opened."));

		// write script file
		f << "set term pdf enhanced font \"Sans, 20\" color solid lw 2.5" << endl;
		f << "set output '" << pdf << "'" << endl;
		f << "set xlabel \"Time (s)\"" << endl;
		f << "set ylabel \"Stimulation (a.u.)\"" << endl;
		f << "set y2range [0:2]" << endl;
		f << "set key horizontal outside top left" << endl;
		f << "unset y2tics" << endl;
		f << "unset y2label" << endl;
		f << "set style fill transparent" << endl;

		int range = stimulation_end - stimulation_start;
		f << "plot [x=" << (stimulation_start-0.1*range-100)*dt << ":" << (stimulation_end+0.1*range+100)*dt << "] "
		  << "'" << txt << "' using 1:2 t 'Det. shape' with lines, "
#if STIM_TYPE == POISSON_STIMULATION
		  << "'" << txt << "' using 1:3 t 'Poisson' with filledcu x1, " 
#elif STIM_TYPE == GAUSS_STIMULATION
		  << "'" << txt << "' using 1:($4-$5):($4+$5) notitle lc rgb '#ccebff' with filledcurves, " 
		  << "'" << txt << "' using 1:4 t 'Gauss' lc rgb '#0099ff' with lines, " 
#elif STIM_TYPE == OU_STIMULATION
		  << "'" << txt << "' using 1:($4-$5):($4+$5) notitle lc rgb '#ccebff' with filledcurves, " 
		  << "'" << txt << "' using 1:4 t 'OU' lc rgb '#0099ff' with lines, " 
#endif
#if STIM_PREPROC == ON
		  << "'" << txt << "' using 1:6 t 'Flag' axes x1y2 with dots"
#endif
		  << endl;
		f.close();

		txt = string("gnuplot ") + gpl;
		system(txt.c_str());
	}
	
	/*** getStimulationIntervals ***
	 * Returns the vector of stimulation intervals *
	 * - return: reference to the vector of stimulation intervals */
	const vector<Interval>& getStimulationIntervals() const
	{
		return intervals;
	}

	/*** getDt ***
	 * Returns the timestep length *
	 * - return: the length of a timestep */
	const double getDt() const
	{
		return dt;
	}

#if STIM_PREPROC == ON
	/*** getStimFlags ***
	 * Returns the vector of pre-processed stimulation flags (specifiers for the existence of stimulation) *
	 * - return: reference to the vector of stimulation flags */
	const vector<bool>& getStimFlags() const
	{
		return stim_flags;
	}
#endif

	/*** getStimulationStart ***
	 * Returns the starting time of the first stimulation interval *
	 * - return: start time of the stimulation */
	const int getStimulationStart() const
	{
		return stimulation_start;
	}

	/*** getStimulationEnd ***
	 * Returns the end time of the last stimulation interval *
	 * - return: end time of the stimulation */
	const int getStimulationEnd() const
	{
		return stimulation_end;
	}

	/*** getStimulationDuration ***
	 * Returns the duration of the whole time span during which stimulation occurs *
	 * - return: duration of stimulation */
	const int getStimulationDuration() const
	{
		if (stimulation_end > stimulation_start)
			return stimulation_end - stimulation_start;
		else
			return 0;
	}
	
	/*** isSet ***
	 * Returns if any stimulation is set or not *
	 * - return: true if stimulation is set, false if not */
	const int isSet() const
	{
		if (stimulation_end > 0)
			return true;
		else
			return false;
	}

	/*** clear ***
	 * Clears all stimulus content */
	void clear()
	{
		vector<Interval>().swap(intervals);
		
		stimulation_start = numeric_limits<int>::max();
		stimulation_end = 0;
	}

	/*** Empty constructor ***
	 * Just for syntactic reasons */
  	Stimulus()
	{
		clear();
	}

	/*** Principal constructor ***
	 * Sets main characteristics of stimulus *
	 * - int _n: total timestep count of one period */
  	Stimulus(double _dt) : dt(_dt)
	{
		clear();
	}

	/*** Copy constructor ***
	 * Reads shape and properties from a Stimulus variable *
	 * and copies them into a new Stimulus (basically does *
	 * the same as assignment operator) *
	 * - Stimulus& st: reference to a variable of type Stimulus */
  	Stimulus(const Stimulus& _st) // : rg(getClockSeed()), u_dist((0.0,1.0))
	{
		clear();

		intervals = _st.getStimulationIntervals();
#if STIM_PREPROC == ON
		stim_flags = _st.getStimFlags();
#endif
		stimulation_start = _st.getStimulationStart();
		stimulation_end = _st.getStimulationEnd();
		dt = _st.getDt();
	}

	/*** Assignment operator ***
	 * Reads shape and properties from a Stimulus variable *
	 * and copies them into a new Stimulus (basically does *
	 * the same as copy constructor) *
	 * - Stimulus& st: reference to a variable of type Stimulus */
	Stimulus& operator=(const Stimulus& _st)
	{
		if (this != &_st)
		{
			clear();
			
			intervals = _st.getStimulationIntervals();
#if STIM_PREPROC == ON
			stim_flags = _st.getStimFlags();
#endif
			stimulation_start = _st.getStimulationStart();
			stimulation_end = _st.getStimulationEnd();
		}
		return *this;
	}


	/*** Destructor ***
	 * Frees reserved memory */
  	~Stimulus()
	{

	}

};

