// Contains the generic functions for simulations
// Namespaces are used to choose the kind of simulation

#include <algorithm>
#include <fstream>
#include <map>
#include <chrono>
#include <thread>
#include <iomanip>
#include "simul.h"
#include "SimulPipe.h"
#include "SimulTonks.h"

// Load the parameters and initialize the distributions
Simulation::Simulation(const Parameters &params) : p(params),
	noise(std::sqrt(2. * p.temperature * p.timestep)),
	distribUnif(0., 1.), distribNormal(0., 1.) {
}

// Run the simulation. Return 1 if an incident happened.
int Simulation::run(std::vector<Observables> &obs, std::mt19937 &rndGen) {
	// Initialization of the positions
	init(rndGen);

	// Thermalization
	for (long j=0 ; j<p.nbItersTh ; ++j) {
		update(rndGen, true);
		if (p.checkOrder && !isOrdered()) {
			return 1;
		}
	}

	// Initial observables
	setInitXTracers();
	initObservables(obs, p);
	computeObservables(obs[0]);

	// Main loop
	for (long j=0 ; j<p.nbIters-1 ; ++j) {
		update(rndGen, false);
		computeObservables(obs[j+1]);
		if (p.checkOrder && !isOrdered()) {
			return 1;
		}
	}
	return 0;
}


// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;
	std::vector< std::vector<Observables> > allSumsObs(p.nbThreads);

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
			nbSimulsPerThread, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	std::vector<Observables> sumObs;
	initObservables(sumObs, p);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		addObservables(sumObs, allSumsObs[k], p);
	}

	int status = exportObservables(sumObs, p);

	return status;
}

// Run several simulations and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   std::vector<Observables> &sumObs,
						   const unsigned int seed) {
    // Random generator
	std::mt19937 rndGen(seed);

	// Initialize sumObs
	initObservables(sumObs, p);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation (of the right type!)
		Simulation *simul;
		if (p.simulName == "pipe") {
			simul = new SimulPipe(p);
		} else if (p.simulName == "tonks") {
			simul = new SimulTonks(p);
		} else {
			return;
		}
		std::vector<Observables> obs;
		int status = simul->run(obs, rndGen);
		delete simul;

		if (status != 0) {
			std::cerr << "Warning: Wrong order of the particles" << std::endl;
		}

		// Add the observables to the sum
		addObservables(sumObs, obs, p);
	}
}

// Initialize a vector of observables
void initObservables(std::vector<Observables> &obs, const Parameters &p) {
	obs.resize(p.nbIters);
	for (long t = 0 ; t < p.nbIters ; ++t) {
		obs[t].moments.resize(p.nbTracers);
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			obs[t].moments[i].assign(p.nbTracers - i, 0);
		}
	}
}

// Add observables o2 to observables o1.
void addObservables(std::vector<Observables> &obs1,
					const std::vector<Observables> &obs2, const Parameters &p)
{
	for (long t = 0 ; t < p.nbIters ; ++t) {
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				obs1[t].moments[i][j] += obs2[t].moments[i][j];
			}
		}
	}
}

// Export the observables to a file.
int exportObservables(const std::vector<Observables> &sumObs,
					  const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			file << " ";
			for (long k = j ; k < j + i + 1 ; ++k) {
				file << "x" << k+1;
			}
		}
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	// Data (we write the average and not the sum)
	for (long k = 0 ; k < p.nbIters ; ++k) {
		file << k * p.timestep;
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				file << " " << sumObs[k].moments[i][j] / p.nbSimuls;
			}
		}
		file << "\n";
	}

	file.close();

	return 0;
}
