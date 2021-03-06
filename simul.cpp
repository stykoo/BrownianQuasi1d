/*
Copyright Université Pierre et Marie Curie (2017)
Contributors: Alexis Poncet

alexis.poncet@ens.fr

This software is a computer program whose purpose is to simulate the
motion of interacting Brownian particles in a quasi-1d environment,
some of them being driven by an external force.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/*
 * BrownianQuasi1d
 * simul.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Contain the implementation of the class Simulation and of the
 * generic functions related to the simulation.
 */

#include <algorithm>
#include <fstream>
#include <map>
#include <chrono>
#include <thread>
#include <iomanip>
#include "simul.h"
#include "Simuls1d.h"
#include "SimulCanal.h"
#include "SimulPipe.h"

// Load the parameters and initialize the distributions
Simulation::Simulation(const Parameters &params) : p(params),
	noise(std::sqrt(2. * p.temperature * p.timestep)),
	distribUnif(0., 1.), distribNormal(0., 1.) {
		initXTracers.assign(p.nbTracers, 0.);
}

// Run the simulation.
// Return 1 if particles cross.
// Return 2 if initialization failed.
// Return 3 if a NaN appears in the observables
int Simulation::run(std::vector<Observables> &obs, std::mt19937 &rndGen) {
	// Initialization of the positions
	int status = init(rndGen);
	if (status) {
		return 2;
	}

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
		if ((j + 1) % p.skip == 0) {
			int status = computeObservables(obs[(j+1)/p.skip]);
			if (status != 0) {
				return 3;
			}

			if (p.checkOrder && !isOrdered()) {
				return 1;
			}
		}
	}
	return 0;
}

// Set the initial positions of the tracers
void Simulation::setInitXTracers() {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		initXTracers[i] = getPosX(p.idTracers[i]);
	}
}

// Compute the observables.
// Return 1 if a NaN appears, 0 otherwise.
int Simulation::computeObservables(Observables &o) {
	if (p.computeProfs) {
		double dx = periodicBC(getPosX(p.idTracers[0]) - initXTracers[0],
				               p.length);
		o.moments1[0] = dx;
		getPosRel(o.profiles[0], p.nbPtsProfs);
		for (long i = 1 ; i < DEFAULT_NB_MOMENTS ; ++i) {
			o.moments1[i] = dx * o.moments1[i-1];
			for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
				o.profiles[i][j] = dx * o.profiles[i-1][j];
			}
		}
	} else {
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			o.pos[i] = getPosX(p.idTracers[i]);
			o.displ[i] = periodicBC(getPosX(p.idTracers[i]) - initXTracers[i],
									p.length);
			if (std::isnan(o.pos[i]) || std::isnan(o.displ[i])) {
				return 1;
			}
		}
	}
	return 0;
}

// Check if the positions of the particles are ordered.
bool Simulation::isOrdered() {
	long c = 0;
	for (long i = 0 ; i < p.nbParticles-1 ; ++i) {
		if (getPosX(i) > getPosX(i+1)) {
			++c;
		}
	}
	if (getPosX(p.nbParticles-1) > getPosX(0)) {
		++c;
	}
	return (c < 2);
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
	bool failed = false;
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
			nbSimulsPerThread, std::ref(allSumsObs[i]), rd(), &failed)); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	if (failed) {
		return 2;
	}

	// Initialize the total sum
	std::vector<Observables> sumObs;
	initObservables(sumObs, p);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		addObservables(sumObs, allSumsObs[k], p);
	}

	int status;
	if (p.computeProfs) {
		status = exportMoments1(sumObs, p);
		if (status)
			return status;
		status = exportProfiles(sumObs, p);
	} else {
		status = exportObservables(sumObs, p);
	}
	return status;
}

// Run several simulations and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   std::vector<Observables> &sumObs,
						   const unsigned int seed, bool *failed) {
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
		} else if (p.simulName == "canal") {
			simul = new SimulCanal(p);
		} else if (p.simulName == "coulomb") {
			simul = new SimulCoulomb(p);
		} else if (p.simulName == "dipole") {
			simul = new SimulDipole(p);
		} else if (p.simulName == "coulCircle") {
			simul = new SimulCoulombCircle(p);
		} else if (p.simulName == "dipCircle") {
			simul = new SimulDipoleCircle(p);
		} else {
			return;
		}
		std::vector<Observables> obs;
		int status = simul->run(obs, rndGen);

		while (status == 1 || status == 3) {
			if (status == 1) {
				std::cerr << "Warning: Wrong order of the particles ";
			} else {
				std::cerr << "Warning: A NaN appeared in the observables ";
			}
			std::cerr << "(thread: " << std::this_thread::get_id()
				<< ", simulation: " << s+1 << "). "
				<< "Running simulation again." << std::endl;
			status = simul->run(obs, rndGen);
		}
		delete simul;

		if (status == 2) {
			std::cerr << "Could not initialize the system properly "
				<< "(thread: " << std::this_thread::get_id()
				<< ", simulation: " << s+1 << "). "
				<< "You should reduce the timestep." << std::endl;
			*failed = true;
			return;
		}

		// Add the observables to the sum
		addObservables(sumObs, obs, p);
	}
}

// Initialize a vector of observables
void initObservables(std::vector<Observables> &obs, const Parameters &p) {
	const long n = p.nbIters / p.skip;
	obs.resize(n);
	for (long t = 0 ; t < n ; ++t) {
		if (p.computeProfs) {
			obs[t].moments1.assign(DEFAULT_NB_MOMENTS, 0.);
			obs[t].profiles.assign(DEFAULT_NB_MOMENTS,
					               std::vector<double>(p.nbPtsProfs, 0.));
		} else {
			obs[t].pos.assign(p.nbTracers, 0.);
			obs[t].displ.assign(p.nbTracers, 0.);
		}
	}
}

// Add observables o2 to observables o1.
void addObservables(std::vector<Observables> &obs1,
					const std::vector<Observables> &obs2, const Parameters &p)
{
	const long n = p.nbIters / p.skip;
	for (long t = 0 ; t < n ; ++t) {
		if (p.computeProfs) {
			for (long i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
				obs1[t].moments1[i] += obs2[t].moments1[i];
				for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
					obs1[t].profiles[i][j] += obs2[t].profiles[i][j];
				}
			}
		} else {
			for (long i = 0 ; i < p.nbTracers ; ++i) {
				obs1[t].pos[i] += obs2[t].pos[i];
				obs1[t].displ[i] += obs2[t].displ[i];
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
		file << " ";
		file << "x" << i+1;
	}
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		file << " ";
		file << "dx" << i+1;
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	const long n = p.nbIters / p.skip;
	// Data (we write the average and not the sum)
	for (long k = 0 ; k < n ; ++k) {
		file << k * p.skip * p.timestep;
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			file << " " << sumObs[k].pos[i] / p.nbSimuls;
		}
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			file << " " << sumObs[k].displ[i] / p.nbSimuls;
		}
		file << "\n";
	}

	file.close();

	return 0;
}

// Export the moments of TP 1
int exportMoments1(const std::vector<Observables> &sumObs,
				   const Parameters &p) {
	std::string fname = p.output + "_moments.dat";

	std::ofstream file(fname);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (size_t i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
		file << " X^" << i+1;
	}
	file << "\n";

	const long n = p.nbIters / p.skip;

	for (long i = 0 ; i < n ; ++i) {
		double t = i * p.skip * p.timestep;
		file << t;
		for (size_t k = 0 ; k < DEFAULT_NB_MOMENTS  ; ++k) {
			file << " " << sumObs[i].moments1[k] / p.nbSimuls;
		}
		file << "\n";
	}
	file.close();
	return 0;
}

// Export the profiles to files.
int exportProfiles(const std::vector<Observables> &sumObs,
				   const Parameters &p) {
	const long n = p.nbIters / p.skip;

	for (long i = 0 ; i < n ; ++i) {
		double t = i * p.skip * p.timestep;
		std::string fname = p.output + "_prof" + std::to_string(t) + ".dat";

		std::ofstream file(fname);
		if (!file.is_open()) {
			return 1;
		}

		// Header
		file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
			<< "): ";
		printParameters(p, file);
		file << "\n# x";
		for (size_t i = 0 ; i < DEFAULT_NB_MOMENTS ; ++i) {
			file << " rho(x)*X^" << i;
		}
		file << "\n";
		for (long j = 0 ; j < p.nbPtsProfs ; ++j) {
			file << j * p.precProfs;
			for (size_t k = 0 ; k < DEFAULT_NB_MOMENTS  ; ++k) {
				file << " "
					 << (sumObs[i].profiles[k][j] * p.nbPtsProfs
						 / p.length / p.nbSimuls);
			}
		  	file << "\n";
		}
		file.close();
	}
	return 0;
}
