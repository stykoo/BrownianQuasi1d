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
 * Contain the definition of the abstract class Simulation
 * and of generic functions related to the simulation.
 */

#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <array>
#include <random>
#include "parameters.h"

#define MAX_ITERS_INIT 100

struct Observables {
	std::vector<double> pos;  // Positions of the tracers (along x)
	std::vector<double> displ;  // Displacements of the tracers (along x)
};

class Simulation {
	public:
		Simulation(const Parameters &p);
		virtual ~Simulation() {}
		int run(std::vector<Observables> &obs, std::mt19937 &rndGen);
	
	protected:
		const Parameters p;
		const double noise;
		// Positions of the tracers after thermalization
		std::vector<double> initXTracers;
		// Distributions of random numbers
		std::uniform_real_distribution<double> distribUnif;
		std::normal_distribution<double> distribNormal;


		void setInitXTracers();
		int computeObservables(Observables &o);
		bool isOrdered();

		virtual int init(std::mt19937 &rndGen) = 0;
		virtual void update(std::mt19937 &rndGen,
				            const bool thermalization = false) = 0;
		virtual double getPosX(const long i) = 0;
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed, bool *failed);
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
   					  std::mt19937 &rndGen);
void initObservables(std::vector<Observables> &obs, const Parameters &p);
void addObservables(std::vector<Observables> &obs1,
		            const std::vector<Observables> &obs2, const Parameters &p);
int exportObservables(const std::vector<Observables> &sumObs,
		              const Parameters &p);

// Translate x into interval [-L/2, L/2[
inline double periodicBC(const double x, const double L) {
    return x - L * std::round(x / L);
}

// Compute a^b with b a positive integer
template<typename T, typename U>
T mypow(const T a, const U b) {
    if (b <= 0) {
        return 1;
	} else if (b % 2 == 1) {
        return a * mypow(a, b - 1);
	} else {
		T c = mypow(a, b / 2);
		return c * c;
	}
}

// Return the sign of a number (+1 if positive or null, -1 if negative)
inline int sign(const double x) {
    return (0. <= x) - (x < 0.);
}

#endif
