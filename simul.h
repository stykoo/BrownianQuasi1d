#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <array>
#include <random>
#include "parameters.h"

struct Observables {
	std::vector< std::vector<double> > moments;  // All the moments
};

class Simulation {
	public:
		Simulation(const Parameters &p);
		virtual ~Simulation() {}
		void run(std::vector<Observables> &obs, std::mt19937 &rndGen);
	
	protected:
		const Parameters p;

		virtual void init(std::mt19937 &rndGen) = 0;
		virtual void setInitXTracers() = 0;
		virtual void update(std::mt19937 &rndGen,
				            const bool thermalization = false) = 0;
		virtual void computeObservables(Observables &o) = 0;
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed);
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
