#ifndef SIMUL_H
#define SIMUL_H

#include <vector>
#include <array>
#include <random>
#include "parameters.h"

struct State {
	// Positions of the particles
	std::vector< std::array<double, DIM> > positions;
	// Old positions in Y/Z
	std::vector< std::array<double, DIM-1> > oldPosYZ;
	// Forces between the particles
	std::vector< std::array<double, DIM> > forces;
	// Positions of the tracers after thermalization
	std::vector<double> initXTracers;
};

struct Observables {
	std::vector< std::vector<double> > moments;  // All the moments
};

int runSimulations(const Parameters &p);
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
		                   std::vector<Observables> &sumObs,
						   const unsigned int seed);
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
   					  std::mt19937 &rndGen);
void initState(State &state, const Parameters &p, std::mt19937 &rndGen);
void setInitXTracers(State &state, const Parameters &p);
void updateState(State &state, const Parameters &p, std::mt19937 &rndGen,
		         const bool thermalization=false);
void calcForcesBetweenParticles(State &state, const Parameters &p);
void keepInChannel(State &state, const Parameters &p);
void initObservables(std::vector<Observables> &obs, const Parameters &p);
void computeObservables(const State &state, const Parameters &p,
		                Observables &o);
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

void reflexionInCircle(const double xIn, const double yIn,
		               const double xOut, const double yOut,
					   const double R,
					   double &xFin, double &yFin);
void findIntersection(const double xIn, const double yIn,
		              const double xOut, const double yOut,
					  const double R,
					  double &xCross, double &yCross);
void basicReflexion(const double ux, const double uy,
		            double normalX, double normalY,
					double &xRefl, double &yRefl);
void solveSecondOrderEq(const double a, const double b, const double c,
                        double &sol1, double &sol2);
int sign(const double x);

#endif
