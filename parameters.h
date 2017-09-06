#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_NB_MOMENTS 2
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_THREADS 1
#define DEFAULT_OUTPUT_FILE "observables.dat"

#include <iostream>
#include <string>
#include <vector>

struct Parameters {
	std::string simulName;  // Type of simulation ('pipe', 'tonks')

	long nbParticles;  // Number of particles
	double density;  // Linear density
	double radExtra;  // Radius of the channel minus radius of particles (1.0)
	double length;  // Length of the channel (from nbParticles and density)

	double temperature;  // Temperature
	double eps;  // Potential strength
	double timestep;  // Time step
	long nbIters;  // Simulation duration
	long nbItersTh;  // Thermalization duration

	long nbTracers;  // Number of tracers
	std::vector<long> idTracers;  // Indices of the tracers
	std::vector<double> forces;  // Forces on the tracers
	
	long nbSimuls;  // Number of simulations
	int nbMoments;  // Number of moments to compute
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file

	bool checkOrder;  // Check the order of the particles at each iteration
	bool verbose;  // Verbose mode
};

int checkParameters(const Parameters &p);
void printParameters(const Parameters &p, std::ostream &stream = std::cout);

// Return 1 and print error if a is negative. Return 0 otherwise.
template<typename T>
int checkPositive(const T a, const std::string label) {
    if (a <= 0.) {
		std::cerr << "Error: " << label << " should be strictly positive."
            << std::endl;
        return 1;
    }
    return 0;
}

#endif
