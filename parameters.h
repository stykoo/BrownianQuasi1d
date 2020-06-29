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
 * parameters.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Contain the definition of the parameters.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#define DEFAULT_NB_MOMENTS 4
#define DEFAULT_OUTPUT_PRECISION 15
#define DEFAULT_THREADS 1
#define DEFAULT_PREC_PROF 0.1
#define DEFAULT_OUTPUT_FILE "observables.dat"

#include <iostream>
#include <string>
#include <vector>

const std::vector<std::string> NAMES = {"tonks", "canal", "pipe", "coulomb",
	                                    "dipole", "coulCircle", "dipCircle" };

struct Parameters {
	// Type of simulation ('pipe', 'tonks', 'canal', 'coulomb', 'dipole')
	std::string simulName;

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
	int nbThreads;  // Number of threads
	std::string output;  // Name of the output file
	long skip;  // Compute observables every given number of iterations

	bool computeProfs; // Compute the profiles
	double precProfs; // Spatial resolution for the profiles
	long nbPtsProfs; // Number of points for the profiles

	bool checkOrder;  // Check the order of the particles at each iteration
	bool verbose;  // Verbose mode
};

int checkParameters(const Parameters &p);
int checkSimulName(const std::string name);
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
