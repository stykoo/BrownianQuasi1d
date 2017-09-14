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
 * parseArguments.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Parse the command-line arguments.
 */

#include <iostream>
#include <exception>
#include <boost/program_options.hpp>
#include "parseArguments.h"

namespace po = boost::program_options;

// Parse command-line arguments and store the values into Parameters.
// Return 1 if displaying help, 2 if a problem occurs, 0 otherwise.
int parseArguments(int argc, char **argv, Parameters &p) {
	po::options_description opts("Options");
	opts.add_options()
		("name", po::value<std::string>(&p.simulName)->required(),
		 "Should be 'tonks', 'coulomb', 'dipole', 'canal' or 'pipe'")
		("particles", po::value<long>(&p.nbParticles)->required(),
		 "Number of particles")
		("density", po::value<double>(&p.density)->required(),
		 "Linear density")
		("radExtra", po::value<double>(&p.radExtra)->default_value(0.1),
		 "Diameter of the channel")

		("temperature", po::value<double>(&p.temperature)->default_value(1.),
		 "Temperature")
		("eps", po::value<double>(&p.eps)->required(),
		 "Strength of the inter-particle potential")
		("timestep", po::value<double>(&p.timestep)->required(), "Timestep")
		("nbIters", po::value<long>(&p.nbIters)->required(),
		 "Number of iterations")
		("nbItersTh", po::value<long>(&p.nbItersTh)->default_value(0),
		 "Number of iterations of thermalization")

		("ids",
		 po::value< std::vector<long> >(&p.idTracers)->multitoken()
		                                              ->required(),
		 "Initial relative positions of the tracers")
		("forces",
		 po::value< std::vector<double> >(&p.forces)->multitoken()->required(),
		 "Forces on the tracers.")

		("simuls", po::value<long>(&p.nbSimuls)->default_value(1),
		 "Number of repetitions of the simulation")
		("threads", po::value<int>(&p.nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output", po::value<std::string>(&p.output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")
		("skip", po::value<long>(&p.skip)->default_value(1),
		 "Compute observables every given number of iterations")

        ("checkOrder", po::bool_switch(&p.checkOrder),
		 "Check the order of the particles at each iteration")
        ("verbose", po::bool_switch(&p.verbose), "Verbose mode")
        ("help", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts,
					                     po::command_line_style::unix_style ^
										 po::command_line_style::allow_short),
				  vars);

        // Display help and exit
        if (vars.count("help")) {
			std::cout << "Usage: " << argv[0] << " options\n";
			std::cout << opts << std::endl;
            return 1;
        }

        po::notify(vars);

		// Compute here the parameters that should be computed.
		p.length = p.nbParticles / p.density;
		p.nbTracers = (long) p.idTracers.size();
	} catch (std::exception &e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 2;
	}

	return 0;
}

