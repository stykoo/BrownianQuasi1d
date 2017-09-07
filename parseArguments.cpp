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
		 "Should be 'pipe', 'tonks' or 'canal'")
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
		("moments", po::value<int>(&p.nbMoments)->default_value(
			DEFAULT_NB_MOMENTS), "Number of moments to compute")
		("threads", po::value<int>(&p.nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output", po::value<std::string>(&p.output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")

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

