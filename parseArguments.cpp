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
		("particles,n", po::value<long>(&p.nbParticles)->required(),
		 "Number of particles")
		("density,r", po::value<double>(&p.density)->required(),
		 "Linear density")
		("radExtra,R", po::value<double>(&p.radExtra)->required(),
		 "Diameter of the channel")

		("temperature,T", po::value<double>(&p.temperature)->required(),
		 "Temperature")
		("eps,e", po::value<double>(&p.eps)->required(),
		 "Strength of the inter-particle potential")
		("timestep,t", po::value<double>(&p.timestep)->required(), "Timestep")
		("nbIters,I", po::value<long>(&p.nbIters)->required(),
		 "Number of iterations")
		("nbItersTh,J", po::value<long>(&p.nbItersTh)->required(),
		 "Number of iterations of thermalization")

		("ids,i",
		 po::value< std::vector<long> >(&p.idTracers)->multitoken()
		                                              ->required(),
		 "Initial relative positions of the tracers")
		("forces,f",
		 po::value< std::vector<double> >(&p.forces)->multitoken()->required(),
		 "Forces on the tracers.")

		("simuls,s", po::value<long>(&p.nbSimuls)->required(),
		 "Number of repetitions of the simulation")
		("moments,M", po::value<int>(&p.nbMoments)->default_value(
			DEFAULT_NB_MOMENTS), "Number of moments to compute")
		("threads,c", po::value<int>(&p.nbThreads)->default_value(
			DEFAULT_THREADS), "Number of threads")
		("output,o", po::value<std::string>(&p.output)->default_value(
			DEFAULT_OUTPUT_FILE), "Output file")

        ("verbose,v", po::bool_switch(&p.verbose), "Verbose mode")
        ("help,h", "Print help message and exit")
		;

	try {
		po::variables_map vars;
		po::store(po::parse_command_line(argc, argv, opts), vars);

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

