#include "parameters.h"

// Check if the parameters are valid. Return 0 if they are, 1 otherwise.
int checkParameters(const Parameters &p) {
	if (p.simulName != "pipe" && p.simulName != "tonks") {
		std::cerr << "Wrong simulation name (" << p.simulName << "). "
			<< "Only 'pipe' and 'tonks' are allowed." << std::endl;
	}

    if (checkPositive(p.nbParticles, "nbParticles") ||
        checkPositive(p.density, "density") ||
        checkPositive(p.radExtra, "radExtra") ||
        checkPositive(p.length, "length") ||
        checkPositive(p.temperature, "temperature") ||
        checkPositive(p.eps, "eps") ||
        checkPositive(p.timestep, "timestep") ||
        checkPositive(p.nbIters, "nbIters") ||
        checkPositive(p.nbTracers, "nbTracers") ||
        checkPositive(p.nbTracers, "nbSimuls") ||
        checkPositive(p.nbTracers, "nbMoments") ||
        checkPositive(p.nbMoments, "nbThreads")) {
        return 1;
    }

	if (p.nbTracers > p.nbParticles) {
		std::cerr << "Error: The number of tracers should be smaller than"
			<< " the number of particles." << std::endl;
		return 1;
	}
	if (p.radExtra >= 1.) {
		std::cerr << "Error: The extra radius should be smaller than 1."
			<< std::endl;
		return 1;
	}
	if ((long) p.idTracers.size() != p.nbTracers) {
		std::cerr << "Critical error: the size of the vector of ids"
			<< " is wrong." << std::endl;
		return 1;
	}
	for (auto po : p.idTracers) {
		if (po < 0 || po >= p.nbParticles) {
			std::cerr << "Error: the positions should in the right range."
				<< std::endl;
			return 1;
		}
	}
	for (long i = 0 ; i < p.nbTracers - 1 ; ++i) {
		if (p.idTracers[i] >= p.idTracers[i+1]) {
			std::cerr << "Error: please sort the positions."
				<< std::endl;
			return 1;
		}
	}
	if ((long) p.forces.size() != p.nbTracers) {
		std::cerr << "Error: the number of forces and the number"
			<< " of tracers should be equal." << std::endl;
		return 1;
	}
	return 0;
}

// Print the parameters to stream.
void printParameters(const Parameters &p, std::ostream &stream) {
	stream << "_" << p.simulName << "_"
		<< ", particles=" << p.nbParticles << ", density=" << p.density
		<< ", radExtra=" << p.radExtra << ", length=" << p.length
		<< ", temperature=" << p.temperature << ", eps=" << p.eps
	   	<< ", timestep=" << p.timestep
		<< ", nbIters=" << p.nbIters << ", nbItersTh=" << p.nbItersTh
		<< ", nbSimuls=" << p.nbSimuls << ", nbMoments=" << p.nbMoments;
	stream << ", idTracers=";
	for (auto po : p.idTracers) {
		stream << po << ":";
	}
	stream << ", forces=";
	for (auto f : p.forces) {
		stream << f << ":";
	}
}

