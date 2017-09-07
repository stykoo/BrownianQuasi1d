#include <algorithm>
#include "simul.h"
#include "SimulTonks.h"

SimulTonks::SimulTonks(const Parameters &p) : Simulation(p) {
}

SimulTonks::~SimulTonks() {
}

// Generate an initial state.
void SimulTonks::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = p.length * (distribUnif(rndGen) - 0.5);
		forces[i] = 0;  // Arbitrary
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end());

	// If the potential is strong, the order of the particles may not be
	// conserved in the first iterations: we do some thermalization.
	update(rndGen, true);
	while (!isOrdered()) {
		std::sort(positions.begin(), positions.end());
		update(rndGen, true);
	}
}

// Implement one step of the time evolution of the system.
void SimulTonks::update(std::mt19937 &rndGen, const bool thermalization) {
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		// Noise
		positions[i] += noise * distribNormal(rndGen);

		// Forces between particles
		positions[i] += p.timestep * forces[i];
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			positions[p.idTracers[i]] += p.timestep * p.forces[i];
		}
	}
	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = periodicBC(positions[i], p.length);
	}
}

// Get position in X of particle i
double SimulTonks::getPosX(const long i) {
	return positions[i];
}

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void SimulTonks::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dx = periodicBC(positions[i] - positions[iPrev], p.length);
		if (dx < 1. && dx > 0.) {
			double f = p.eps * (1. - dx);
			forces[i] += f;
			forces[iPrev] -= f;
		}
	}
}
