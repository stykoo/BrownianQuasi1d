#include <algorithm>
#include "simul.h"
#include "SimulTonks.h"

SimulTonks::SimulTonks(const Parameters &p) : Simulation(p) {
}

// Generate an initial state.
void SimulTonks::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

	std::uniform_real_distribution<double> distrib(0., 1.);

	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i] = p.length * (distrib(rndGen) - 0.5);
		forces[i] = 0;  // Arbitrary
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end());
	for (long i=0 ; i<p.nbParticles ; ++i) {
		std::cout << positions[i] << " ";
	}
	std::cout << std::endl;

	setInitXTracers();
}

// Set initial x-position of the tracers.
void SimulTonks::setInitXTracers() {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		initXTracers[i] = positions[p.idTracers[i]];
	}
}

// Implement one step of the time evolution of the system.
void SimulTonks::update(std::mt19937 &rndGen, const bool thermalization) {
	std::normal_distribution<double> rndForNoise(0., sqrt(2. * p.temperature
														  * p.timestep));
	
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		// Noise
		positions[i] += rndForNoise(rndGen);

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
		std::cout << positions[i] << " "; // TODO
	}
	std::cout << std::endl;
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

// Compute the observables.
void SimulTonks::computeObservables(Observables &o) {
	std::vector<double> xsPer(p.nbTracers);
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		// TODO
		// xsPer[i] = periodicBC(positions[i] - initXTracers[i], p.length);
		xsPer[i] = periodicBC(positions[i], p.length);
	}
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			o.moments[i][j] = 1;
			for (long k = j ; k < j + i + 1 ; ++k) {
				o.moments[i][j] *= xsPer[k];
			}
		}
	}
}
