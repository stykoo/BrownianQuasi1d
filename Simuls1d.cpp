#include <algorithm>
#include "simul.h"
#include "Simuls1d.h"

// Generate an initial state.
void Simul1d::init(std::mt19937 &rndGen) {
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
void Simul1d::update(std::mt19937 &rndGen, const bool thermalization) {
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
double Simul1d::getPosX(const long i) {
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

// Compute the forces between the particles.
void SimulCoulomb::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double dx = periodicBC(positions[i] - positions[j], p.length);
			if (dx != 0.) {
				double f = sign(dx) * p.eps / (dx * dx);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}

// Compute the forces between the particles.
void SimulDipole::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		forces[i] = 0;
	}

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (long j=i+1 ; j<p.nbParticles ; ++j) {
			double dx = periodicBC(positions[i] - positions[j], p.length);
			if (dx != 0.) {
				double f = sign(dx) * 3. * p.eps / mypow(dx, 4);
				forces[i] += f;
				forces[j] -= f;
			}
		}
	}
}
