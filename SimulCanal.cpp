#include <algorithm>
#include "simul.h"
#include "SimulCanal.h"

SimulCanal::SimulCanal(const Parameters &p) : Simulation(p) {
}

SimulCanal::~SimulCanal() {
}

// Generate an initial state.
void SimulCanal::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i][0] = p.length * (distribUnif(rndGen) - 0.5);
		positions[i][1] = p.radExtra * (distribUnif(rndGen) - 0.5);
		forces[i][0] = 0;  // Arbitrary
		forces[i][1] = 0;  // Arbitrary
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end(),
			  [](auto const &a, auto const &b) {
				 return a.front() < b.front();
			  });

	// If the potential is strong, the order of the particles may not be
	// conserved in the first iterations: we do some thermalization.
	update(rndGen, true);
	while (!isOrdered()) {
		std::sort(positions.begin(), positions.end(),
				  [](auto const &a, auto const &b) {
					 return a.front() < b.front();
				  });
		update(rndGen, true);
	}

	setInitXTracers();
}

// Implement one step of the time evolution of the system.
void SimulCanal::update(std::mt19937 &rndGen, const bool thermalization) {
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM_CANAL ; ++a) {
			// Noise
			positions[i][a] += noise * distribNormal(rndGen);

			// Forces between particles
			positions[i][a] += p.timestep * forces[i][a];
		}
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			positions[p.idTracers[i]][0] += p.timestep * p.forces[i];
		}
	}

	keepInCanal();
}

// Get position in X of particle i
double SimulCanal::getPosX(const long i) {
	return positions[i][0];
}

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void SimulCanal::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dr[DIM_CANAL];
		double distsq = 0.;
		for (int a=0 ; a<DIM_CANAL ; ++a) {
			dr[a] = positions[i][a] - positions[iPrev][a];
			distsq += dr[a] * dr[a];
			forces[i][a] = 0;
		}
		if (distsq < 1. && distsq > 0.) {
			for (int a=0 ; a<DIM_CANAL ; ++a) {
				double f = p.eps * (1. / sqrt(distsq) - 1.) * dr[a];
				forces[i][a] += f;
				forces[iPrev][a] -= f;
			}
		}
	}
}

void SimulCanal::keepInCanal() {
	// PBC in x
	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i][0] = periodicBC(positions[i][0], p.length);
	}

	// Reflexion in y
	for (long i=0 ; i<p.nbParticles ; ++i) {
		// This should do an arbitrary number of reflexions (to be checked!)
		double t = (positions[i][1] + p.radExtra) / (2.*p.radExtra);
		int k = ((int) std::floor(t)) % 2;
		positions[i][1] = k * periodicBC(positions[i][1], 2.*p.radExtra);
	}
}
