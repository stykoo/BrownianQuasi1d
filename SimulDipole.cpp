#include <algorithm>
#include "simul.h"
#include "SimulDipole.h"

SimulDipole::SimulDipole(const Parameters &p) : Simulation(p) {
}

SimulDipole::~SimulDipole() {
}

// Generate an initial state.
void SimulDipole::init(std::mt19937 &rndGen) {
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

	setInitXTracers();
}

// Set initial x-position of the tracers.
void SimulDipole::setInitXTracers() {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		initXTracers[i] = positions[p.idTracers[i]];
	}
}

// Implement one step of the time evolution of the system.
void SimulDipole::update(std::mt19937 &rndGen, const bool thermalization) {
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

// Compute the observables.
void SimulDipole::computeObservables(Observables &o) {
	std::vector<double> xsPer(p.nbTracers);
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		xsPer[i] = periodicBC(positions[p.idTracers[i]], p.length);
		// xsPer[i] = periodicBC(positions[p.idTracers[i]] - initXTracers[i],
		//                       p.length);
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

bool SimulDipole::isOrdered() {
	long c = 0;
	for (long i = 0 ; i < p.nbParticles-1 ; ++i) {
		if (positions[i] > positions[i+1]) {
			++c;
		}
	}
	if (positions[p.nbParticles-1] > positions[0]) {
		++c;
	}
	/* if (c >= 2) {
		std::cout << c << std::endl;
		for (long i=0 ; i<p.nbParticles ; ++i) {
			std::cout << positions[i] << " ";
		}
		std::cout << std::endl;
	} */
	return (c < 2);
}
