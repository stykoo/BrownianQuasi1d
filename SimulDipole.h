#ifndef SIMUL_DIPOLE_H
#define SIMUL_DIPOLE_H

#include "parameters.h"
#include "simul.h"

class SimulDipole : public Simulation {
	public:
		SimulDipole(const Parameters &p);
		~SimulDipole();

	protected:
		// Positions of the particles
		std::vector<double> positions;
		// Forces between the particles
		std::vector<double> forces;

		// Methods to implement from Simulation
		void init(std::mt19937 &rndGen) override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
		double getPosX(const long i) override;

		// New methods
		void calcForcesBetweenParticles();
};

#endif
