#ifndef SIMUL_TONKS_H
#define SIMUL_TONKS_H

#include "parameters.h"
#include "simul.h"

class SimulTonks : public Simulation {
	public:
		SimulTonks(const Parameters &p);
		~SimulTonks();

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
