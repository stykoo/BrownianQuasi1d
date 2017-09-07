#ifndef SIMUL_CANAL_H
#define SIMUL_CANAL_H

#include "parameters.h"
#include "simul.h"

#define DIM_CANAL 2

class SimulCanal : public Simulation {
	public:
		SimulCanal(const Parameters &p);
		~SimulCanal();

	protected:
		// Positions of the particles
		std::vector< std::array<double, DIM_CANAL> > positions;
		// Forces between the particles
		std::vector< std::array<double, DIM_CANAL> > forces;

		// Methods to implement from Simulation
		void init(std::mt19937 &rndGen) override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
		double getPosX(const long i) override;

		// New methods
		void calcForcesBetweenParticles();
		void keepInCanal();
};

#endif
