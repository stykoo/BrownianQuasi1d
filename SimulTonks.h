#ifndef SIMUL_TONKS_H
#define SIMUL_TONKS_H

#include "parameters.h"
#include "simul.h"

class SimulTonks : public Simulation {
	public:
		SimulTonks(const Parameters &p);
		~SimulTonks(){}

	protected:
		// Positions of the particles
		std::vector<double> positions;
		// Forces between the particles
		std::vector<double> forces;
		// Positions of the tracers after thermalization
		std::vector<double> initXTracers;

		// Methods to implement from Simulation
		void init(std::mt19937 &rndGen) override;
		void setInitXTracers() override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
		void computeObservables(Observables &o) override;

		// New methods
		void calcForcesBetweenParticles();
};

#endif
