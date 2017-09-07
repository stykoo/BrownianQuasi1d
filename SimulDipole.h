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
		// Positions of the tracers after thermalization
		std::vector<double> initXTracers;


		// Methods to implement from Simulation
		void init(std::mt19937 &rndGen) override;
		void setInitXTracers() override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
		void computeObservables(Observables &o) override;
		bool isOrdered() override;

		// New methods
		void calcForcesBetweenParticles();
};

#endif
