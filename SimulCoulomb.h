#ifndef SIMUL_COULOMB_H
#define SIMUL_COULOMB_H

#include "parameters.h"
#include "simul.h"

class SimulCoulomb : public Simulation {
	public:
		SimulCoulomb(const Parameters &p);
		~SimulCoulomb();

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
