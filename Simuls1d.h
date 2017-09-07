#ifndef SIMULS_1D_H
#define SIMULS_1D_H

#include "parameters.h"
#include "simul.h"

class Simul1d : public Simulation {
	public:
		Simul1d(const Parameters &p) : Simulation(p) {}
		virtual ~Simul1d() {}

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
		virtual void calcForcesBetweenParticles() = 0;
};

class SimulTonks : public Simul1d {
	public:
		SimulTonks(const Parameters &p) : Simul1d(p) {}
		~SimulTonks() {}

	protected:
		void calcForcesBetweenParticles() override;
};

class SimulCoulomb : public Simul1d {
	public:
		SimulCoulomb(const Parameters &p) : Simul1d(p) {}
		~SimulCoulomb() {}

	protected:
		void calcForcesBetweenParticles() override;
};

class SimulDipole : public Simul1d {
	public:
		SimulDipole(const Parameters &p) : Simul1d(p) {}
		~SimulDipole() {}

	protected:
		void calcForcesBetweenParticles() override;
};

#endif
