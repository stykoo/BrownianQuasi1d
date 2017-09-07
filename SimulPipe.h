#ifndef SIMUL_PIPE_H
#define SIMUL_PIPE_H

#define DIM_PIPE 3

#include "parameters.h"
#include "simul.h"

class SimulPipe : public Simulation {
	public:
		SimulPipe(const Parameters &p);
		~SimulPipe(){}

	protected:
		// Positions of the particles
		std::vector< std::array<double, DIM_PIPE> > positions;
		// Old positions in Y/Z
		std::vector< std::array<double, DIM_PIPE-1> > oldPosYZ;
		// Forces between the particles
		std::vector< std::array<double, DIM_PIPE> > forces;

		// Methods to implement from Simulation
		void init(std::mt19937 &rndGen) override;
		void update(std::mt19937 &rndGen, const bool thermalization) override;
		double getPosX(const long i) override;

		// New methods
		void calcForcesBetweenParticles();
		void keepInChannel();
};

// Other functions that are helpful
void reflexionInCircle(const double xIn, const double yIn,
		               const double xOut, const double yOut,
					   const double R,
					   double &xFin, double &yFin);
void findIntersection(const double xIn, const double yIn,
		              const double xOut, const double yOut,
					  const double R,
					  double &xCross, double &yCross);
void basicReflexion(const double ux, const double uy,
		            double normalX, double normalY,
					double &xRefl, double &yRefl);
void solveSecondOrderEq(const double a, const double b, const double c,
                        double &sol1, double &sol2);

#endif
