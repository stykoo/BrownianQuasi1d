/*
Copyright Université Pierre et Marie Curie (2017)
Contributors: Alexis Poncet

alexis.poncet@ens.fr

This software is a computer program whose purpose is to simulate the
motion of interacting Brownian particles in a quasi-1d environment,
some of them being driven by an external force.

This software is governed by the CeCILL license under French law and
abiding by the rules of distribution of free software.  You can  use,
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info".

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability.

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or
data to be ensured and,  more generally, to use and operate it in the
same conditions as regards security.

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/
/*
 * BrownianQuasi1d
 * SimulPipe.h
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Definitions related to the simulation in a pipe.
 */

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
