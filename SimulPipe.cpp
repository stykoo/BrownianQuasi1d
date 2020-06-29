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
 * SimulPipe.cpp
 *
 * Author: Alexis Poncet
 * Email: alexis.poncet@ens.fr
 *
 * Implementation of the simulation in a pipe
 * (3d system confined into a cylinder).
 */

#include <algorithm>
#include "simul.h"
#include "SimulPipe.h"

SimulPipe::SimulPipe(const Parameters &p) : Simulation(p) {
}

// Generate an initial state.
int SimulPipe::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM_PIPE ; ++a) {
			forces[i][a] = 0;
		}
	}

	std::uniform_real_distribution<double> distrib(0., 1.);

	for (auto &pos : positions) {
		pos[0] = p.length * (distrib(rndGen) - 0.5);
		double r = p.radExtra * distrib(rndGen);
		double theta = 2. * M_PI * distrib(rndGen);
		pos[1] = r * std::cos(theta);
		pos[2] = r * std::sin(theta);
	}

	// Sort by x value
	std::sort(positions.begin(), positions.end(),
			  [](auto const &a, auto const &b) {
				 return a.front() < b.front();
			  });

	// If the potential is strong, the order of the particles may not be
	// conserved in the first iterations: we do some thermalization.
	update(rndGen, true);
	for (int i = 0 ; i < MAX_ITERS_INIT ; ++i) {
		if (isOrdered()) {
			return 0;
		}
		std::sort(positions.begin(), positions.end(),
				  [](auto const &a, auto const &b) {
					 return a.front() < b.front();
				  });
		update(rndGen, true);
	}
	return 1;
}

// Implement one step of the time evolution of the system.
void SimulPipe::update(std::mt19937 &rndGen, const bool thermalization) {
	std::normal_distribution<double> rndForNoise(0., sqrt(2. * p.temperature
														  * p.timestep));
	
	// Old position
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM_PIPE-1 ; ++a) {
			oldPosYZ[i][a] = positions[i][a+1];
		}
	}
	
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM_PIPE ; ++a) {
			// Noise
			positions[i][a] += rndForNoise(rndGen);

			// Forces between particles
			positions[i][a] += p.timestep * forces[i][a];
		}
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			positions[p.idTracers[i]][0] += p.timestep * p.forces[i];
		}
	}

	// Keep all the particles in the channel
	keepInChannel();
}

// Get position in X of particle i
double SimulPipe::getPosX(const long i) {
	return positions[i][0];
}

// Get profile
void SimulPipe::getPosRel(std::vector<double> &posr, const long nbPts) {
	double x;
	long k;
	for (long j = 0 ; j < nbPts ; ++j) {
		posr[j] = 0.0;
	}
	for (long i = 1 ; i < p.nbParticles ; ++i) {
		x = periodicBCpos(positions[i][0] - positions[0][0], p.length);
		k = (long) (x / p.length * nbPts);
		posr[k] += 1.0;
	}
}

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void SimulPipe::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dr[DIM_PIPE];
		double distsq = 0.;
		for (int a=0 ; a<DIM_PIPE ; ++a) {
			dr[a] = positions[i][a] - positions[iPrev][a];
			distsq += dr[a] * dr[a];
			forces[i][a] = 0;
		}
		if (distsq < 1. && distsq > 0.) {
			for (int a=0 ; a<DIM_PIPE ; ++a) {
				double f = p.eps * (1. / sqrt(distsq) - 1.) * dr[a];
				forces[i][a] += f;
				forces[iPrev][a] -= f;
			}
		}
	}
}

void SimulPipe::keepInChannel() {
	// Periodic boundary conditions in X
	for (long i=0 ; i<p.nbParticles ; ++i) {
		positions[i][0] = periodicBC(positions[i][0], p.length);
	}

	// Circular wall in Y/Z
	for (long i=0 ; i<p.nbParticles ; ++i) {
		double radsq = (positions[i][1] * positions[i][1]
				        + positions[i][2] * positions[i][2]);

		if(radsq > p.radExtra * p.radExtra) {
			double yNew, zNew;
			reflexionInCircle(oldPosYZ[i][0], oldPosYZ[i][1],
					          positions[i][1], positions[i][2],
							  p.radExtra, yNew, zNew);
			positions[i][1] = yNew;
			positions[i][2] = zNew;
		}
	}
}

// A ray going from (xIn, yIn) to (xOut, yOut) is reflected inside
// a circle of center the origin and of radius R.
// This computes (xFin, yFin) the final point of the ray inside the circle.
// WE ASSUME (xIn, yIn) IS INSIDE THE CIRCLE AND (xOut, yOut) IS OUTSIDE
void reflexionInCircle(const double xIn, const double yIn,
		               const double xOut, const double yOut,
					   const double R,
					   double &xFin, double &yFin) {
	// Find the point of intersection between the ray and the circle.
	double xCross = xIn, yCross = yIn;  // Arbitrary
	findIntersection(xIn, yIn, xOut, yOut, R, xCross, yCross);

	// Do the reflexion with respect to the tangent to the circle.
	double xRefl = xCross, yRefl = yCross;  // Arbitrary
	basicReflexion(xOut-xCross, yOut-yCross, xCross, yCross, xRefl, yRefl);

	// Add the reflected point to the intersection point.
	xFin = xCross + xRefl;
	yFin = yCross + yRefl;
}

// Find the point of intersection between the ray and the circle.
// WE ASSUME IT EXISTS AND IS UNIQUE
void findIntersection(const double xIn, const double yIn,
		              const double xOut, const double yOut,
					  const double R,
					  double &xCross, double &yCross) {
	double dx = xOut - xIn;
	double dy = yOut - yIn;
	double a = dx*dx + dy*dy;
	double b = 2 * (xIn*dx + yIn*dy);
	double c = xIn*xIn + yIn*yIn - R*R;

	double t1 = 0, t2 = 0;  // Arbitrary
	solveSecondOrderEq(a, b, c, t1, t2);

	if (t1 >= 0 && t1 <= 1) {
		xCross = xIn + t1 * dx;
		yCross = yIn + t1 * dy;
	} else if (t2 >= 0 && t2 <= 1) {
		xCross = xIn + t2 * dx;
		yCross = yIn + t2 * dy;
	} else {
		return;
	}
}

// Do the reflexion of a vector (ux, uy) given a normal vector to the line
// (not necessarily normalized).
void basicReflexion(const double ux, const double uy,
		            double normalX, double normalY,
					double &xRefl, double &yRefl) {
	double nn = std::sqrt(normalX*normalX + normalY*normalY);
	normalX /= nn;
	normalY /= nn;

	double tgtX = -normalY;
	double tgtY = normalX;

	double prodScalNorm = ux*normalX + uy*normalY;
	double prodScalTgt = ux*tgtX + uy*tgtY;

	xRefl = -prodScalNorm * normalX + prodScalTgt * tgtX;
	yRefl = -prodScalNorm * normalY + prodScalTgt * tgtY;
}

// Solve the equation ax^2 + bx + c = 0 ASSUMING A SOLUTION EXISTS
void solveSecondOrderEq(const double a, const double b, const double c,
                        double &sol1, double &sol2) {
	if (a == 0.) {
		sol1 = -c/b;
		sol2 = -c/b;
	} else {
		double d = b*b - 4*a*c;
		if(d < 0) {
			return;
		}
		// Trick for numerical stability
		double tmp = -0.5 * (b + sign(b) * std::sqrt(d));
		sol1 = tmp / a;
		sol2 = c / tmp;
	}
}
