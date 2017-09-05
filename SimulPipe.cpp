#include <algorithm>
#include "simul.h"
#include "SimulPipe.h"

SimulPipe::SimulPipe(const Parameters &p) : Simulation(p) {
}

// Generate an initial state.
void SimulPipe::init(std::mt19937 &rndGen) {
	positions.resize(p.nbParticles);
	forces.resize(p.nbParticles);
	initXTracers.resize(p.nbTracers);

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

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<SIMUL_PIPE_DIM ; ++a) {
			forces[i][a] = 0;
		}
	}

	setInitXTracers();
}

// Set initial x-position of the tracers.
void SimulPipe::setInitXTracers() {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		initXTracers[i] = positions[p.idTracers[i]][0];
	}
}

// Implement one step of the time evolution of the system.
void SimulPipe::update(std::mt19937 &rndGen, const bool thermalization) {
	std::normal_distribution<double> rndForNoise(0., sqrt(2. * p.temperature
														  * p.timestep));
	
	// Old position
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<SIMUL_PIPE_DIM-1 ; ++a) {
			oldPosYZ[i][a] = positions[i][a+1];
		}
	}
	
	calcForcesBetweenParticles();
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<SIMUL_PIPE_DIM ; ++a) {
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

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void SimulPipe::calcForcesBetweenParticles() {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dr[SIMUL_PIPE_DIM];
		double distsq = 0.;
		for (int a=0 ; a<SIMUL_PIPE_DIM ; ++a) {
			dr[a] = positions[i][a] - positions[iPrev][a];
			distsq += dr[a] * dr[a];
			forces[i][a] = 0;
		}
		if (distsq < 1. && distsq > 0.) {
			for (int a=0 ; a<SIMUL_PIPE_DIM ; ++a) {
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

// Compute the observables.
void SimulPipe::computeObservables(Observables &o) {
	std::vector<double> xsPer(p.nbTracers);
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		xsPer[i] = periodicBC(positions[i][0] - initXTracers[i],
							  p.length);
	}
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			o.moments[i][j] = 1;
			for (long k = j ; k < j + i + 1 ; ++k) {
				o.moments[i][j] *= xsPer[k];
			}
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
