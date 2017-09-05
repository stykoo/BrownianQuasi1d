#include <algorithm>
#include <fstream>
#include <map>
#include <chrono>
#include <thread>
#include <iomanip>
#include "simul.h"

// Run all the simulations and store the moments.
int runSimulations(const Parameters &p) {
	std::random_device rd;
	std::vector<std::thread> threads;
	std::vector< std::vector<Observables> > allSumsObs(p.nbThreads);

	long nbSimulsPerThread = p.nbSimuls / p.nbThreads;
	
	if (p.nbSimuls % p.nbThreads != 0) {
		std::cerr << "Warning: nbSimuls is not a multiple of nbThreads."
			<< nbSimulsPerThread * p.nbThreads << " simulations will be done."
			<< std::endl;
	}

	// Threads
	for (int i = 0 ; i < p.nbThreads ; ++i) {
		threads.push_back(std::thread(runMultipleSimulations, p,
			nbSimulsPerThread, std::ref(allSumsObs[i]), rd())); 
	}
	
	// Wait for everyone
	for (auto &th : threads) {
		th.join();
	}

	// Initialize the total sum
	std::vector<Observables> sumObs;
	initObservables(sumObs, p);

	// Add the observables to the total sum
	for (int k = 0 ; k < p.nbThreads ; ++k) {
		addObservables(sumObs, allSumsObs[k], p);
	}

	int status = exportObservables(sumObs, p);

	return status;
}

// Run several simulations and store the sum of the observables.
// This function is usually called as a thread.
void runMultipleSimulations(const Parameters &p, const long nbSimuls,
						   std::vector<Observables> &sumObs,
						   const unsigned int seed) {
    // Random generator
	std::mt19937 rndGen(seed);

	// Initialize sumObs
	initObservables(sumObs, p);

	for (long s = 0 ; s < nbSimuls ; ++s) {
		if (p.verbose) {
			std::cout << "Simulation " << s + 1 
				<< " on thread " << std::this_thread::get_id() << std::endl;
		}

		// Run a single simulation
		std::vector<Observables> obs;
		runOneSimulation(p, obs, rndGen);

		// Add the observables to the sum
		addObservables(sumObs, obs, p);
	}
}

// Run a single simulation and compute the moments.
void runOneSimulation(const Parameters &p, std::vector<Observables> &obs,
					 std::mt19937 &rndGen) {
	// Positions of the tracers
	State state;
	initState(state, p, rndGen);

	// Thermalization
	for (long j=0 ; j<p.nbItersTh ; ++j) {
		updateState(state, p, rndGen, true);
	}

	// Initial observables
	setInitXTracers(state, p);
	initObservables(obs, p);
	computeObservables(state, p, obs[0]);

	// Main loop
	for (long j=0 ; j<p.nbIters-1 ; ++j) {
		updateState(state, p, rndGen, false);
		computeObservables(state, p, obs[j+1]);
	}
}

// Generate an initial state.
void initState(State &state, const Parameters &p, std::mt19937 &rndGen) {
	state.positions.resize(p.nbParticles);
	state.forces.resize(p.nbParticles);
	state.initXTracers.resize(p.nbTracers);

	std::uniform_real_distribution<double> distrib(0., 1.);

	for (auto &pos : state.positions) {
		pos[0] = p.length * (distrib(rndGen) - 0.5);
		double r = p.radExtra * distrib(rndGen);
		double theta = 2. * M_PI * distrib(rndGen);
		pos[1] = r * std::cos(theta);
		pos[2] = r * std::sin(theta);
	}

	// Sort by x value
	std::sort(state.positions.begin(), state.positions.end(),
			  [](auto const &a, auto const &b) {
				 return a.front() < b.front();
			  });

	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM ; ++a) {
			state.forces[i][a] = 0;
		}
	}

	setInitXTracers(state, p);
}

// Set initial x-position of the tracers.
void setInitXTracers(State &state, const Parameters &p) {
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		state.initXTracers[i] = state.positions[p.idTracers[i]][0];
	}
}

// Implement one step of the time evolution of the system.
void updateState(State &state, const Parameters &p, std::mt19937 &rndGen,
				 const bool thermalization) {
	std::normal_distribution<double> rndForNoise(0., sqrt(2. * p.temperature
														  * p.timestep));
	
	// Old position
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM-1 ; ++a) {
			state.oldPosYZ[i][a] = state.positions[i][a+1];
		}
	}
	
	calcForcesBetweenParticles(state, p);
	
	for (long i=0 ; i<p.nbParticles ; ++i) {
		for (int a=0 ; a<DIM ; ++a) {
			// Noise
			state.positions[i][a] += rndForNoise(rndGen);

			// Forces between particles
			state.positions[i][a] += p.timestep * state.forces[i][a];
		}
	}

	// Forces on the tracers
	if (!thermalization) {
		for (long i=0 ; i<p.nbTracers ; ++i) {
			state.positions[p.idTracers[i]][0] += p.timestep * p.forces[i];
		}
	}

	// Keep all the particles in the channel
	keepInChannel(state, p);
}

// Compute the forces between the particles.
// WE ASSUME THAT THE PARTICLES ARE ORDERED AND NEVER CROSS.
void calcForcesBetweenParticles(State &state, const Parameters &p) {
	for (long i=0 ; i<p.nbParticles ; ++i) {
		long iPrev = (i + p.nbParticles - 1) % p.nbParticles;

		double dr[DIM];
		double distsq = 0.;
		for (int a=0 ; a<DIM ; ++a) {
			dr[a] = state.positions[i][a] - state.positions[iPrev][a];
			distsq += dr[a] * dr[a];
			state.forces[i][a] = 0;
		}
		if (distsq < 1. && distsq > 0.) {
			for (int a=0 ; a<DIM ; ++a) {
				double f = p.eps * (1. / sqrt(distsq) - 1.) * dr[a];
			state.forces[i][a] += f;
			state.forces[iPrev][a] -= f;
			}
		}
	}
}

void keepInChannel(State &state, const Parameters &p) {
	// Periodic boundary conditions in X
	for (long i=0 ; i<p.nbParticles ; ++i) {
		state.positions[i][0] = periodicBC(state.positions[i][0], p.length);
	}

	// Circular wall in Y/Z
	for (long i=0 ; i<p.nbParticles ; ++i) {
		double radsq = (state.positions[i][1] * state.positions[i][1]
				        + state.positions[i][2] * state.positions[i][2]);

		if(radsq > p.radExtra * p.radExtra) {
			double yNew, zNew;
			reflexionInCircle(state.oldPosYZ[i][0], state.oldPosYZ[i][1],
					          state.positions[i][1], state.positions[i][2],
							  p.radExtra, yNew, zNew);
			state.positions[i][1] = yNew;
			state.positions[i][2] = zNew;
		}
	}
}

// Initialize a vector of observables
void initObservables(std::vector<Observables> &obs, const Parameters &p) {
	obs.resize(p.nbIters);
	for (long t = 0 ; t < p.nbIters ; ++t) {
		obs[t].moments.resize(p.nbTracers);
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			obs[t].moments[i].assign(p.nbTracers - i, 0);
		}
	}
}

// Compute the observables (except those for the occupations).
void computeObservables(const State &state, const Parameters &p,
						Observables &o) {
	std::vector<double> xsPer(p.nbTracers);
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		xsPer[i] = periodicBC(state.positions[i][0] - state.initXTracers[i],
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

// Add observables o2 to observables o1.
void addObservables(std::vector<Observables> &obs1,
					const std::vector<Observables> &obs2, const Parameters &p)
{
	for (long t = 0 ; t < p.nbIters ; ++t) {
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				obs1[t].moments[i][j] += obs2[t].moments[i][j];
			}
		}
	}
}

// Export the observables to a file.
int exportObservables(const std::vector<Observables> &sumObs,
					  const Parameters &p) {
	std::ofstream file(p.output);
	if (!file.is_open()) {
		return 1;
	}

	// Header
	file << "# BrownianQuasi1d  (" << __DATE__ <<  ", " << __TIME__
		<< "): ";
	printParameters(p, file);
	file << "\n# t";
	for (long i = 0 ; i < p.nbTracers ; ++i) {
		for (long j = 0 ; j < p.nbTracers - i ; ++j) {
			file << " ";
			for (long k = j ; k < j + i + 1 ; ++k) {
				file << "x" << k+1;
			}
		}
	}
	file << "\n";

	file << std::scientific << std::setprecision(DEFAULT_OUTPUT_PRECISION);

	// Data (we write the average and not the sum)
	for (long k = 0 ; k < p.nbIters ; ++k) {
		file << k * p.timestep;
		for (long i = 0 ; i < p.nbTracers ; ++i) {
			for (long j = 0 ; j < p.nbTracers - i ; ++j) {
				file << " " << sumObs[k].moments[i][j] / p.nbSimuls;
			}
		}
		file << "\n";
	}

	file.close();

	return 0;
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

// Return the sign of a number (+1 if positive or null, -1 if negative)
int sign(const double x) {
    return (0. <= x) - (x < 0.);
}
