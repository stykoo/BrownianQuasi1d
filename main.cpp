#include <iostream>
#include "parameters.h"
#include "parseArguments.h"
#include "simul.h"

int main(int argc, char **argv) {
	// Load parameters
    Parameters params;
    int status = parseArguments(argc, argv, params);
	if (status) {
		return status - 1;  // Displaying help returns 0
	}
	status = checkParameters(params);
	if (status) {
		return status;
	}

    if (params.verbose) {
		std::cout << "--- This is " << argv[0] << ", compiled on " << __DATE__
			<< " at " << __TIME__ << " ---" << std::endl;
        printParameters(params);
		std::cout << std::endl;
    }

	status = runSimulations(params);

	return status;
}
