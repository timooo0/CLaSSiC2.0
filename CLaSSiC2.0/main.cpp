#include <iostream>
#include "simulation.hpp"


int main() {
	Simulation sim;

	for (int i = 0; i < constants::nSims; i++) {
		//sim.initialize((0.0001+(float)i) / constants::nSims * 0.5 * constants::pi);
		sim.run(i, ((float)i) / constants::nSims*100);
		// sim.run(i, 0);
		sim.initialize();
	}
	
	return 0;
}