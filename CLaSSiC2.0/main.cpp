#include <iostream>
#include "simulation.hpp"

int main() {
	/*
	Simulation is initialized here and called with paramters such as temperature.
	TODO: Implement a way to run different structures and parameters from shell command
	*/
	Simulation sim;

	for (int i = 0; i < constants::nSims; i++) {
		// sim.initialize((0.0001+(float)i) / constants::nSims * 0.5 * constants::pi);
		// sim.run(i, ((float)i) / constants::nSims*100);
		// sim.initialize(0.872665);
		sim.run(i, 5);
	}
	
	return 0;
}