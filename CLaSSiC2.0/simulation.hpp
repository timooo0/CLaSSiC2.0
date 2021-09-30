#include <vector>
#include "integrator.hpp"

class Simulation {
private:
	double spin[constants::nAtoms*3];
	double position[constants::nAtoms * 3];
	double randomField[constants::nAtoms * 3];
	std::vector<std::vector<double>> neighbours = std::vector<std::vector<double>>(constants::nAtoms, std::vector<double>());
	Integrator integrator;

public:
	Simulation();
	void load();
	void initialize();
	void normalize();
	void run();
};