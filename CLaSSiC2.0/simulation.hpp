#include <vector>
#include "integrator.hpp"

class Simulation
{
private:
	std::vector<double> spin = std::vector<double>(constants::nAtoms * 3);
	std::vector<double> position;
	std::vector<double> randomField = std::vector<double>(constants::nAtoms * 3, 0);
	std::vector<std::vector<double>> neighbours = std::vector<std::vector<double>>(constants::nAtoms, std::vector<double>());
	Integrator integrator;

public:
	Simulation();
	void load();
	void initialize(double angle = 0.25 * constants::pi);
	void run(int iterator = 0, double temperature = 0);
	void normalize();
	void writeConstants(std::ofstream &file, double &temperature);
};