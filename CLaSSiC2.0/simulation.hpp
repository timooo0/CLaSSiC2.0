#include <vector>
#include "integrator.hpp"

class Simulation
{
private:
	std::vector<double> spin = std::vector<double>(constants::nAtoms * 3);
	std::vector<double> position;
	std::vector<double> randomField = std::vector<double>(constants::nAtoms * 3, 0);
	std::vector<std::vector<int>> neighbours = std::vector<std::vector<int>>(constants::nAtoms, std::vector<int>());
	Integrator integrator;

public:
	Simulation();
	void load();
	void initialize();
	void run();
	void normalize();
	void writeConstants(std::ofstream &file);

	std::vector<double>* getSpin();
	std::vector<double>* getPosition();
	std::vector<double>* getRandomField();
	std::vector<std::vector<int>>* getNeighbours();
};