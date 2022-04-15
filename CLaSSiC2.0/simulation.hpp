#include <vector>
#include "integrator.hpp"

class Simulation
{
private:
	Vector3D spin =  Vector3D(constants::nAtoms, true);
	std::vector<float> position;
	std::vector<float> randomField = std::vector<float>(constants::nAtoms * 3, 0);
	std::vector<std::vector<int>> neighbours = std::vector<std::vector<int>>(constants::nAtoms, std::vector<int>());
	std::vector<float> totalEnergy = std::vector<float>((int)constants::steps/100);
	Integrator integrator;

public:
	Simulation();
	void load();
	void initialize();
	void run();
	void normalize();
	std::string addFileNumber(std::string input);
	void writeConstants(std::ofstream &file);
	void writePositions();
	float distance(std::vector<float>::iterator a, std::vector<float>::iterator b);
	void addNeighbours(std::vector<float> a, std::vector<float> b, int i, int j);

	Vector3D* getSpin();
	std::vector<float>* getPosition();
	std::vector<float>* getRandomField();
	std::vector<std::vector<int>>* getNeighbours();
};