#include <string>
#include <vector>
#include <cmath>

namespace constants {
	//Physical constants
	const double boltzmann = 1.38064852e-23;
	const double gFactor = -2.002;
	const double bohrMagneton = 9.274009994e-24;
	const double hBar = 1.054571800e-34;
	const double gamma = -1.760859644e11;
	const double pi = 3.14159265359;

	//System constants
	const double spinSize = 3.5;
	const int baseAtoms = 1;
	const int nUnitCells = 1;
	const int nDimensions = 1;
	const int nAtoms = baseAtoms * std::pow(nUnitCells, nDimensions);
	const std::vector<std::vector<double>> unitVectors = {{baseAtoms, 0, 0}, {0, 1, 0}, {0, 0, 0}};

	//Simulation settings;
	const int spinInit = 1; //0: Z-direction, 1: angle/45 degrees
	const double minDistance = 1.1;
	const int nSims = 1;

	//Simulation parameters
	const double dt = 1e-15;
	const int steps = 1e6;
	const double anisotropy[3] = { 0, 0, 1 }; // Should be normalized
	const double anisotropyStrength = 0;
	const double magneticField[3] = { 0, 0, 20 };
	const double J = 2*boltzmann;
	const double lambda = 0;

	//File constants
	const std::string inputFile = "spin.csv";
	const std::string outputFile = "data//data.dat";
	const double offset = 12; //offset, dt, J, lambda, 3 magneticField, 3 anisotropy, temperature, length
	const double length = 3 * nAtoms;
	
	//Prefactors
	const double exchangePrefactor = -2 * J * spinSize / (gFactor * bohrMagneton);
	const double temperatureSigma =2. * lambda * boltzmann * hBar * dt / (gFactor * gFactor * bohrMagneton * bohrMagneton*spinSize);
	const double anisotropyMatrix[9] = { 
		anisotropyStrength * anisotropy[0] * anisotropy[0],
		anisotropyStrength * anisotropy[0] * anisotropy[1],
		anisotropyStrength * anisotropy[0] * anisotropy[2],
		anisotropyStrength * anisotropy[1] * anisotropy[0],
		anisotropyStrength * anisotropy[1] * anisotropy[1],
		anisotropyStrength * anisotropy[1] * anisotropy[2],
		anisotropyStrength * anisotropy[2] * anisotropy[0],
		anisotropyStrength * anisotropy[2] * anisotropy[1],
		anisotropyStrength * anisotropy[2] * anisotropy[2] };

}