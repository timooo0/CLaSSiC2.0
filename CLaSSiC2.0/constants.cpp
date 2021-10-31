#include <string>
#include <vector>
#include <cmath>
#include "constants.hpp"

namespace constants
{
	// Physical constants
    double boltzmann = 1.38064852e-23;
    double gFactor = -2.002;
    double bohrMagneton = 9.274009994e-24;
    double hBar = 1.054571817e-34;
    double gamma = -1.760859644e11;
    double pi = 3.14159265359;

	// System constants
	double spinSize = 3.5;
	int baseAtoms = 1;
	int nUnitCells = 1;
	int nDimensions = 1;
	int nAtoms = baseAtoms * std::pow(nUnitCells, nDimensions);
	std::vector<std::vector<double>> unitVectors = {{(double)baseAtoms, 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};

	// Simulation settings;
	int spinInit = 1; // 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
	double angle = 0.25 * constants::pi;
	double minDistance = 1.1;
	int nSims = 1;

	// Simulation parameters
	double dt = 1e-14;
	int steps = 1e6;
	double anisotropy[3] = {0, 0, 1}; // Should be normalized
	double anisotropyStrength = 0;
	double magneticField[3] = {0, 0, 50};
	double J = 0 * boltzmann;
	double temperature = 0;
	double lambda = 1e-2;

	// File constants
	std::string inputFile = "spin.csv";
	std::string outputFile = "data/data.dat";
	double offset = 14; // offset, nAtoms, dt, steps. J, lambda, 3 magneticField, 3 anisotropy, temperature, length
	double length = 3 * nAtoms;

	// Prefactors
	double exchangePrefactor = -2 * J * spinSize / (gFactor * bohrMagneton);
	double temperatureSigma = 2. * lambda * boltzmann * hBar * dt / (gFactor * gFactor * bohrMagneton * bohrMagneton * spinSize);
	double anisotropyMatrix[9] = {
		anisotropyStrength * anisotropy[0] * anisotropy[0],
		anisotropyStrength * anisotropy[0] * anisotropy[1],
		anisotropyStrength * anisotropy[0] * anisotropy[2],
		anisotropyStrength * anisotropy[1] * anisotropy[0],
		anisotropyStrength * anisotropy[1] * anisotropy[1],
		anisotropyStrength * anisotropy[1] * anisotropy[2],
		anisotropyStrength * anisotropy[2] * anisotropy[0],
		anisotropyStrength * anisotropy[2] * anisotropy[1],
		anisotropyStrength * anisotropy[2] * anisotropy[2]};

}

// Enables printing for debug purposes
bool printPositions = false;
bool printNeighbours = false;
bool printInitialize = false;
