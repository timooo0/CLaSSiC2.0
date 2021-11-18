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
	int nUnitCells = 12;
	int nDimensions = 1;
	int nAtoms;
	std::vector<std::vector<double>> unitVectors;

	// Simulation settings;
	int spinInit = 2; // 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
	double angle = 0.25 * constants::pi;
	int mode = 0;
	double minDistance = 1.1;

	// Simulation parameters
	double dt = 1e-15;
	int steps = 1e6;
	double anisotropy[3] = {0, 0, 1}; // Should be normalized
	double anisotropyStrength = 0;
	double anisotropyAxis = 0;
	double anisotropyPlane = 0;
	double magneticField[3] = {0, 0, 20};
	double J = 2 * boltzmann;
	double temperature = 1;
	double lambda = 1e-3;

	// File constants
	std::string inputFile = "spin.csv";
	std::string outputFile = "data/data.dat";
	double offset = 15; // offset, nAtoms, dt, steps. J, lambda, 3 magneticField, 3 anisotropy, temperature, length, anisotropyStrength
	double length;

	// Prefactors
	double exchangePrefactor;
	double temperatureSigma;
	double anisotropyMatrix[9];

}

// Enables printing for debug purposes
bool printPositions = false;
bool printNeighbours = false;
bool printInitialize = false;

void setDerivatives(){
	constants::nAtoms = constants::baseAtoms * std::pow(constants::nUnitCells, constants::nDimensions);
	constants::unitVectors = {{(double)constants::baseAtoms, 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};

	constants::length = 3 * constants::nAtoms;

	constants::exchangePrefactor = -2 * constants::J * constants::spinSize / (constants::gFactor * constants::bohrMagneton);
	constants::temperatureSigma = 2. * constants::lambda * constants::boltzmann * constants::hBar * constants::dt / 
		(constants::gFactor * constants::gFactor * constants::bohrMagneton * constants::bohrMagneton * constants::spinSize);
	double anisM[9] = {
		constants::anisotropyStrength * constants::anisotropy[0] * constants::anisotropy[0],
		constants::anisotropyStrength * constants::anisotropy[0] * constants::anisotropy[1],
		constants::anisotropyStrength * constants::anisotropy[0] * constants::anisotropy[2],
		constants::anisotropyStrength * constants::anisotropy[1] * constants::anisotropy[0],
		constants::anisotropyStrength * constants::anisotropy[1] * constants::anisotropy[1],
		constants::anisotropyStrength * constants::anisotropy[1] * constants::anisotropy[2],
		constants::anisotropyStrength * constants::anisotropy[2] * constants::anisotropy[0],
		constants::anisotropyStrength * constants::anisotropy[2] * constants::anisotropy[1],
		constants::anisotropyStrength * constants::anisotropy[2] * constants::anisotropy[2]};
	memcpy(constants::anisotropyMatrix, anisM, sizeof(anisM));
}