#include <string>
#include <vector>
#include <iostream>
#include <cmath>
#include "constants.hpp"

namespace constants
{
	// Physical constants
    double boltzmann = 1.38064852e-23;
    double gFactor = 2.002;
    double bohrMagneton = 9.274009994e-24;
    double hBar = 1.054571817e-34;
    double gamma = -1.760859644e11;
    double pi = 3.14159265359;

	// System constants
	float spinSize = 3.5;
	int nUnitCells = 8;
	int nDimensions = 1;
	int nAtoms;
	int nNeighbours;
	std::vector<std::vector<float>> unitVectors = {{1., 0., 0.}};
	std::vector<std::vector<float>> basisPosition = {{0., 0., 0.}};

	// Simulation settings;
	int spinInit = 2; // 0: Z-direction, 1: angle with z-axis, 2: small angle with z-axis, 3: rotor, 4: Z-direction alternating
	int geometry = 1;
	float angle = 0.25 * constants::pi;
	int mode = 1;
	float minDistance = 1.1;
	bool periodicBoundary = true;
	bool stabilize = false;

	// Simulation parameters
	float dt = 1e-15;
	int steps = (int) 1e6;
	int burnInSteps = 0;
	float anisotropy[3] = {0, 0, 1}; // Should be normalized
	float anisotropyStrength = 0;
	float anisotropyAxis = 0;
	float anisotropyPlane = 0;
	float magneticField[3] = {0, 0, 0};
	float J = -2 * boltzmann;
	float temperature = 0;
	float lambda = 0;

	// File constants
	std::string positionFile = "data/position.csv";
	std::string outputFile = "data/data.dat";
	std::string energyFile = "data/energy.dat";
	
	 // offset, nAtoms, dt, steps. J, lambda, 3 magneticField, 3 anisotropy, temperature, length, geometry, nUnitCells, anisotropyStrength
	float offset = 17.0f;
	float length;

	// Prefactors
	float exchangePrefactor;
	float temperatureSigma;
	float anisotropyMatrix[9];

}

// Enables printing for debug purposes
bool printPositions = false;
bool printNeighbours = false;
bool printInitialize = false;

void setDerivatives(){
	constants::nAtoms = constants::basisPosition.size() * (int)std::pow(constants::nUnitCells, constants::nDimensions);
	std::cout << "nAtoms: " << constants::nAtoms << std::endl;
	// constants::unitVectors = {{(float)constants::baseAtoms, 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};

	constants::length = 3 * constants::nAtoms;
	float sub = (constants::gFactor * constants::gFactor * constants::bohrMagneton * constants::bohrMagneton * constants::spinSize);
	std::cout << "yo: " << 2.0f * constants::lambda * constants::boltzmann * constants::hBar * constants::dt * constants::temperature << ", " << sub << std::endl;


	constants::exchangePrefactor = 2 * constants::J * constants::spinSize / (constants::gFactor * constants::bohrMagneton);
	constants::temperatureSigma = sqrtf(2.0f * constants::lambda * constants::boltzmann * constants::hBar * constants::dt * constants::temperature/ 
		(constants::gFactor * constants::gFactor * constants::bohrMagneton * constants::bohrMagneton * constants::spinSize));
	constants::anisotropyStrength = constants::anisotropyAxis;
	float anisM[9] = {
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