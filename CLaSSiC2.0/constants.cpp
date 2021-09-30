#include <string>



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
	const int nAtoms = 1;

	//Simulation settings;
	const int spinInit = 1; //0: Z-direction, 1: 45 degrees
	const double minDistance = 1.1;

	//Simulation parameters
	const double dt = 1e-15;
	const int steps = 1e6;
	const double anisotropy[3] = { 0, 0, 0 };
	const double magneticField[3] = { 0, 0, 5 };
	const double J = 0;
	const double lambda = 0;

	//File constants
	const std::string inputFile = "spin.csv";
	const std::string outputFile = "data.dat";
	const double offset = 12; //offset, dt, J, lambda, 3 magneticField, 3 anisotropy, temperature, length
	const double length = 3 * nAtoms;
	
	//Prefactors
	const double exchangePrefactor = -2 * J * spinSize / (gFactor * bohrMagneton);
}