#include <iostream>
#include <string>
#include "simulation.hpp"
#include "constants.hpp"

int main(int argc, char *argv[])
{
	/*
	Get all the parameters from the program call
	Simulation is initialized here and run.
	*/
	for (int i = 1; i < argc; i++)
	{
		if (std::string(argv[i]) == "-dt"){
			constants::dt = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-steps"){
			constants::steps = (int)std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-J"){
			constants::J = std::stod(argv[i+1])*constants::boltzmann;
		}
		else if (std::string(argv[i]) == "-lambda"){
			constants::lambda = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-B"){
			constants::magneticField[2] = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-anisotropy"){
			constants::anisotropyStrength = std::stod(argv[i+1]);
			std::cout << "anisotropyStrength: " << constants::anisotropyStrength << std::endl;
		}
		else if (std::string(argv[i]) == "-T"){
			constants::temperature = std::stod(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-init"){
			constants::spinInit = std::stoi(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-angle"){
			constants::angle = std::stod(argv[i+1])*constants::pi/180.;
			std::cout << "angle: " << constants::angle << std::endl;
		}
		else if (std::string(argv[i]) == "-mode"){
			constants::mode = std::stoi(argv[i+1]);
		}
		else if (std::string(argv[i]) == "-structure"){
			if (std::string(argv[i+1]) == "single"){
				constants::nDimensions = 1;
			}
			if (std::string(argv[i+1]) == "line"){
				constants::nDimensions = 1;
				constants::unitVectors = {{1., 0., 0.}, {0., 0., 0.}, {0., 0., 0.}};
			}
			else if (std::string(argv[i+1]) == "square") {
				constants::nDimensions = 2;
				constants::unitVectors = {{1., 0., 0.}, {0., 1., 0.}, {0., 0., 0.}};
			} else {
				std::cerr << "-structure only has the following options: single and square\n";
			}
		}
		else if (std::string(argv[i]) == "-nCellsX"){
			constants::nUnitCells = (int)std::stod(argv[i+1]);
		}

		i++;
	}
	setDerivatives();
	Simulation sim;
	sim.run();
	
	return 0;
}