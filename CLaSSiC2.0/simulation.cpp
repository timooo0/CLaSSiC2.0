#include "simulation.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <algorithm>
#include <filesystem>
#include <chrono>

Simulation::Simulation()
{
	/*
	Constructer for the simulation class.
		- Calls load to build the crystal
		- Calls initialize to set to start orientation of the spins
	*/

	load();
	initialize();
}

void Simulation::load()
{
	/*
	Builds the crystal from the ground up.
		1. Makes an unit cell
		2. Extends the unit cell in specified dimensions
		3. Determine nearest neighbours
		4. Apply periodic boundary conditions
	*/

	// Create the basis atom
	position.push_back(0); // x
	position.push_back(0); // y
	position.push_back(0); // z

	// Add more unit cells
	std::vector<double> basePosition;
	for (int dim = 0; dim < constants::nDimensions; dim++)
	{
		basePosition = position;
		for (int cell = 1; cell < constants::nUnitCells; cell++)
		{
			for (int i = 0; i < basePosition.size() / 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					position.push_back(basePosition[3 * i + j] + cell * (constants::unitVectors[dim][j]));
				}
			}
		}
	}

	// Calculate nearest neighbours
	double distance;
	for (int i = 0; i < constants::nAtoms; i++)
	{
		for (int j = 0; j < constants::nAtoms; j++)
		{
			if (i != j)
			{
				distance = std::sqrt((position[3 * i] - position[3 * j]) * (position[3 * i] - position[3 * j]) +
									 (position[3 * i + 1] - position[3 * j + 1]) * (position[3 * i + 1] - position[3 * j + 1]) +
									 (position[3 * i + 2] - position[3 * j + 2]) * (position[3 * i + 2] - position[3 * j + 2]));
				if (distance < constants::minDistance)
				{
					neighbours[i].push_back(j);
				}
			}
		}
	}

	// Periodic boundary conditions
	bool applyBoundry = true;
	for (int dim = 0; dim < constants::nDimensions; dim++)
	{
		for (int i = 0; i < constants::nAtoms; i++)
		{
			for (int j = 0; j < constants::nAtoms; j++)
			{
				if (i != j)
				{
					applyBoundry = true;
					for (int k = 0; k < 3; k++)
					{
						if (constants::unitVectors[dim][k] != 0)
						{
							if (position[3 * i + k] + constants::nUnitCells * constants::unitVectors[dim][k] - 1 != position[3 * j + k])
							{
								applyBoundry = false;
							}
						}
						else
						{
							if (position[3 * i + k] != position[3 * j + k])
							{
								applyBoundry = false;
							}
						}
					}

					if (applyBoundry && std::find(neighbours[i].begin(), neighbours[i].end(), j) == neighbours[i].end())
					{
						neighbours[i].push_back(j);
						neighbours[j].push_back(i);
					}
				}
			}
		}
	}

	// Printing

	std::cout << "Atom positions: \n";
	for (int i = 0; i < position.size() / 3; i++)
	{
		std::cout << "x: " << position[3 * i] << " y: " << position[3 * i + 1] << " z: " << position[3 * i + 2] << std::endl;
	}


	std::cout << "Neighbours: \n";
	for (int i = 0; i < constants::nAtoms; i++)
	{
		for (int j = 0; j < neighbours[i].size(); j++)
		{
			std::cout << neighbours[i][j] << " ";
		}
		std::cout << std::endl;
	}

}

void Simulation::initialize()
{
	/*
	Initializing the spin orientations
	*/

	int swap = 1;
	double delta = 1;
	switch (constants::spinInit)
	{
	case 0:
		// Alligned with the z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0;
			spin[3 * i + 1] = 0;
			spin[3 * i + 2] = 1;
		}
		break;
	case 1:
		// Angle with z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = std::cos((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			spin[3 * i + 1] = std::sin((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(constants::angle);
			spin[3 * i + 2] = std::cos(constants::angle);
		}
		break;
	case 2:
		// Small z angle for spin waves
		for (int i = 0; i < (int)constants::nAtoms; i++)
		{
			if (constants::J < 0)
				{
					if (i % 2 == 0)
					{
						swap = 1;
					}
					else
					{
						swap = -1;
					}
				}
				spin[3 * i] = swap * 0.001 * std::cos((float)i /constants::nAtoms * constants::mode * 2 * constants::pi);
				spin[3 * i + 1] = swap * 0.001 * std::sin((float)i /constants::nAtoms * constants::mode * 2 * constants::pi);
				spin[3 * i + 2] = swap * 1;
		}
		normalize();
		break;
	case 3:
		if (constants::J < 0)
		{
			swap = -1;
		}
		// Rotor mode for two spins
		spin[0] = std::cos(constants::angle);
		spin[1] = std::sin(constants::angle);
		spin[2] = std::cos(constants::angle);
		spin[3] = swap * std::cos(constants::angle);
		spin[4] = std::sin(constants::angle); 
		spin[5] = swap * std::cos(constants::angle);
		break;
	case 4:
		// 2D spin waves
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0;
			spin[3 * i + 1] = 0;
			if ((i + i / constants::nUnitCells) % 2 == 0)
			{
				spin[3 * i + 2] = 1;
			}
			else
			{
				spin[3 * i + 2] = -1;
			}
		}
	case 5:
		// Triangle
		std::cout << "Triangle:\n";
		spin[0] = delta*std::cos(30./180.*constants::pi);
		spin[1] = delta*std::sin(30./180.*constants::pi);
		spin[2] = 1;

		spin[3] = -delta*std::cos(30./180.*constants::pi);
		spin[4] = delta*std::sin(30./180.*constants::pi);
		spin[5] = 1;
		
		spin[6] = delta*0;
		spin[7] = -delta*1;
		spin[8] = -1;
	break;
	}
	// Printing

	std::cout << "Spin initialization: \n";
	for (int i = 0; i < constants::nAtoms; i++)
	{
		std::cout << "Atom: " << i / constants::nUnitCells << ", " << i % constants::nUnitCells << "|x: " << spin[3 * i] << " y: " << spin[3 * i + 1] << " z: " << spin[3 * i + 2] << std::endl;
	}
}

void Simulation::run()
{
	/*
	Actually runs the simulation
		1. Write constants to file
		2. Determine temperature fluctuations
		3. Stopwatch to track time
		4. Actual run
		5. close file and print time
	*/

	// Open and write to file
	int i = 0;
	std::ofstream file;
	std::string dataFile = constants::outputFile;
	dataFile.insert(dataFile.find("."), std::to_string(i));
	while (std::filesystem::exists(dataFile))
	{
		i++;
		dataFile = constants::outputFile;
		dataFile.insert(dataFile.find("."), std::to_string(i));
	}
	file.open(dataFile, std::ios::binary);
	writeConstants(file);

	// Temperature fluctatuations
	std::random_device seed;
	std::default_random_engine engine{seed()};
	std::normal_distribution<double> temperatureDistribution;
	if (constants::temperatureSigma != 0 && constants::temperature != 0)
	{
		temperatureDistribution = std::normal_distribution<double>(0, std::sqrt(constants::temperatureSigma * constants::temperature));
	}

	// Stopwatch
	auto start = std::chrono::high_resolution_clock::now();

	// Actual run
	for (int i = 0; i < constants::steps; i++)
	{
		if (constants::temperatureSigma != 0 && constants::temperature != 0)
		{
			for (int i = 0; i < constants::nAtoms * 3; i++)
			{
				randomField[i] = temperatureDistribution(engine);
			}
		}

		integrator.integrate(neighbours, spin, randomField);
		normalize();

		if (i % 1000 == 0)
		{
			std::cout << "Progress: " << (float)i / constants::steps * 100 << "%\r";
		}
		if (i % 100 == 0)
		{
			file.write((char *)&spin[0], sizeof(double) * constants::nAtoms * 3);
		}
	}

	// Wrap up
	file.close();
	auto msInt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
	std::cout << "\nDuration: " << (float)msInt.count() / 1000 << "seconds\n"
			  << std::endl;
}

void Simulation::normalize()
{
	/*
	Helper function to normalize all of the spins
	*/
	double invLength;
	for (int i = 0; i < constants::nAtoms; i++)
	{
		invLength = 1 / std::sqrt(spin[3 * i] * spin[3 * i] + spin[3 * i + 1] * spin[3 * i + 1] + spin[3 * i + 2] * spin[3 * i + 2]);
		spin[3 * i] *= invLength;
		spin[3 * i + 1] *= invLength;
		spin[3 * i + 2] *= invLength;
	}
}

void Simulation::writeConstants(std::ofstream &f)
{
	/*
	Writes all constants and system parameters to the start of the file
	*/
	double saveAtoms = constants::nAtoms;
	double saveSteps = constants::steps;
	f.write((char *)&constants::offset, sizeof(double));
	f.write((char *)&saveAtoms, sizeof(double));
	f.write((char *)&constants::dt, sizeof(double));
	f.write((char *)&saveSteps, sizeof(double));
	f.write((char *)&constants::J, sizeof(double));
	f.write((char *)&constants::lambda, sizeof(double));
	f.write((char *)&constants::magneticField[0], 3 * sizeof(double));
	f.write((char *)&constants::anisotropy[0], 3 * sizeof(double));
	f.write((char *)&constants::temperature, sizeof(double));
	f.write((char *)&constants::length, sizeof(double));
	f.write((char *)&constants::anisotropyStrength, sizeof(double));
}
