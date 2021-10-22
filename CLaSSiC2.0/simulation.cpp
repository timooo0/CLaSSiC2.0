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
		- Removes all previous data
		- Calls load to build the crystal
		- Calls initialize to set to start orientation of the spins
	*/

	int i = 0;
	int succes = 0;
	std::filesystem::remove_all("data/");
	std::filesystem::create_directory("data");

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

		TODO: remove the file read since crystals are now constructed from scratch
	*/

	std::fstream file;
	std::string line;
	size_t pos = 0;
	int i = 0;
	file.open(constants::inputFile);
	if (!file.is_open())
	{
		std::cout << "Could not open file.\n";
	}

	// Get the positions
	while (std::getline(file, line) && i / 3 < constants::baseAtoms)
	{
		while ((pos = line.find(", ")) != std::string::npos)
		{
			position.push_back(std::stod(line.substr(0, pos)));
			line.erase(0, pos + 1);
			i++;
		}
		if (!line.empty())
		{
			position.push_back(std::stod(line));
			i++;
		}
	}

	file.close();

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
#if printPositions
	std::cout << "Atom positions: \n";
	for (int i = 0; i < position.size() / 3; i++)
	{
		std::cout << "x: " << position[3 * i] << " y: " << position[3 * i + 1] << " z: " << position[3 * i + 2] << std::endl;
	}
#endif

#if printNeighbours
	std::cout << "Neighbours: \n";
	for (int i = 0; i < constants::nAtoms; i++)
	{
		for (int j = 0; j < neighbours[i].size(); j++)
		{
			std::cout << neighbours[i][j] << " ";
		}
		std::cout << std::endl;
	}
#endif
}

void Simulation::initialize(double angle)
{
	/*
	Initializing the spin orientations
	*/
	float mode = 4;
	switch (constants::spinInit)
	{
	case 0:
		//Alligned with the z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0;
			spin[3 * i + 1] = 0;
			spin[3 * i + 2] = 1;
		}
		break;
	case 1:
		//Angle of 45 degrees with z-axis
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = std::cos((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(angle);
			spin[3 * i + 1] = std::sin((float)i / constants::nAtoms * 2 * constants::pi) * std::sin(angle);
			spin[3 * i + 2] = std::cos(angle);
		}
		break;
	case 2:
		//Small z angle for spin waves
		for (int i = 0; i < constants::nAtoms; i++)
		{
			spin[3 * i] = 0.001 * std::cos((float)i / constants::nAtoms * mode * 2 * constants::pi) * std::sin(angle);
			spin[3 * i + 1] = 0.001 * std::sin((float)i / constants::nAtoms * mode * 2 * constants::pi) * std::sin(angle);
			spin[3 * i + 2] = 1;
		}
		normalize();
		break;
	case 3:
		//Rotor mode for two spins
		spin[0] = std::cos(0) * std::sin(angle);
		spin[1] = 0;
		spin[2] = std::cos(angle);
		spin[3] = -spin[0];
		spin[4] = 0;
		spin[5] = spin[2];
		break;
	}

//Printing
# if printInitialize
	std::cout << "Spin initialization: \n";
	for (int i = 0; i < constants::nAtoms; i++) {
		std::cout << "x: " << spin[3 * i] << " y: " << spin[3 * i + 1] << " z: " << spin[3 * i + 2] << std::endl;
	}
# endif
}

void Simulation::run(int iterator, double temperature)
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
	std::ofstream file;
	std::string dataFile = constants::outputFile;
	dataFile.insert(dataFile.find("."), std::to_string(iterator));
	file.open(dataFile, std::ios::binary);
	writeConstants(file, temperature);

	// Temperature fluctatuations
	std::random_device seed;
	std::default_random_engine engine{seed()};
	std::normal_distribution<double> temperatureDistribution;
	if (constants::temperatureSigma != 0 && temperature != 0)
	{
		temperatureDistribution = std::normal_distribution<double>(0, std::sqrt(constants::temperatureSigma * temperature));
	}

	// Stopwatch
	auto start = std::chrono::high_resolution_clock::now();

	// Actual run
	for (int i = 0; i < constants::steps; i++)
	{
		if (constants::temperatureSigma != 0 && temperature != 0)
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
	std::cout << "\nDuration: " << (float)msInt.count() / 1000 << "seconds" << std::endl;
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

void Simulation::writeConstants(std::ofstream &f, double &temperature)
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
	f.write((char *)&temperature, sizeof(double));
	f.write((char *)&constants::length, sizeof(double));
}


