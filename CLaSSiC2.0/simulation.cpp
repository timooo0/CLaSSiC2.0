#include "simulation.hpp"
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <chrono>

Simulation::Simulation()
{
	int i = 0;
	int succes = 0;
	while (succes == 0) {
		std::string dataFile = constants::outputFile;
		dataFile.insert(dataFile.find("."), std::to_string(i));
		succes = std::remove(constants::outputFile.c_str());
		i++;
	}

	load();
	initialize();
}

void Simulation::load()
{
	std::fstream file;
	std::string line;
	size_t pos = 0;
	int i = 0;
	file.open(constants::inputFile);

	//Get the positions
	while (std::getline(file, line)&& i/3<constants::nAtoms) {
		while ((pos = line.find(', ')) != std::string::npos) {
			position[i] = std::stod(line.substr(0, pos));
			line.erase(0, pos + 1);
			i++;
		}
		if (!line.empty()) {
			position[i] = std::stod(line);
			i++;
		}
	}
	
	file.close();

	std::cout << "Atom positions: \n";
	for (int i = 0; i < constants::nAtoms; i++) {
		std::cout << "x: " << position[3 * i] << " y: " << position[3 * i + 1] << " z: " << position[3 * i + 2] << std::endl;
	}
	std::cout << std::endl;

	//Calculate nearest neighbours
	double distance;
	for (int i = 0; i < constants::nAtoms; i++) {
		for (int j = 0; j < constants::nAtoms; j++) {
			if (i != j) {
				distance = std::sqrt((position[3*i] - position[3 * j]) * (position[3 * i] - position[3 * j]) +
					(position[3 * i+1] - position[3 * j+1])*(position[3 * i+1] - position[3 * j+1]) +
					(position[3 * i+2] - position[3 * j+2])*(position[3 * i+2] - position[3 * j+2]));
				if (distance < constants::minDistance) {
					neighbours[i].push_back(j);
				}
			}
		}
	}
	
	//Periodic boundary conditions
	bool applyBoundry = true;
	for (int i = 0; i < constants::nAtoms; i++) {
		for (int j = 0; j < constants::nAtoms; j++) {
			if (i != j) {
				applyBoundry = true;
				for (int k = 0; k < 3; k++) {
					for (int l = 0; l<constants::nDimensions;  l++) {
						if (position[3*i+k]+ constants::unitVectors[l][k] != position[3*j+k]) {
							applyBoundry = false;
						}
					}
				}
				if (applyBoundry) {
					neighbours[i].push_back(j);
					neighbours[j].push_back(i);
				}
			}
		}
	}
	std::cout << "Neighbours: \n";
	for (int i = 0; i < constants::nAtoms; i++) {
		for (int j = 0; j < neighbours[i].size(); j++) {
			std::cout << neighbours[i][j] << " ";
		}
		std::cout << std::endl;
	}
}

void Simulation::initialize(double angle)
{
	switch (constants::spinInit) {
		case 0:
			for (int i = 0; i < constants::nAtoms; i++) {
				spin[3 * i] = 0;
				spin[3 * i + 1] = 0;
				spin[3 * i + 2] = 1;

			}
			break;
		case 1:
			for (int i = 0; i < constants::nAtoms; i++) {
				spin[3 * i]     = 0.001*std::cos((float)i / constants::nAtoms * 2 * constants::pi)*std::sin(angle);
				spin[3 * i + 1] = 0.001*std::sin((float)i / constants::nAtoms * 2 * constants::pi)*std::sin(angle);
				spin[3 * i + 2] = std::cos(angle);
			}
			break;
	}

	//std::cout << "Spin initialization: \n";
	//for (int i = 0; i < constants::nAtoms; i++) {
	//	std::cout << "x: " << spin[3 * i] << " y: " << spin[3 * i + 1] << " z: " << spin[3 * i + 2] << std::endl;
	//}
}

void Simulation::normalize()
{
	double invLength;
	for (int i = 0; i < constants::nAtoms; i++) {
		invLength = 1 / std::sqrt(spin[3 * i] * spin[3 * i] + spin[3 * i + 1] * spin[3 * i + 1] + spin[3 * i + 2] * spin[3 * i + 2]);
		spin[3 * i] *= invLength;
		spin[3 * i + 1] *= invLength;
		spin[3 * i + 2] *= invLength;
	}
}

void Simulation::run(int iterator, double temperature)
{
	//Function timing
	auto start = std::chrono::high_resolution_clock::now();
	
	//Temperature fluctatuations
	std::random_device seed;
	std::default_random_engine engine{ seed() };
	std::normal_distribution<double> temperatureDistribution;
	if (constants::temperatureSigma != 0 && temperature != 0) {
		temperatureDistribution = std::normal_distribution<double>(0, std::sqrt(constants::temperatureSigma*temperature));
	}
	for (int i = 0; i < constants::nAtoms*3; i++) {
		randomField[0] = 0;
	}


	std::ofstream file;
	std::string dataFile = constants::outputFile;
	dataFile.insert(dataFile.find("."),std::to_string(iterator));
	file.open(dataFile, std::ios::binary);

	file.write((char*)&constants::offset, sizeof(double));
	file.write((char*)&constants::dt, sizeof(double));
	file.write((char*)&constants::J, sizeof(double));
	file.write((char*)&constants::lambda, sizeof(double));
	file.write((char*)&constants::magneticField[0], 3*sizeof(double));
	file.write((char*)&constants::anisotropy[0], 3*sizeof(double));
	file.write((char*)&temperature, sizeof(double));
	file.write((char*)&constants::length, sizeof(double));

	for (int i = 0; i < constants::steps; i++) {
		if (constants::temperatureSigma != 0 && temperature != 0) {
			for (int i = 0; i < constants::nAtoms * 3; i++) {
				randomField[i] = temperatureDistribution(engine);
			}
		}

		integrator.integrate(neighbours, spin, randomField);
		normalize();
		//std::cout << "x: " << spin[0] << " y: " << spin[1] << " z: " << spin[2] << std::endl;
		file.write((char*)&spin[0], sizeof(double)*constants::nAtoms*3);
	}
	file.close();

	auto msInt = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
	std::cout << "Duration: " << (float)msInt.count() / 1000 << "seconds" << std::endl;
}
