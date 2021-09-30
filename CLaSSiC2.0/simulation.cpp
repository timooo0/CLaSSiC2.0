#include "simulation.hpp"
#include <fstream>
#include <iostream>
#include <cmath>

Simulation::Simulation()
{
	std::remove(constants::outputFile.c_str());
	load();
	initialize();
}

void Simulation::load()
{
	std::fstream file;
	std::string line;
	size_t pos = 0;
	unsigned int i = 0;
	file.open(constants::inputFile);

	while (std::getline(file, line)) {
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

	double distance;
	for (int i = 0; i < constants::nAtoms; i++) {
		for (int j = 0; j < constants::nAtoms; j++) {
			if (i != j) {
				distance = std::sqrt(position[i] * position[j] + position[i + 1] * position[j + 1] + position[i + 2] * position[j + 2]);
				if (distance < constants::minDistance) {
					neighbours[i].push_back(j);
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

void Simulation::initialize()
{
	switch (constants::spinInit) {
		case 0:
			for (int i = 0; i < constants::nAtoms*3; i++) {
				spin[i] = 0;
			}
			break;
		case 1:
			for (int i = 0; i < constants::nAtoms; i++) {
				spin[3 * i] = std::sin(0.25*constants::pi);
				spin[3 * i + 1] = 0;
				spin[3*i+2] = std::cos(0.25*constants::pi);
			}
			break;
	}

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

void Simulation::run()
{
	//temperature effects
	double temperature = 0;
	for (int i = 0; i < constants::nAtoms*3; i++) {
		randomField[0] = 0;
	}


	std::ofstream file;
	file.open(constants::outputFile, std::ios::binary);
	file.write((char*)&constants::offset, sizeof(double));
	file.write((char*)&constants::dt, sizeof(double));
	file.write((char*)&constants::J, sizeof(double));
	file.write((char*)&constants::lambda, sizeof(double));
	file.write((char*)&constants::magneticField[0], 3*sizeof(double));
	file.write((char*)&constants::anisotropy[0], 3*sizeof(double));
	file.write((char*)&temperature, sizeof(double));
	file.write((char*)&constants::length, sizeof(double));

	for (int i = 0; i < constants::steps; i++) {
		file.write((char*)&spin[0], sizeof(spin));
		integrator.integrate(neighbours, spin, randomField);
	}
	file.close();
}
