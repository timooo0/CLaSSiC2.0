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
		- make sure the data directory exists
		- Calls load to build the crystal
		- Calls initialize to set to start orientation of the spins
	*/
	if (!std::filesystem::is_directory("data") || !std::filesystem::exists("data")) {
    	std::filesystem::create_directory("data"); // create src folder
	}

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

	// Create the basis
	for (int i=0; i < constants::basisPosition.size(); i++){
		for (int j = 0; j < 3; j++){
			position.push_back(constants::basisPosition[i][j]);
		}
	}


	// Add more unit cells
	int offset = 0;
	for (int dim = 0; dim < constants::nDimensions; dim++)
	{
		for (int cell = 1; cell < constants::nUnitCells; cell++)
		{
			if (constants::nDimensions>1){
				if (position.size()%(3*constants::basisPosition.size()*2*constants::nUnitCells)==0 && constants::unitVectors[1][0]!=0){
					offset++;
					std::cout << "Position size: " << position.size()/3 << std::endl;
				}
			}
			// std::cout <<  std::pow(constants::nUnitCells, dim) << std::endl;
			for (int i = 0; i < std::pow(constants::nUnitCells, dim); i++)
			{
				// std::cout <<  constants::basisPosition.size() << std::endl;
				for (int k=0; k<constants::basisPosition.size(); k++){
					for (int j = 0; j < 3; j++)
					{
						// std::cout << "vec: "<<  constants::unitVectors[dim][j] << std::endl;
						// position.push_back(position[3 * i * dim + j + 3*(k*(1+i))] + (constants::unitVectors[dim][j]) * cell-(constants::unitVectors[0][j])*offset);
						// std::cout << position.size()/3 << ": "<<position[3 * k + 3 * i + j] << " | " <<  (constants::unitVectors[dim][j]) * cell << " | " << -(constants::unitVectors[1][j])*offset << std::endl;
						if (j==0){
							position.push_back(position[3 * k + constants::basisPosition.size() * 3 * i + j]+constants::unitVectors[dim][j] * cell-(constants::unitVectors[0][j])*offset);
						} else {
							position.push_back(position[3 * k + constants::basisPosition.size() * 3 * i + j]+constants::unitVectors[dim][j] * cell);
						}
					}
				}
			}
		}
	}
	writePositions();
	std::cout << "positions complete!\n";

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

	std::cout << "nearest neighbours complete!\n";
	// Periodic boundary conditions
	// TODO find a pretty solution
	bool applyBoundry = true;
	if (constants::geometry==3){
		for (int i = 0; i < constants::nAtoms; i++)
		{
			for (int j = 0; j < constants::nAtoms; j++)
			{
				if (i != j)
				{
					applyBoundry = false;
					if ((std::abs(position[3 * j] - position[3 * i])+0.000001 >= (constants::nUnitCells - 1) * constants::unitVectors[0][0] &&
					std::abs(position[3 * j] - position[3 * i])-0.000001 <= (constants::nUnitCells) * constants::unitVectors[0][0])
					&& ( std::abs(position[3 * j + 1]-position[3 * i + 1])-0.000001 <= constants::unitVectors[1][1])
					)
					{
						applyBoundry = true;
					}
					if ((std::abs(position[3 * j + 1] - position[3 * i + 1])+0.000001 >= (constants::nUnitCells - 1) * constants::unitVectors[1][1] &&
					std::abs(position[3 * j]-position[3 * i])-0.000001 <= constants::unitVectors[1][0])
					)
					{
						applyBoundry = true;
					}
					
					if (i==0 && j==constants::nAtoms-1){
						applyBoundry = true;
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
	else if (constants::geometry == 4) {
		for (int i = 0; i < constants::nAtoms; i++)
		{
			for (int j = 0; j < constants::nAtoms; j++)
			{
				if (i != j)
				{
					applyBoundry = false;
					if (i==1 && j==38) {
						std::cout << "pos: " <<position[3 * i+1] << ", "<< position[3 * j+1] <<std::endl;
						std::cout << "nn: " << std::abs(std::abs(position[3 * j] - position[3 * i])-(constants::basisPosition[2][0])) <<
						", " <<std::abs(std::abs(position[3 * j + 1] - position[3 * i + 1])-(constants::nUnitCells * constants::unitVectors[1][1]-constants::basisPosition[2][1])) << std::endl;
					}
					// Horizontal boundary conditions
					if ((std::abs(std::abs(position[3 * j] - position[3 * i])-(constants::nUnitCells * constants::unitVectors[0][0] - constants::basisPosition[2][0])) < 0.000001 &&
					std::abs(std::abs(position[3 * j + 1] - position[3 * i + 1])-constants::basisPosition[2][1]) < 0.000001) ||
					std::abs(std::abs(position[3 * j] - position[3 * i])-(constants::nUnitCells * constants::unitVectors[0][0]- constants::basisPosition[1][0])) < 0.000001 &&
					std::abs(position[3 * j + 1] - position[3 * i + 1]) < 0.000001)
					{
						applyBoundry = true;
						std::cout << "horizontal succes: " << i << ", " << j << std::endl;
					}
					// Vertical boundary conditions
					if (std::abs(std::abs(position[3 * j] - position[3 * i])-(constants::basisPosition[2][0])) < 0.000001 &&
					std::abs(std::abs(position[3 * j + 1] - position[3 * i + 1])-(constants::nUnitCells * constants::unitVectors[1][1]-constants::basisPosition[2][1])) < 0.000001)
					{
						applyBoundry = true;
						// std::cout << "vertical succes: " << i << ", " << j << std::endl;
					}
					// Corner boundary conditions
					if (std::abs(std::abs(position[3 * j] - position[3 * i])-(constants::nUnitCells * constants::unitVectors[0][0] - constants::basisPosition[2][0])) < 0.000001 &&
					std::abs(std::abs(position[3 * j + 1] - position[3 * i + 1])-(constants::nUnitCells * constants::unitVectors[1][1]-constants::basisPosition[2][1])) < 0.000001)
					{
						applyBoundry = true;
						// std::cout << "Corner succes: " << i << ", " << j << std::endl;
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
	else {
		for (int i = 0; i < constants::nAtoms; i++)
		{
			for (int j = 0; j < constants::nAtoms; j++)
			{
				if (i != j)
				{
					applyBoundry = false;
					if (std::abs(position[3 * j] - position[3 * i]) + 0.000001 >= (constants::nUnitCells - 1) * constants::unitVectors[0][0] &&
					std::abs(position[3 * j+1]-position[3 * i+1]) - 0.000001 < constants::unitVectors[0][1])
					{
						applyBoundry = true;
					}
					if (constants::nDimensions > 1){
						if (std::abs(position[3 * j+1]-position[3 * i+1]) + 0.000001 >= (constants::nUnitCells - 1) * constants::unitVectors[1][1] &&
						std::abs(position[3 * j]-position[3 * i]) - 0.000001 <= constants::unitVectors[1][0])
						{
							applyBoundry = true;
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
	std::cout << "periodic boundaries complete!\n";

	// Printing

	std::cout << "Atom positions: \n";
	for (int i = 0; i < position.size() / 3; i++)
	{
		std::cout << "{" << position[3 * i] << " , " << position[3 * i + 1] << "}, ";
	}


	std::cout << "Neighbours: \n";
	for (int i = 0; i < constants::nAtoms; i++)
	{
		std::cout << i << "(" << neighbours[i].size() << ")"<<  ": ";
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
		spin[3] = -std::cos(constants::angle);
		spin[4] = -std::sin(constants::angle); 
		spin[5] = std::cos(constants::angle);
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
	break;
	case 5:
		// Triangle
		std::cout << constants::spinInit;	
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
		case 6:
		// Small z angle for spin waves
		spin[0] = 1;
		spin[1] = 0;
		spin[2] = 0.2;
		for (int i = 1; i < (int)constants::nAtoms; i++)
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
				spin[3 * i] = 0;
				spin[3 * i + 1] = 0;
				spin[3 * i + 2] = swap * 1;
		}
		normalize();
	break;
		case 7:
		for (int i = 0; i < constants::nAtoms; i++){
			spin[3*i] = std::sin(2*constants::pi*i/3);
			spin[3*i+1] = std::cos(2*constants::pi*i/3);
			spin[3*i+2] = 0;
		}
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
	std::ofstream file;
	file.open(addFileNumber(constants::outputFile), std::ios::binary);
	writeConstants(file);

	// Temperature fluctatuations
	std::random_device seed;
	std::default_random_engine engine{seed()};
	std::normal_distribution<double> temperatureDistribution;
	if (constants::temperatureSigma != 0 && constants::temperature != 0)
	{
		temperatureDistribution = std::normal_distribution<double>(0, constants::temperatureSigma);
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
std::string Simulation::addFileNumber(std::string input){
	int i=0;

	std::string fileString = input;
	fileString.insert(input.find("."), std::to_string(i));
	while (std::filesystem::exists(fileString))
	{
		i++;
		fileString = constants::outputFile;
		fileString.insert(input.find("."), std::to_string(i));
	}

	return fileString;
}


void Simulation::writeConstants(std::ofstream &f)
{
	/*
	Writes all constants and system parameters to the start of the file
	*/
	double saveAtoms = constants::nAtoms;
	double saveSteps = constants::steps;
	double saveGeometry = constants::geometry;
	double saveNUnitcells = constants::nUnitCells;
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
	f.write((char *)&saveGeometry, sizeof(double));
	f.write((char *)&saveNUnitcells, sizeof(double));
	f.write((char *)&constants::anisotropyStrength, sizeof(double));
}

void Simulation::writePositions()
{

	std::ofstream positionFile;
	positionFile.open(addFileNumber(constants::positionFile));

	std::string seperator;
	for (int i=0; i < position.size();i++){
		positionFile << position[i];
		seperator = ((i+1)%3!=0) ? ", " : "\n";
		positionFile << seperator;
	}
	positionFile << "\n";
	positionFile.close();
}


std::vector<double>* Simulation::getSpin(){
	return &spin;
}
std::vector<double>* Simulation::getPosition(){
	return &position;
}
std::vector<double>* Simulation::getRandomField(){
	return &randomField;
}
std::vector<std::vector<int>>* Simulation::getNeighbours(){
	return &neighbours;
}