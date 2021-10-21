#include "integrator.hpp"
#include <iostream>

void Integrator::calculateEffectiveField(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin)
{
	/*
	Calculates the effective by combining: external field, anisotropy and nearest neighbour interaction
	*/
	for (int i = 0; i < constants::nAtoms; i++)
	{
		// External field
		memcpy(&effectiveField[3 * i], &constants::magneticField, sizeof(constants::magneticField)); // set to external field

		// Anisotropy
		effectiveField[3 * i] -= dotProduct(&constants::anisotropyMatrix[0], &spin[3 * i]);
		effectiveField[3 * i + 1] -= dotProduct(&constants::anisotropyMatrix[3], &spin[3 * i]);
		effectiveField[3 * i + 2] -= dotProduct(&constants::anisotropyMatrix[6], &spin[3 * i]);

		// Nearest neighbours
		for (int j : neighbours[i])
		{
			effectiveField[3 * i] += constants::exchangePrefactor * spin[3 * j];		 // x
			effectiveField[3 * i + 1] += constants::exchangePrefactor * spin[3 * j + 1]; // y
			effectiveField[3 * i + 2] += constants::exchangePrefactor * spin[3 * j + 2]; // z
		}
	}
}

void Integrator::evaluate(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin)
{
	/*
	function evaluation
	*/
	calculateEffectiveField(neighbours, spin);
	for (int i = 0; i < constants::nAtoms; i++)
	{

		rk[3 * i] = constants::gamma * constants::dt * (spin[3 * i + 1] * effectiveField[3 * i + 2] - spin[3 * i + 2] * effectiveField[3 * i + 1] + constants::lambda * (effectiveField[3 * i + 1] * spin[3 * i] * spin[3 * i + 1] - effectiveField[3 * i] * spin[3 * i + 1] * spin[3 * i + 1] + effectiveField[3 * i + 2] * spin[3 * i] * spin[3 * i + 2] - effectiveField[3 * i] * spin[3 * i + 2] * spin[3 * i + 2]));
		rk[3 * i + 1] = constants::gamma * constants::dt * (spin[3 * i + 2] * effectiveField[3 * i] - spin[3 * i] * effectiveField[3 * i + 2] + constants::lambda * (-effectiveField[3 * i + 1] * spin[3 * i] * spin[3 * i] + effectiveField[3 * i] * spin[3 * i] * spin[3 * i + 1] + effectiveField[3 * i + 2] * spin[3 * i + 1] * spin[3 * i + 2] - effectiveField[3 * i + 1] * spin[3 * i + 2] * spin[3 * i + 2]));
		rk[3 * i + 2] = constants::gamma * constants::dt * (spin[3 * i] * effectiveField[3 * i + 1] - spin[3 * i + 1] * effectiveField[3 * i] + constants::lambda * (-effectiveField[3 * i + 2] * spin[3 * i] * spin[3 * i] - effectiveField[3 * i + 2] * spin[3 * i + 1] * spin[3 * i + 1] + effectiveField[3 * i] * spin[3 * i] * spin[3 * i + 2] + effectiveField[3 * i + 1] * spin[3 * i + 1] * spin[3 * i + 2]));
	}
}

void Integrator::integrate(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin, std::vector<double> &randomField)
{
	evaluate(neighbours, spin);
	for (int i = 0; i < constants::nAtoms * 3; i++)
	{
		rkPos[i] = spin[i] + rk[i] / 2;
	}
	evaluate(neighbours, rkPos);

	double temperatureSpin[3];
	for (int i = 0; i < constants::nAtoms; i++)
	{
		temperatureSpin[0] = constants::gamma * (spin[3 * i + 1] * randomField[3 * i + 2] - spin[3 * i + 2] * randomField[3 * i + 1]);
		temperatureSpin[1] = constants::gamma * (spin[3 * i + 2] * randomField[3 * i] - spin[3 * i] * randomField[3 * i + 2]);
		temperatureSpin[2] = constants::gamma * (spin[3 * i] * randomField[3 * i + 1] - spin[3 * i + 1] * randomField[3 * i]);

		spin[3 * i] += rk[3 * i] + temperatureSpin[0];
		spin[3 * i + 1] += rk[3 * i + 1] + temperatureSpin[1];
		spin[3 * i + 2] += rk[3 * i + 2] + temperatureSpin[2];
	}
}

double Integrator::dotProduct(const double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}