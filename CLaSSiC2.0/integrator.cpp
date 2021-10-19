#include "integrator.hpp"
#include <iostream>

void Integrator::calculateEffectiveField(std::vector<std::vector<double>>& neighbours, std::vector<double>& spin)
{
	for (int i = 0; i < constants::nAtoms; i++) {
		memcpy(&effectiveField[3 * i], &constants::magneticField, sizeof(constants::magneticField)); //set to external field
		effectiveField[3 * i]     -= dotProduct(&constants::anisotropyMatrix[0], &spin[3 * i]);
		effectiveField[3 * i + 1] -= dotProduct(&constants::anisotropyMatrix[3], &spin[3 * i]);
		effectiveField[3 * i + 2] -= dotProduct(&constants::anisotropyMatrix[6], &spin[3 * i]);

		for (int j : neighbours[i]) {
			effectiveField[3 * i]     += constants::exchangePrefactor * spin[3 * j];		//x
			effectiveField[3 * i + 1] += constants::exchangePrefactor * spin[3 * j + 1];	//y
			effectiveField[3 * i + 2] += constants::exchangePrefactor * spin[3 * j + 2];	//z
		}
	}
}

void Integrator::rungeKutta(std::vector<std::vector<double>>& neighbours, std::vector<double>& spin)
{
	double pDotp;
	double pDotB;

	calculateEffectiveField(neighbours, spin);
	for (int i = 0; i < constants::nAtoms; i++) {
		pDotp = dotProduct(&spin[3 * i], &spin[3 * i]);
		pDotB = dotProduct(&spin[3 * i], &effectiveField[3 * i]);

		rk[3 * i] = constants::gamma * constants::dt * (spin[3 * i + 1] * effectiveField[3 * i + 2] - spin[3 * i + 2] * effectiveField[3 * i + 1]
			// + constants::lambda*(effectiveField[3*i]*pDotp - spin[3*i]*pDotB));
			+ constants::lambda*(
				effectiveField[3*i+1] * spin[3*i] * spin[3*i+1]
				- effectiveField[3*i] * spin[3*i+1] * spin[3*i+1]
				+ effectiveField[3*i+2] * spin[3*i] * spin[3*i+2]
				- effectiveField[3*i] * spin[3*i+2] * spin[3*i+2]));
		rk[3*i + 1] = constants::gamma * constants::dt * (spin[3 * i + 2] * effectiveField[3 * i] - spin[3 * i] * effectiveField[3 * i + 2]
			// + constants::lambda * (effectiveField[3 * i + 1] * pDotp - spin[3 * i + 1] * pDotB));
			+constants::lambda*(
				- effectiveField[3*i+1] * spin[3*i]*spin[3*i]
				+ effectiveField[3*i] * spin[3*i] * spin[3*i+1]
				+ effectiveField[3*i+2] * spin[3*i+1] * spin[3*i+2]
				- effectiveField[3*i+1] * spin[3*i+2]*spin[3*i+2]));
		rk[3 * i + 2] =  constants::gamma * constants::dt *(spin[3 * i] * effectiveField[3 * i + 1] - spin[3 * i + 1] * effectiveField[3 * i]
			// + constants::lambda * (effectiveField[3 * i + 2] * pDotp - spin[3 * i + 2] * pDotB));
			+constants::lambda*(
				- effectiveField[3*i+2] * spin[3*i]*spin[3*i]
				- effectiveField[3*i+2] * spin[3*i+1]*spin[3*i+1]
				+ effectiveField[3*i] * spin[3*i]*spin[3*i+2]
				+ effectiveField[3*i+1] * spin[3*i+1]*spin[3*i+2]));

		// if (rk[3*i+2]!=0){
		// 	std::cout.precision(std::numeric_limits<double>::digits10+2);
		// 	std::cout << "Unequal to zero: \n";
		// 	std::cout << spin[3 * i] * effectiveField[3 * i + 1] << std::endl;
		// 	std::cout << spin[3 * i + 1] * effectiveField[3 * i] << std::endl;
		// 	std::cout << spin[3 * i] * effectiveField[3 * i + 1] - spin[3 * i + 1] * effectiveField[3 * i] << std::endl;
		// }
	}
}

void Integrator::integrate(std::vector<std::vector<double>>& neighbours, std::vector<double>& spin, std::vector<double>& randomField)
{
	rungeKutta(neighbours, spin);
	rkSave = rk;
	for (int i = 0; i < constants::nAtoms * 3; i++) {
		rkPos[i] = spin[i] + rk[i]/2;
	}
	rungeKutta(neighbours, rkPos);

	double temperatureSpin[3];
	for (int i = 0; i < constants::nAtoms; i++) {
		temperatureSpin[0] = constants::gamma*(spin[3*i+1]*randomField[3*i+2]-spin[3*i+2]*randomField[3*i+1]);
		temperatureSpin[1] = constants::gamma * (spin[3 * i + 2] * randomField[3 * i] - spin[3 * i] * randomField[3 * i + 2]);
		temperatureSpin[2] = constants::gamma * (spin[3 * i] * randomField[3 * i + 1] - spin[3 * i + 1] * randomField[3 * i]);

		// spin[3*i] += 0.5 * rk[3*i] * constants::dt + 0.5 * rkSave[3*i + 0]*constants::dt+ temperatureSpin[0];
		// spin[3 * i+1] +=  0.5 * rk[3*i+1] * constants::dt + 0.5 * rkSave[3*i + 1]*constants::dt+ temperatureSpin[1];
		// spin[3 * i+2] +=  0.5 * rk[3*i+2] * constants::dt + 0.5 * rkSave[3*i + 2]*constants::dt+ temperatureSpin[2];
		spin[3*i] 	  += rk[3*i]   + temperatureSpin[0];
		spin[3 * i+1] += rk[3*i+1] + temperatureSpin[1];
		spin[3 * i+2] += rk[3*i+2] + temperatureSpin[2];
	}
}

double Integrator::dotProduct(double a[3], double b[3])
{
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double Integrator::dotProduct(const double a[3], double b[3])
{
	return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}