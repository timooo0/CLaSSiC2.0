#include "constants.cpp"
#include <vector>

class Integrator {
private:
	double rk[3*constants::nAtoms];
	double rkPos[3 * constants::nAtoms];
	double effectiveField[3 * constants::nAtoms];
public:
	void calculateEffectiveField(std::vector<std::vector<double>>& neighbours, double spin[constants::nAtoms * 3]);
	void rungeKutta(std::vector<std::vector<double>>& neighbours, double spin[constants::nAtoms*3]);
	void integrate(std::vector<std::vector<double>>& neighbours, double spin[constants::nAtoms * 3], double randomField[constants::nAtoms * 3]);
	double dotProduct(double a[3], double b[3]);
	
};