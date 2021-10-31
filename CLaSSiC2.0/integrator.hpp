#include "constants.hpp"
#include <vector>

class Integrator
{
private:
	std::vector<double> rk = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> rkPos = std::vector<double>(3 * constants::nAtoms);
	std::vector<double> effectiveField = std::vector<double>(3 * constants::nAtoms);

public:
	void calculateEffectiveField(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin);
	void evaluate(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin);
	void integrate(std::vector<std::vector<double>> &neighbours, std::vector<double> &spin, std::vector<double> &randomField);
	double dotProduct(const double a[3], double b[3]);
};