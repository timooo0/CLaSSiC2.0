#include "constants.hpp"
#include <vector>

class Integrator
{
private:
	Vector3D rk = Vector3D(constants::nAtoms, true);
	Vector3D rkPos = Vector3D(constants::nAtoms, true);
	Vector3D effectiveField = Vector3D(constants::nAtoms, true);
	Vector3D stabilizerField = Vector3D(constants::nAtoms, true);

public:
	void calculateEffectiveField(std::vector<std::vector<int>> &neighbours, Vector3D &spin);
	void evaluate(std::vector<std::vector<int>> &neighbours, Vector3D &spin);
	void integrate(std::vector<std::vector<int>> &neighbours, Vector3D &spin, std::vector<float> &randomField);
	float calculateEnergy(Vector3D &spin);
	float dotProduct(const float a[3], float b[3]);

	Vector3D* getEffectiveField();
	void setStabilizerField();
};