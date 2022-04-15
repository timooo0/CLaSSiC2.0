#include <vector>
#include <string>

class Vector3D{
public:
	std::vector<float> x;
	std::vector<float> y;
	std::vector<float> z;

    Vector3D (int size, bool pad){
        // size += pad && size%8!=0 ? 8-size%8 : 0;
        x = std::vector<float>(size);
        y = std::vector<float>(size);
        z = std::vector<float>(size);
    }

};

namespace constants
{
    // Physical constants
    extern double boltzmann;
    extern double gFactor;
    extern double bohrMagneton;
    extern double hBar;
    extern double gamma;
    extern double pi;

    // System constants
    extern float spinSize;
    extern int nUnitCells;
    extern int nDimensions;
    extern int nAtoms;
    extern std::vector<std::vector<float>> unitVectors;
    extern std::vector<std::vector<float>> basisPosition;

    // Simulation settings;
    extern int spinInit;
    extern int geometry;
    extern float angle;
    extern int mode;
    extern int nNeighbours;
    extern float minDistance;
    extern bool periodicBoundary;
    extern bool stabilize;

    // Simulation parameters
    extern float dt;
    extern int steps;
    extern int burnInSteps;
    extern float anisotropy[3];
    extern float anisotropyStrength;
    extern float anisotropyAxis;
	extern float anisotropyPlane;
    extern float magneticField[3];
    extern float J;
    extern float lambda;
    extern float temperature;

    // File constants
    extern std::string positionFile;
    extern std::string outputFile;
    extern std::string energyFile;
    extern float offset;
    extern float length;

    // Prefactors
    extern float exchangePrefactor;
    extern float temperatureSigma;
    extern float anisotropyMatrix[9];
}

void setDerivatives();
// Enables printing for debug purposes
// bool printPositions = false;
// bool printNeighbours = false;
// bool printInitialize = false;