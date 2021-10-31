#include <vector>
#include <string>

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
    extern double spinSize;
    extern int baseAtoms;
    extern int nUnitCells;
    extern int nDimensions;
    extern int nAtoms;
    extern std::vector<std::vector<double>> unitVectors;

    // Simulation settings;
    extern  int spinInit;
    extern double angle;
    extern double minDistance;
    extern int nSims;

    // Simulation parameters
    extern double dt;
    extern int steps;
    extern double anisotropy[3];
    extern double anisotropyStrength;
    extern double magneticField[3];
    extern double J;
    extern double lambda;
    extern double temperature;

    // File constants
    extern std::string inputFile;
    extern std::string outputFile;
    extern double offset;
    extern double length;

    // Prefactors
    extern double exchangePrefactor;
    extern double temperatureSigma;
    extern double anisotropyMatrix[9];
}

// Enables printing for debug purposes
// bool printPositions = false;
// bool printNeighbours = false;
// bool printInitialize = false;