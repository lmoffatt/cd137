#ifndef RESULTSSIMULATOR_H
#define RESULTSSIMULATOR_H
#include "Results.h"
#include "SimParameters.h"
#include "Cell_simulator.h"

class ResultsSimulator {
public:
    Results Simulate(const SimParameters& simPar,
                     const Results& results);

    double SumSquare(const SimParameters& SimPar,
                     const Results& results);



};


#endif // RESULTSSIMULATOR_H
